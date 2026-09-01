/**
 * \file cuda_helpers.cu
 * \brief Standalone CUDA implementations of SAXS operations.
 */

#include "cuda_helpers.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace foxs_cuda {
namespace saxs {
namespace internal {
namespace {

constexpr std::size_t MAX_THREADS = 512;

void check_cuda(cudaError_t status, const char* operation) {
  if (status == cudaSuccess) return;
  std::ostringstream message;
  message << operation << ": " << cudaGetErrorString(status);
  throw std::runtime_error(message.str());
}

template <class T>
class DeviceBuffer {
 public:
  explicit DeviceBuffer(std::size_t size) : pointer_(nullptr), size_(size) {
    if (size_ > 0) {
      check_cuda(cudaMalloc(reinterpret_cast<void**>(&pointer_),
                            size_ * sizeof(T)),
                 "cudaMalloc");
    }
  }

  ~DeviceBuffer() {
    if (pointer_) cudaFree(pointer_);
  }

  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;

  T* get() { return pointer_; }
  const T* get() const { return pointer_; }

  void copy_from_host(const T* source) {
    if (size_ > 0) {
      check_cuda(cudaMemcpy(pointer_, source, size_ * sizeof(T),
                            cudaMemcpyHostToDevice),
                 "cudaMemcpy host to device");
    }
  }

  void copy_to_host(T* destination) const {
    if (size_ > 0) {
      check_cuda(cudaMemcpy(destination, pointer_, size_ * sizeof(T),
                            cudaMemcpyDeviceToHost),
                 "cudaMemcpy device to host");
    }
  }

 private:
  T* pointer_;
  std::size_t size_;
};

std::vector<double> flatten(
    const std::vector<std::vector<double>>& channels,
    std::size_t particle_count) {
  std::vector<double> values;
  values.reserve(channels.size() * particle_count);
  for (const std::vector<double>& channel : channels) {
    if (channel.size() != particle_count) {
      throw std::invalid_argument("CUDA form-factor channel has wrong size");
    }
    values.insert(values.end(), channel.begin(), channel.end());
  }
  return values;
}

std::size_t output_channel_count(std::size_t factor_channel_count) {
  if (factor_channel_count == 1) return 1;
  if (factor_channel_count == 2) return 3;
  if (factor_channel_count == 3) return 6;
  throw std::invalid_argument(
      "CUDA distance distribution requires one, two, or three form-factor channels");
}

double maximum_squared_distance(const std::vector<double>& coordinates1,
                                const std::vector<double>& coordinates2,
                                bool same_particles) {
  const std::vector<double>& coordinates_b =
      same_particles ? coordinates1 : coordinates2;
  double bound = 0.0;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    double min1 = coordinates1[axis], max1 = coordinates1[axis];
    double min2 = coordinates_b[axis], max2 = coordinates_b[axis];
    for (std::size_t i = axis; i < coordinates1.size(); i += 3) {
      min1 = std::min(min1, coordinates1[i]);
      max1 = std::max(max1, coordinates1[i]);
    }
    for (std::size_t i = axis; i < coordinates_b.size(); i += 3) {
      min2 = std::min(min2, coordinates_b[i]);
      max2 = std::max(max2, coordinates_b[i]);
    }
    const double separation =
        std::max(std::abs(max1 - min2), std::abs(max2 - min1));
    bound += separation * separation;
  }
  return bound;
}

__device__ double atomic_add_double(double* address, double value) {
#if __CUDA_ARCH__ >= 600
  return atomicAdd(address, value);
#else
  unsigned long long int* address_as_ull =
      reinterpret_cast<unsigned long long int*>(address);
  unsigned long long int old = *address_as_ull;
  unsigned long long int assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(value + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
#endif
}

__device__ void add_pair_weights(double* output, std::size_t bin,
                                 std::size_t bin_count,
                                 const double* factors1,
                                 const double* factors2,
                                 std::size_t i, std::size_t j,
                                 std::size_t count1, std::size_t count2,
                                 std::size_t factor_channels,
                                 bool diagonal) {
  const double a0 = factors1[i];
  const double b0 = factors2[j];
  const double pair_scale = diagonal ? 1.0 : 2.0;
  atomic_add_double(output + bin, pair_scale * a0 * b0);
  if (factor_channels < 2) return;

  const double a1 = factors1[count1 + i];
  const double b1 = factors2[count2 + j];
  atomic_add_double(output + bin_count + bin, pair_scale * a1 * b1);
  atomic_add_double(output + 2 * bin_count + bin,
                    diagonal ? 2.0 * a0 * a1
                             : 2.0 * (a0 * b1 + b0 * a1));
  if (factor_channels < 3) return;

  const double a2 = factors1[2 * count1 + i];
  const double b2 = factors2[2 * count2 + j];
  atomic_add_double(output + 3 * bin_count + bin, pair_scale * a2 * b2);
  atomic_add_double(output + 4 * bin_count + bin,
                    diagonal ? 2.0 * a0 * a2
                             : 2.0 * (a0 * b2 + b0 * a2));
  atomic_add_double(output + 5 * bin_count + bin,
                    diagonal ? 2.0 * a2 * a1
                             : 2.0 * (a2 * b1 + b2 * a1));
}

__global__ void make_distance_distributions(
    const double* coordinates1, const double* factors1, std::size_t count1,
    const double* coordinates2, const double* factors2, std::size_t count2,
    bool same_particles, std::size_t factor_channels, double inverse_bin_size,
    std::size_t bin_count, double* output) {
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;
  if (i >= count1 || j >= count2) return;
  if (same_particles && j < i) return;

  const double dx = coordinates1[3 * i] - coordinates2[3 * j];
  const double dy = coordinates1[3 * i + 1] - coordinates2[3 * j + 1];
  const double dz = coordinates1[3 * i + 2] - coordinates2[3 * j + 2];
  const double squared_distance = dx * dx + dy * dy + dz * dz;
  std::size_t bin = static_cast<std::size_t>(
      floor(squared_distance * inverse_bin_size + 0.5));
  if (bin >= bin_count) bin = bin_count - 1;

  add_pair_weights(output, bin, bin_count, factors1, factors2, i, j,
                   count1, count2, factor_channels,
                   same_particles && i == j);
}

template <class T>
__device__ T square(T value) {
  return value * value;
}

__device__ double sinc_pi(double value) {
  return fabs(value) < 1e-6 ? 1.0 : sin(value) / value;
}

__global__ void make_profile(const double* radial_distribution, const double* q,
                             const double* distances, double* intensity,
                             double modulation_function_parameter,
                             std::size_t radial_size, std::size_t q_size) {
  __shared__ double partial[MAX_THREADS];
  const std::size_t k = blockIdx.x;
  if (k >= q_size) return;
  partial[threadIdx.x] = 0.0;
  for (std::size_t r = threadIdx.x; r < radial_size; r += blockDim.x) {
    partial[threadIdx.x] +=
        radial_distribution[r] * sinc_pi(distances[r] * q[k]);
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    double total = 0.0;
    for (std::size_t i = 0; i < blockDim.x; ++i) total += partial[i];
    intensity[k] =
        total * exp(-modulation_function_parameter * square(q[k]));
  }
}

}  // namespace

bool cuda_device_available(std::string* reason) {
  int count = 0;
  const cudaError_t status = cudaGetDeviceCount(&count);
  if (status != cudaSuccess) {
    if (reason) *reason = cudaGetErrorString(status);
    return false;
  }
  if (count == 0) {
    if (reason) *reason = "no CUDA devices found";
    return false;
  }
  if (reason) reason->clear();
  return true;
}

void distance_distributions_cuda(
    const std::vector<double>& coordinates1,
    const std::vector<std::vector<double>>& form_factors1,
    const std::vector<double>& coordinates2,
    const std::vector<std::vector<double>>& form_factors2,
    bool same_particles, double bin_size,
    std::vector<std::vector<double>>& distributions) {
  if (coordinates1.size() % 3 != 0 || coordinates2.size() % 3 != 0) {
    throw std::invalid_argument("CUDA coordinates must be packed xyz triples");
  }
  const std::size_t count1 = coordinates1.size() / 3;
  const std::size_t count2 = coordinates2.size() / 3;
  if (same_particles && count1 != count2) {
    throw std::invalid_argument("CUDA same-particle inputs have different sizes");
  }
  if (form_factors1.size() != form_factors2.size()) {
    throw std::invalid_argument("CUDA inputs have different channel counts");
  }
  if (!(bin_size > 0.0)) {
    throw std::invalid_argument("CUDA distribution bin size must be positive");
  }

  const std::size_t factor_channels = form_factors1.size();
  const std::size_t output_channels = output_channel_count(factor_channels);
  if (count1 == 0 || count2 == 0) {
    distributions.assign(output_channels, std::vector<double>(1, 0.0));
    return;
  }
  const std::vector<double> flat_factors1 = flatten(form_factors1, count1);
  const std::vector<double> flat_factors2 = flatten(form_factors2, count2);
  const double max_distance =
      maximum_squared_distance(coordinates1, coordinates2, same_particles);
  const std::size_t bin_count =
      static_cast<std::size_t>(std::floor(max_distance / bin_size + 0.5)) + 2;

  DeviceBuffer<double> device_coordinates1(coordinates1.size());
  DeviceBuffer<double> device_factors1(flat_factors1.size());
  DeviceBuffer<double> device_coordinates2(coordinates2.size());
  DeviceBuffer<double> device_factors2(flat_factors2.size());
  DeviceBuffer<double> device_output(output_channels * bin_count);
  device_coordinates1.copy_from_host(coordinates1.data());
  device_factors1.copy_from_host(flat_factors1.data());
  device_coordinates2.copy_from_host(coordinates2.data());
  device_factors2.copy_from_host(flat_factors2.data());
  check_cuda(cudaMemset(device_output.get(), 0,
                        output_channels * bin_count * sizeof(double)),
             "cudaMemset");

  const dim3 threads(16, 16);
  const dim3 blocks((count1 + threads.x - 1) / threads.x,
                    (count2 + threads.y - 1) / threads.y);
  make_distance_distributions<<<blocks, threads>>>(
      device_coordinates1.get(), device_factors1.get(), count1,
      device_coordinates2.get(), device_factors2.get(), count2,
      same_particles, factor_channels, 1.0 / bin_size, bin_count,
      device_output.get());
  check_cuda(cudaGetLastError(), "distance-distribution kernel launch");
  check_cuda(cudaDeviceSynchronize(), "distance-distribution kernel");

  std::vector<double> flat_output(output_channels * bin_count);
  device_output.copy_to_host(flat_output.data());
  distributions.assign(output_channels, std::vector<double>(bin_count));
  for (std::size_t channel = 0; channel < output_channels; ++channel) {
    std::copy(flat_output.begin() + channel * bin_count,
              flat_output.begin() + (channel + 1) * bin_count,
              distributions[channel].begin());
  }
}

void squared_distribution_2_profile_cuda(
    const double* radial_distribution, const double* q,
    const double* distances, double* intensity,
    double modulation_function_parameter, std::size_t radial_size,
    std::size_t q_size) {
  DeviceBuffer<double> device_distribution(radial_size);
  DeviceBuffer<double> device_distances(radial_size);
  DeviceBuffer<double> device_q(q_size);
  DeviceBuffer<double> device_intensity(q_size);
  device_distribution.copy_from_host(radial_distribution);
  device_distances.copy_from_host(distances);
  device_q.copy_from_host(q);

  const std::size_t thread_count = std::min(MAX_THREADS, radial_size);
  make_profile<<<q_size, thread_count>>>(
      device_distribution.get(), device_q.get(), device_distances.get(),
      device_intensity.get(), modulation_function_parameter, radial_size,
      q_size);
  check_cuda(cudaGetLastError(), "profile kernel launch");
  check_cuda(cudaDeviceSynchronize(), "profile kernel");
  device_intensity.copy_to_host(intensity);
}

}  // namespace internal
}  // namespace saxs
}  // namespace foxs_cuda
