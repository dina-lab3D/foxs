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

constexpr unsigned int PAIR_TILE = 32;
constexpr unsigned int PAIR_THREADS = 256;
constexpr unsigned int PAIR_ROWS_PER_ITERATION = PAIR_THREADS / PAIR_TILE;
constexpr unsigned int PAIRS_PER_THREAD =
    PAIR_TILE * PAIR_TILE / PAIR_THREADS;
constexpr unsigned int PROFILE_THREADS = 256;
constexpr unsigned int MAX_PROFILE_CHANNELS = 6;

static_assert(PAIR_THREADS % PAIR_TILE == 0,
              "Pair thread count must be divisible by the tile width");
static_assert(PAIR_TILE % PAIR_ROWS_PER_ITERATION == 0,
              "Pair threads must cover a whole number of tile rows");
static_assert(PAIRS_PER_THREAD == 4,
              "The 32 by 32 pair tile must assign four pairs per thread");

void check_cuda(cudaError_t status, const char* operation) {
  if (status == cudaSuccess) return;
  std::ostringstream message;
  message << operation << ": " << cudaGetErrorString(status);
  throw std::runtime_error(message.str());
}

template <class T>
class ReusableDeviceBuffer {
 public:
  ReusableDeviceBuffer() : pointer_(nullptr), capacity_(0) {}
  ~ReusableDeviceBuffer() {
    if (pointer_) cudaFree(pointer_);
  }

  ReusableDeviceBuffer(const ReusableDeviceBuffer&) = delete;
  ReusableDeviceBuffer& operator=(const ReusableDeviceBuffer&) = delete;

  void reserve(std::size_t size) {
    if (size <= capacity_) return;
    if (pointer_) {
      check_cuda(cudaFree(pointer_), "cudaFree while growing CUDA workspace");
      pointer_ = nullptr;
      capacity_ = 0;
    }
    check_cuda(cudaMalloc(reinterpret_cast<void**>(&pointer_),
                          size * sizeof(T)),
               "cudaMalloc for CUDA workspace");
    capacity_ = size;
  }

  T* get() { return pointer_; }
  const T* get() const { return pointer_; }

  void copy_from_host(const T* source, std::size_t size) {
    reserve(size);
    if (size > 0) {
      check_cuda(cudaMemcpy(pointer_, source, size * sizeof(T),
                            cudaMemcpyHostToDevice),
                 "cudaMemcpy host to device");
    }
  }

  void copy_to_host(T* destination, std::size_t size) const {
    if (size > capacity_) {
      throw std::logic_error("CUDA workspace copy exceeds buffer capacity");
    }
    if (size > 0) {
      check_cuda(cudaMemcpy(destination, pointer_, size * sizeof(T),
                            cudaMemcpyDeviceToHost),
                 "cudaMemcpy device to host");
    }
  }

 private:
  T* pointer_;
  std::size_t capacity_;
};

struct CudaWorkspace {
  ReusableDeviceBuffer<double> coordinates1;
  ReusableDeviceBuffer<double> coordinates2;
  ReusableDeviceBuffer<double> factors1;
  ReusableDeviceBuffer<double> factors2;
  ReusableDeviceBuffer<double> distributions;
  ReusableDeviceBuffer<double> q;
  ReusableDeviceBuffer<double> profiles;
  ReusableDeviceBuffer<double> distances;
};

thread_local CudaWorkspace workspace;

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

struct DistributionShape {
  std::size_t count1;
  std::size_t count2;
  std::size_t factor_channels;
  std::size_t output_channels;
  std::size_t bin_count;
};

DistributionShape validate_and_get_shape(
    const std::vector<double>& coordinates1,
    const std::vector<std::vector<double>>& form_factors1,
    const std::vector<double>& coordinates2,
    const std::vector<std::vector<double>>& form_factors2,
    bool same_particles, double bin_size) {
  if (coordinates1.size() % 3 != 0 || coordinates2.size() % 3 != 0) {
    throw std::invalid_argument("CUDA coordinates must be packed xyz triples");
  }
  if (!(bin_size > 0.0)) {
    throw std::invalid_argument("CUDA distribution bin size must be positive");
  }
  DistributionShape shape;
  shape.count1 = coordinates1.size() / 3;
  shape.count2 = coordinates2.size() / 3;
  if (same_particles && shape.count1 != shape.count2) {
    throw std::invalid_argument("CUDA same-particle inputs have different sizes");
  }
  if (form_factors1.size() != form_factors2.size()) {
    throw std::invalid_argument("CUDA inputs have different channel counts");
  }
  shape.factor_channels = form_factors1.size();
  shape.output_channels = output_channel_count(shape.factor_channels);
  if (shape.count1 == 0 || shape.count2 == 0) {
    shape.bin_count = 1;
  } else {
    const double maximum =
        maximum_squared_distance(coordinates1, coordinates2, same_particles);
    shape.bin_count =
        static_cast<std::size_t>(std::floor(maximum / bin_size + 0.5)) + 2;
  }
  return shape;
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

__device__ void load_tile(const double* coordinates, const double* factors,
                          std::size_t count, std::size_t base,
                          std::size_t factor_channels,
                          double* tile_coordinates, double* tile_factors,
                          unsigned int lane) {
  if (lane >= PAIR_TILE) return;
  const std::size_t index = base + lane;
  if (index < count) {
    tile_coordinates[3 * lane] = coordinates[3 * index];
    tile_coordinates[3 * lane + 1] = coordinates[3 * index + 1];
    tile_coordinates[3 * lane + 2] = coordinates[3 * index + 2];
    for (std::size_t channel = 0; channel < factor_channels; ++channel) {
      tile_factors[channel * PAIR_TILE + lane] =
          factors[channel * count + index];
    }
  }
}

__device__ void process_pair(const double* tile_coordinates1,
                             const double* tile_factors1,
                             const double* tile_coordinates2,
                             const double* tile_factors2,
                             unsigned int local1, unsigned int local2,
                             std::size_t factor_channels,
                             double inverse_bin_size, std::size_t bin_count,
                             bool diagonal, double* output) {
  const double dx = tile_coordinates1[3 * local1] -
                    tile_coordinates2[3 * local2];
  const double dy = tile_coordinates1[3 * local1 + 1] -
                    tile_coordinates2[3 * local2 + 1];
  const double dz = tile_coordinates1[3 * local1 + 2] -
                    tile_coordinates2[3 * local2 + 2];
  const double squared_distance = dx * dx + dy * dy + dz * dz;
  std::size_t bin = static_cast<std::size_t>(
      floor(squared_distance * inverse_bin_size + 0.5));
  if (bin >= bin_count) bin = bin_count - 1;

  // Tile factors are packed using PAIR_TILE rather than the global counts.
  add_pair_weights(output, bin, bin_count, tile_factors1, tile_factors2,
                   local1, local2, PAIR_TILE, PAIR_TILE, factor_channels,
                   diagonal);
}

__global__ void make_same_distance_distributions_tiled(
    const double* coordinates, const double* factors, std::size_t count,
    std::size_t factor_channels, double inverse_bin_size,
    std::size_t bin_count, double* output) {
  __shared__ double coordinates1[3 * PAIR_TILE];
  __shared__ double coordinates2[3 * PAIR_TILE];
  __shared__ double factors1[3 * PAIR_TILE];
  __shared__ double factors2[3 * PAIR_TILE];

  const unsigned long long block = blockIdx.x;
  unsigned long long tile2 = static_cast<unsigned long long>(
      floor((sqrt(8.0 * static_cast<double>(block) + 1.0) - 1.0) * 0.5));
  while ((tile2 + 1) * (tile2 + 2) / 2 <= block) ++tile2;
  while (tile2 * (tile2 + 1) / 2 > block) --tile2;
  const unsigned long long tile1 = block - tile2 * (tile2 + 1) / 2;

  const unsigned int thread = threadIdx.x;
  load_tile(coordinates, factors, count, tile1 * PAIR_TILE, factor_channels,
            coordinates1, factors1, thread);
  load_tile(coordinates, factors, count, tile2 * PAIR_TILE, factor_channels,
            coordinates2, factors2, thread);
  __syncthreads();

  const unsigned int local2 = thread % PAIR_TILE;
  const std::size_t global2 = tile2 * PAIR_TILE + local2;
  if (global2 >= count) return;

#pragma unroll
  for (unsigned int iteration = 0; iteration < PAIRS_PER_THREAD;
       ++iteration) {
    const unsigned int local1 =
        thread / PAIR_TILE + iteration * PAIR_ROWS_PER_ITERATION;
    const std::size_t global1 = tile1 * PAIR_TILE + local1;
    if (global1 >= count) break;
    if (tile1 == tile2 && local2 < local1) continue;
    process_pair(coordinates1, factors1, coordinates2, factors2,
                 local1, local2, factor_channels, inverse_bin_size,
                 bin_count, tile1 == tile2 && local1 == local2, output);
  }
}

__global__ void make_cross_distance_distributions_tiled(
    const double* input_coordinates1, const double* input_factors1,
    std::size_t count1, const double* input_coordinates2,
    const double* input_factors2, std::size_t count2,
    std::size_t factor_channels, double inverse_bin_size,
    std::size_t bin_count, double* output) {
  __shared__ double coordinates1[3 * PAIR_TILE];
  __shared__ double coordinates2[3 * PAIR_TILE];
  __shared__ double factors1[3 * PAIR_TILE];
  __shared__ double factors2[3 * PAIR_TILE];

  const unsigned int thread = threadIdx.x;
  const std::size_t base1 = blockIdx.x * PAIR_TILE;
  const std::size_t base2 = blockIdx.y * PAIR_TILE;
  load_tile(input_coordinates1, input_factors1, count1, base1,
            factor_channels, coordinates1, factors1, thread);
  load_tile(input_coordinates2, input_factors2, count2, base2,
            factor_channels, coordinates2, factors2, thread);
  __syncthreads();

  const unsigned int local2 = thread % PAIR_TILE;
  const std::size_t global2 = base2 + local2;
  if (global2 >= count2) return;

#pragma unroll
  for (unsigned int iteration = 0; iteration < PAIRS_PER_THREAD;
       ++iteration) {
    const unsigned int local1 =
        thread / PAIR_TILE + iteration * PAIR_ROWS_PER_ITERATION;
    const std::size_t global1 = base1 + local1;
    if (global1 >= count1) break;
    process_pair(coordinates1, factors1, coordinates2, factors2,
                 local1, local2, factor_channels, inverse_bin_size,
                 bin_count, false, output);
  }
}

__device__ double sinc(double value) {
  return fabs(value) < 1e-8 ? 1.0 : sin(value) / value;
}

__global__ void transform_distance_distributions(
    const double* distributions, std::size_t bin_count,
    std::size_t channel_count, double bin_size, const double* q,
    std::size_t q_count, double modulation_function_parameter,
    double* profiles) {
  __shared__ double partial[MAX_PROFILE_CHANNELS][PROFILE_THREADS];
  const std::size_t q_index = blockIdx.x;
  if (q_index >= q_count) return;

  double sums[MAX_PROFILE_CHANNELS] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (std::size_t bin = threadIdx.x; bin < bin_count;
       bin += blockDim.x) {
    const double kernel = sinc(sqrt(bin * bin_size) * q[q_index]);
    for (std::size_t channel = 0; channel < channel_count; ++channel) {
      sums[channel] += distributions[channel * bin_count + bin] * kernel;
    }
  }
  for (std::size_t channel = 0; channel < channel_count; ++channel) {
    partial[channel][threadIdx.x] = sums[channel];
  }
  __syncthreads();

  for (unsigned int stride = blockDim.x / 2; stride >= 32; stride >>= 1) {
    if (threadIdx.x < stride) {
      for (std::size_t channel = 0; channel < channel_count; ++channel) {
        partial[channel][threadIdx.x] +=
            partial[channel][threadIdx.x + stride];
      }
    }
    __syncthreads();
  }

  if (threadIdx.x < 32) {
    for (std::size_t channel = 0; channel < channel_count; ++channel) {
      double value = partial[channel][threadIdx.x];
      value += __shfl_down_sync(0xffffffff, value, 16);
      value += __shfl_down_sync(0xffffffff, value, 8);
      value += __shfl_down_sync(0xffffffff, value, 4);
      value += __shfl_down_sync(0xffffffff, value, 2);
      value += __shfl_down_sync(0xffffffff, value, 1);
      if (threadIdx.x == 0) {
        profiles[channel * q_count + q_index] =
            value * exp(-modulation_function_parameter * q[q_index] *
                        q[q_index]);
      }
    }
  }
}

__global__ void transform_precomputed_distances(
    const double* distribution, const double* distances,
    std::size_t bin_count, const double* q, std::size_t q_count,
    double modulation_function_parameter, double* profile) {
  __shared__ double partial[PROFILE_THREADS];
  const std::size_t q_index = blockIdx.x;
  if (q_index >= q_count) return;
  double sum = 0.0;
  for (std::size_t bin = threadIdx.x; bin < bin_count;
       bin += blockDim.x) {
    sum += distribution[bin] * sinc(distances[bin] * q[q_index]);
  }
  partial[threadIdx.x] = sum;
  __syncthreads();
  for (unsigned int stride = blockDim.x / 2; stride >= 32; stride >>= 1) {
    if (threadIdx.x < stride) partial[threadIdx.x] += partial[threadIdx.x + stride];
    __syncthreads();
  }
  if (threadIdx.x < 32) {
    double value = partial[threadIdx.x];
    value += __shfl_down_sync(0xffffffff, value, 16);
    value += __shfl_down_sync(0xffffffff, value, 8);
    value += __shfl_down_sync(0xffffffff, value, 4);
    value += __shfl_down_sync(0xffffffff, value, 2);
    value += __shfl_down_sync(0xffffffff, value, 1);
    if (threadIdx.x == 0) {
      profile[q_index] =
          value * exp(-modulation_function_parameter * q[q_index] * q[q_index]);
    }
  }
}

DistributionShape calculate_distributions_on_device(
    const std::vector<double>& coordinates1,
    const std::vector<std::vector<double>>& form_factors1,
    const std::vector<double>& coordinates2,
    const std::vector<std::vector<double>>& form_factors2,
    bool same_particles, double bin_size) {
  const DistributionShape shape = validate_and_get_shape(
      coordinates1, form_factors1, coordinates2, form_factors2,
      same_particles, bin_size);
  const std::vector<double> flat_factors1 = flatten(form_factors1, shape.count1);
  workspace.coordinates1.copy_from_host(coordinates1.data(), coordinates1.size());
  workspace.factors1.copy_from_host(flat_factors1.data(), flat_factors1.size());

  const double* device_coordinates2 = workspace.coordinates1.get();
  const double* device_factors2 = workspace.factors1.get();
  if (!same_particles) {
    const std::vector<double> flat_factors2 =
        flatten(form_factors2, shape.count2);
    workspace.coordinates2.copy_from_host(coordinates2.data(), coordinates2.size());
    workspace.factors2.copy_from_host(flat_factors2.data(), flat_factors2.size());
    device_coordinates2 = workspace.coordinates2.get();
    device_factors2 = workspace.factors2.get();
  }

  const std::size_t output_size = shape.output_channels * shape.bin_count;
  workspace.distributions.reserve(output_size);
  check_cuda(cudaMemset(workspace.distributions.get(), 0,
                        output_size * sizeof(double)),
             "cudaMemset distance distributions");
  if (shape.count1 == 0 || shape.count2 == 0) return shape;

  if (same_particles) {
    const unsigned long long tile_count =
        (shape.count1 + PAIR_TILE - 1) / PAIR_TILE;
    const unsigned long long block_count = tile_count * (tile_count + 1) / 2;
    make_same_distance_distributions_tiled<<<
        static_cast<unsigned int>(block_count), PAIR_THREADS>>>(
        workspace.coordinates1.get(), workspace.factors1.get(), shape.count1,
        shape.factor_channels, 1.0 / bin_size, shape.bin_count,
        workspace.distributions.get());
  } else {
    const dim3 blocks((shape.count1 + PAIR_TILE - 1) / PAIR_TILE,
                      (shape.count2 + PAIR_TILE - 1) / PAIR_TILE);
    make_cross_distance_distributions_tiled<<<blocks, PAIR_THREADS>>>(
        workspace.coordinates1.get(), workspace.factors1.get(), shape.count1,
        device_coordinates2, device_factors2, shape.count2,
        shape.factor_channels, 1.0 / bin_size, shape.bin_count,
        workspace.distributions.get());
  }
  check_cuda(cudaGetLastError(), "distance-distribution kernel launch");
  return shape;
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
  const DistributionShape shape = calculate_distributions_on_device(
      coordinates1, form_factors1, coordinates2, form_factors2,
      same_particles, bin_size);
  std::vector<double> flattened(shape.output_channels * shape.bin_count);
  workspace.distributions.copy_to_host(flattened.data(), flattened.size());
  distributions.assign(shape.output_channels,
                       std::vector<double>(shape.bin_count));
  for (std::size_t channel = 0; channel < shape.output_channels; ++channel) {
    std::copy(flattened.begin() + channel * shape.bin_count,
              flattened.begin() + (channel + 1) * shape.bin_count,
              distributions[channel].begin());
  }
}

void distance_distributions_to_profiles_cuda(
    const std::vector<double>& coordinates1,
    const std::vector<std::vector<double>>& form_factors1,
    const std::vector<double>& coordinates2,
    const std::vector<std::vector<double>>& form_factors2,
    bool same_particles, double bin_size, const std::vector<double>& q,
    double modulation_function_parameter,
    std::vector<std::vector<double>>& profiles) {
  const DistributionShape shape = calculate_distributions_on_device(
      coordinates1, form_factors1, coordinates2, form_factors2,
      same_particles, bin_size);
  workspace.q.copy_from_host(q.data(), q.size());
  const std::size_t profile_size = shape.output_channels * q.size();
  workspace.profiles.reserve(profile_size);
  if (!q.empty()) {
    transform_distance_distributions<<<q.size(), PROFILE_THREADS>>>(
        workspace.distributions.get(), shape.bin_count, shape.output_channels,
        bin_size, workspace.q.get(), q.size(), modulation_function_parameter,
        workspace.profiles.get());
    check_cuda(cudaGetLastError(), "distance-to-profile kernel launch");
  }
  std::vector<double> flattened(profile_size);
  workspace.profiles.copy_to_host(flattened.data(), flattened.size());
  profiles.assign(shape.output_channels, std::vector<double>(q.size()));
  for (std::size_t channel = 0; channel < shape.output_channels; ++channel) {
    std::copy(flattened.begin() + channel * q.size(),
              flattened.begin() + (channel + 1) * q.size(),
              profiles[channel].begin());
  }
}

void squared_distribution_2_profile_cuda(
    const double* radial_distribution, const double* q,
    const double* distances, double* intensity,
    double modulation_function_parameter, std::size_t radial_size,
    std::size_t q_size) {
  workspace.distributions.copy_from_host(radial_distribution, radial_size);
  workspace.distances.copy_from_host(distances, radial_size);
  workspace.q.copy_from_host(q, q_size);
  workspace.profiles.reserve(q_size);
  if (q_size > 0) {
    transform_precomputed_distances<<<q_size, PROFILE_THREADS>>>(
        workspace.distributions.get(), workspace.distances.get(), radial_size,
        workspace.q.get(), q_size, modulation_function_parameter,
        workspace.profiles.get());
    check_cuda(cudaGetLastError(), "profile kernel launch");
  }
  workspace.profiles.copy_to_host(intensity, q_size);
}

}  // namespace internal
}  // namespace saxs
}  // namespace foxs_cuda
