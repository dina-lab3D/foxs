#ifndef FOXS_CONFIG_H
#define FOXS_CONFIG_H

#include "Vector3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace foxs {

template <class T>
using Vector = std::vector<T>;

template <class T>
constexpr T square(const T& value) {
  return value * value;
}

constexpr double PI = 3.141592653589793238462643383279502884;

namespace algebra {

using Vector3D = ::Vector3;
using Vector3Ds = std::vector<Vector3D>;
class Vector2D {
 public:
  Vector2D(double x = 0.0, double y = 0.0) : values_{{x, y}} {}
  double operator[](std::size_t index) const { return values_[index]; }

 private:
  std::array<double, 2> values_;
};
using Vector2Ds = std::vector<Vector2D>;

inline double get_squared_distance(const Vector3D& first,
                                   const Vector3D& second) {
  return first.dist2(second);
}

inline double get_distance(const Vector3D& first, const Vector3D& second) {
  return std::sqrt(get_squared_distance(first, second));
}

inline unsigned int get_rounded(double value) {
  return static_cast<unsigned int>(std::lround(value));
}

class LinearFit2D {
 public:
  LinearFit2D(const Vector2Ds& points, const std::vector<double>& errors) {
    double sw = 0.0, sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
      const double error = i < errors.size() ? errors[i] : 1.0;
      const double weight = error != 0.0 ? 1.0 / square(error) : 1.0;
      const double x = points[i][0];
      const double y = points[i][1];
      sw += weight;
      sx += weight * x;
      sy += weight * y;
      sxx += weight * x * x;
      sxy += weight * x * y;
    }
    const double denominator = sw * sxx - sx * sx;
    if (points.empty() || std::abs(denominator) < 1e-30) {
      throw std::runtime_error("Cannot fit a line to the supplied points");
    }
    a_ = (sw * sxy - sx * sy) / denominator;
    b_ = (sy - a_ * sx) / sw;
  }

  double get_a() const { return a_; }
  double get_b() const { return b_; }

 private:
  double a_ = 0.0;
  double b_ = 0.0;
};

class ParabolicFit2D {
 public:
  explicit ParabolicFit2D(const Vector2Ds& points) {
    double matrix[3][4] = {};
    for (const Vector2D& point : points) {
      const double x = point[0];
      const double y = point[1];
      const double powers[5] = {1.0, x, x * x, x * x * x,
                                x * x * x * x};
      for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column)
          matrix[row][column] += powers[row + column];
        matrix[row][3] += y * powers[row];
      }
    }
    for (int pivot = 0; pivot < 3; ++pivot) {
      int best = pivot;
      for (int row = pivot + 1; row < 3; ++row)
        if (std::abs(matrix[row][pivot]) > std::abs(matrix[best][pivot]))
          best = row;
      if (std::abs(matrix[best][pivot]) < 1e-30)
        throw std::runtime_error("Cannot fit a parabola to the supplied points");
      for (int column = pivot; column < 4; ++column)
        std::swap(matrix[pivot][column], matrix[best][column]);
      const double divisor = matrix[pivot][pivot];
      for (int column = pivot; column < 4; ++column)
        matrix[pivot][column] /= divisor;
      for (int row = 0; row < 3; ++row) {
        if (row == pivot) continue;
        const double factor = matrix[row][pivot];
        for (int column = pivot; column < 4; ++column)
          matrix[row][column] -= factor * matrix[pivot][column];
      }
    }
    c_ = matrix[0][3];
    b_ = matrix[1][3];
    a_ = matrix[2][3];
  }

  double get_a() const { return a_; }
  double get_b() const { return b_; }
  double get_c() const { return c_; }

 private:
  double a_ = 0.0;
  double b_ = 0.0;
  double c_ = 0.0;
};

}  // namespace algebra

inline double get_squared_distance(const algebra::Vector3D& first,
                                   const algebra::Vector3D& second) {
  return algebra::get_squared_distance(first, second);
}

inline double get_distance(const algebra::Vector3D& first,
                           const algebra::Vector3D& second) {
  return algebra::get_distance(first, second);
}
}  // namespace foxs

#endif
