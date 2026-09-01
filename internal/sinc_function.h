/**
 * \file internal/sinc_function.h
 * \brief caching of sinc values
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef FOXS_INTERNAL_SINC_FUNCTION_H
#define FOXS_INTERNAL_SINC_FUNCTION_H

#include "foxs_config.h"
#include <boost/math/special_functions/sinc.hpp>
#include <cmath>

namespace foxs { namespace internal {

class SincFunction : public std::vector<float> {
 public:
  // Constructor
  SincFunction(float max_value, float bin_size) {
    bin_size_ = bin_size;
    one_over_bin_size_ = 1.0 / bin_size_;
    max_value_ = max_value;
    unsigned int size = value2index(max_value_) + 1;
    reserve(size);
    for (unsigned int i = 0; i <= size; i++) {
      float x = i * bin_size_;
      push_back(boost::math::sinc_pi(x));
    }
  }

  unsigned int value2index(float value) const {
    return static_cast<unsigned int>(std::lround(value * one_over_bin_size_));
  }

  // get sinc value for x, compute values if they weren't computed yet
  float sinc(float x) {
    unsigned int index = value2index(x);
    if (index >= size()) {
      reserve(index);
      for (unsigned int i = size(); i <= index; i++) {
        float x = i * bin_size_;
        push_back(boost::math::sinc_pi(x));
      }
    }
    return (*this)[index];
  }

 private:
  float bin_size_, one_over_bin_size_;  // resolution of discretization
  float max_value_;
};

class SincCosFunction : public std::vector<float> {
 public:
  // Constructor
  SincCosFunction(float max_value, float bin_size) {
    bin_size_ = bin_size;
    one_over_bin_size_ = 1.0 / bin_size_;
    max_value_ = max_value;
    unsigned int size = value2index(max_value_) + 1;
    reserve(size);
    for (unsigned int i = 0; i <= size; i++) {
      float x = i * bin_size_;
      push_back((boost::math::sinc_pi(x) - std::cos(x)) / x);
    }
  }

  unsigned int value2index(float value) const {
    return static_cast<unsigned int>(std::lround(value * one_over_bin_size_));
  }

  // get value for x, compute values if they weren't computed yet
  float sico(float x) {
    unsigned int index = value2index(x);
    if (index >= size()) {
      reserve(index);
      for (unsigned int i = size(); i <= index; i++) {
        float x = i * bin_size_;
        push_back((boost::math::sinc_pi(x) - std::cos(x)) / x);
      }
    }
    return (*this)[index];
  }

 private:
  float bin_size_, one_over_bin_size_;  // resolution of discretization
  float max_value_;
};

} }  // namespace foxs::internal

#endif /* FOXS_INTERNAL_SINC_FUNCTION_H */
