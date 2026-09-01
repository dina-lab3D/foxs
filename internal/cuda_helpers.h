/**
 * \file cuda_helpers.h
 * \brief GPU implementations of some SAXS operations
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 *
 */

#ifndef FOXS_INTERNAL_CUDA_HELPERS_H
#define FOXS_INTERNAL_CUDA_HELPERS_H

#include <cstddef>
#include <string>
#include <vector>

namespace foxs_cuda {
namespace saxs {
namespace internal {

//! Return true when at least one CUDA device is usable.
bool cuda_device_available(std::string* reason = nullptr);

//! Calculate one, three, or six weighted squared-distance histograms.
/**
 * Coordinates are packed xyz triples. Form factors are channel-major, with
 * one channel for a normal profile and two or three channels for a partial
 * profile. When same_particles is true, self terms are included and each
 * distinct pair is counted twice, matching the CPU Debye calculation.
 */
void distance_distributions_cuda(
    const std::vector<double>& coordinates1,
    const std::vector<std::vector<double>>& form_factors1,
    const std::vector<double>& coordinates2,
    const std::vector<std::vector<double>>& form_factors2,
    bool same_particles, double bin_size,
    std::vector<std::vector<double>>& distributions);

void squared_distribution_2_profile_cuda(
           const double *r_dist, const double *q,
           const double *distances, double *intensity,
           double modulation_function_parameter, size_t r_size, size_t q_size);

} } }

#endif /* FOXS_INTERNAL_CUDA_HELPERS_H */
