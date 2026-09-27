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

//! Compute per-atom solvent-accessible surface areas on the GPU.
/**
 * Coordinates are packed xyz triples. A sparse uniform grid restricts each
 * surface-point intersection test to nearby atoms.
 */
void solvent_accessible_surface_areas_cuda(
    const std::vector<float>& coordinates, const std::vector<float>& radii,
    float probe_radius, float density, std::vector<float>& areas);

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

//! Calculate distance distributions and transform them to SAXS profiles
//! without copying the intermediate histograms back to the host.
void distance_distributions_to_profiles_cuda(
    const std::vector<double>& coordinates1,
    const std::vector<std::vector<double>>& form_factors1,
    const std::vector<double>& coordinates2,
    const std::vector<std::vector<double>>& form_factors2,
    bool same_particles, double bin_size, const std::vector<double>& q,
    double modulation_function_parameter,
    std::vector<std::vector<double>>& profiles);

void squared_distribution_2_profile_cuda(
           const double *r_dist, const double *q,
           const double *distances, double *intensity,
           double modulation_function_parameter, size_t r_size, size_t q_size);

} } }

#endif /* FOXS_INTERNAL_CUDA_HELPERS_H */
