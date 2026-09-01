/**
 * \file cuda_helpers.h
 * \brief GPU implementations of some SAXS operations
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 *
 */

#ifndef FOXS_INTERNAL_CUDA_HELPERS_H
#define FOXS_INTERNAL_CUDA_HELPERS_H

namespace foxs_cuda {
namespace saxs {
namespace internal {

void squared_distribution_2_profile_cuda(
           const double *r_dist, const double *q,
           const double *distances, double *intensity,
           double modulation_function_parameter, size_t r_size, size_t q_size);

} } }

#endif /* FOXS_INTERNAL_CUDA_HELPERS_H */
