/**
 *  \file ChiScore.cpp   \brief Basic SAXS scoring
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#include "ChiScore.h"

namespace foxs {

double ChiScore::compute_score(const Profile* exp_profile,
                               const Profile* model_profile,
                               bool use_offset) const {
  double offset = 0.0;
  if (use_offset) offset = compute_offset(exp_profile, model_profile);
  double c = compute_scale_factor(exp_profile, model_profile, offset);

  double chi_square = 0.0;
  unsigned int profile_size =
    std::min(model_profile->size(), exp_profile->size());

  const std::vector<double>& errors = exp_profile->get_errors();
  const std::vector<double>& exp_intensities = exp_profile->get_intensities();
  const std::vector<double>& model_intensities = model_profile->get_intensities();

  for (unsigned int i = 0; i < profile_size; ++i) {
    const double delta =
        exp_intensities[i] - c * model_intensities[i] + offset;
    // Exclude the uncertainty originated from limitation of floating number
    if (fabs(delta / exp_intensities[i]) >= 1.0e-15)
      chi_square += square(delta) / square(errors[i]);
  }

  chi_square /= profile_size;
  return chi_square; //sqrt(chi_square);
}

double ChiScore::compute_scale_factor(const Profile* exp_profile,
                                      const Profile* model_profile,
                                      const double offset) const {

  const std::vector<double>& errors = exp_profile->get_errors();
  const std::vector<double>& exp_intensities = exp_profile->get_intensities();
  const std::vector<double>& model_intensities = model_profile->get_intensities();
  const unsigned int profile_size =
      std::min(model_profile->size(), exp_profile->size());
  double sum_imod2 = 0.0, sum_imod_iexp = 0.0, sum_imod = 0.0;
  for (unsigned int k = 0; k < profile_size; ++k) {
    const double weight = 1.0 / square(errors[k]);
    sum_imod2 += weight * square(model_intensities[k]);
    sum_imod_iexp += weight * model_intensities[k] * exp_intensities[k];
    sum_imod += weight * model_intensities[k];
  }
  double c =  sum_imod_iexp / sum_imod2;

  if(std::fabs(offset) > 0.0000000001) {
    double constant =  sum_imod /sum_imod2;
    c += offset*constant;
  }
  return c;
}

double ChiScore::compute_offset(const Profile* exp_profile,
                                const Profile* model_profile) const {

  const std::vector<double>& errors = exp_profile->get_errors();
  const std::vector<double>& exp_intensities = exp_profile->get_intensities();
  const std::vector<double>& model_intensities = model_profile->get_intensities();
  const unsigned int profile_size =
      std::min(model_profile->size(), exp_profile->size());
  double sum_imod = 0.0, sum_imod2 = 0.0, sum_imod_iexp = 0.0;
  for (unsigned int k = 0; k < profile_size; ++k) {
    const double weight = 1.0 / square(errors[k]);
    sum_imod += weight * model_intensities[k];
    sum_imod2 += weight * square(model_intensities[k]);
    sum_imod_iexp += weight * model_intensities[k] * exp_intensities[k];
  }
  double constant = sum_imod /sum_imod2;

  // compute scaling
  double c =  sum_imod_iexp / sum_imod2;

  // compute offset
  double sum1 = 0.0, sum2 = 0.0;
  for (unsigned int k = 0; k < profile_size; ++k) {
    const double weight = 1.0 / square(errors[k]);
    const double delta = exp_intensities[k] - c * model_intensities[k];
    const double delta2 = 1.0 - constant * model_intensities[k];
    sum1 += weight * delta * delta2;
    sum2 += weight * square(delta2);
  }
  double offset = -sum1/sum2;

  return offset;
}

}  // namespace foxs
