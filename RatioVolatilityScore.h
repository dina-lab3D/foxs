/**
 * \file RatioVolatilityScore.h \brief Vr score implementation
 *
 * Hura et al. Nature Methods 2013
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef FOXS_RATIO_VOLATILITY_SCORE_H
#define FOXS_RATIO_VOLATILITY_SCORE_H

#include "foxs_config.h"
#include "Profile.h"

namespace foxs {

/**
   Basic implementation of RatioVolatility scoring
*/
class RatioVolatilityScore {
 public:
  explicit RatioVolatilityScore(double dmax = 400) : dmax_(dmax) {}

  double compute_score(const Profile* exp_profile, const Profile* model_profile,
                       bool use_offset = false) const;

  double compute_scale_factor(const Profile* exp_profile,
                              const Profile* model_profile,
                              double offset = 0.0) const;

  double compute_offset(const Profile* exp_profile,
                        const Profile* model_profile) const {
    // not implemented for now
    (void)exp_profile;
    (void)model_profile;
    return 0.0;
  }
 private:
  double dmax_;
};

}  // namespace foxs

#endif /* FOXS_RATIO_VOLATILITY_SCORE_H */
