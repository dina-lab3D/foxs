/**
 * \file ChiScore.h \brief Basic chi score implementation
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef FOXS_CHI_SCORE_H
#define FOXS_CHI_SCORE_H

#include "foxs_config.h"
#include "Profile.h"

namespace foxs {

/**
   Basic implementation of Chi scoring
*/
class ChiScore {
 public:
  ChiScore() = default;

  // returns Chi_square score
  double compute_score(const Profile* exp_profile, const Profile* model_profile,
                       bool use_offset = false) const;

  double compute_scale_factor(const Profile* exp_profile,
                              const Profile* model_profile,
                              double offset = 0.0) const;

  double compute_offset(const Profile* exp_profile,
                        const Profile* model_profile) const;

};

}  // namespace foxs

#endif /* FOXS_CHI_SCORE_H */
