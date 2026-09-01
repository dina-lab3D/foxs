/**
 * \file ChiScoreLog.h \brief scoring with log intensity
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef FOXS_CHI_SCORE_LOG_H
#define FOXS_CHI_SCORE_LOG_H

#include "foxs_config.h"
#include "Profile.h"

namespace foxs {

/**
   chi scoring on log intensities
*/
class ChiScoreLog {
 public:
  ChiScoreLog() = default;

  double compute_score(const Profile* exp_profile, const Profile* model_profile,
                       bool use_offset = false) const;

  double compute_score(const Profile* exp_profile, const Profile* model_profile,
                       double min_q, double max_q) const;

  double compute_scale_factor(const Profile* exp_profile,
                              const Profile* model_profile,
                              double offset = 0.0) const;

  double compute_offset(const Profile* exp_profile,
                        const Profile* model_profile) const {
    // not implemented as no straightforward solution to the equations
    (void)exp_profile;
    (void)model_profile;
    return 0.0;
  }

};

}  // namespace foxs

#endif /* FOXS_CHI_SCORE_LOG_H */
