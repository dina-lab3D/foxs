/**
 * \file ChiFreeScore.h \brief Chi free score implementation
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 * TODO: add reference
 */

#ifndef FOXS_CHI_FREE_SCORE_H
#define FOXS_CHI_FREE_SCORE_H

#include "foxs_config.h"
#include "Profile.h"

namespace foxs {
/**
   Implementation of chi free score
   Accurate assessment of mass, models and resolution by small-angle scattering.
   Rambo RP, Tainer JA. Nature. 2013
*/
class ChiFreeScore {
 public:
  ChiFreeScore(unsigned int ns, unsigned int k) :
            ns_(ns), K_(k) {
    if (K_ % 2 == 0) K_++;  // make sure it is odd for median
    last_scale_updated_ = false;
  }

  double compute_score(const Profile* exp_profile, const Profile* model_profile,
                       bool use_offset = false) const;

  double compute_scale_factor(const Profile* exp_profile,
                              const Profile* model_profile,
                              double offset = 0.0) const;

  double compute_offset(const Profile* exp_profile,
                        const Profile* model_profile) const;

 private:
  unsigned int ns_;  // number of Shannon channels
  unsigned int K_;
  double last_scale_;
  bool last_scale_updated_;
};

}  // namespace foxs

#endif /* FOXS_CHI_FREE_SCORE_H */
