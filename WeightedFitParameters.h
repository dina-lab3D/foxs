/**
 * \file WeightedFitParameters.h
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef FOXS_WEIGHTED_FIT_PARAMETERS_H
#define FOXS_WEIGHTED_FIT_PARAMETERS_H

#include "foxs_config.h"

#include "FitParameters.h"
#include <iostream>

namespace foxs {

//! Parameters of a weighted fit, from WeightedProfileFitter.
class WeightedFitParameters : public FitParameters {
 public:
  WeightedFitParameters() : FitParameters() {}

  WeightedFitParameters(double chi_square, double c1, double c2,
                        const Vector<double>& weights = Vector<double>())
      : FitParameters(chi_square, c1, c2), weights_(weights) {}

  WeightedFitParameters(const FitParameters& fp) :
    FitParameters(fp) {}

  const Vector<double>& get_weights() const { return weights_; }

  void set_weights(const Vector<double>& weights) { weights_ = weights; }

  void show(std::ostream& s) const {
    s << "Chi^2 = " << chi_square_ << " c1 = " << c1_ << " c2 = " << c2_
      << " default chi^2 = " << default_chi_square_ << std::endl;
  }

 private:
  Vector<double> weights_;
};

}  // namespace foxs

#endif /* FOXS_WEIGHTED_FIT_PARAMETERS_H */
