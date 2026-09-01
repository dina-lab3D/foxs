/**
 *  \file RadiusOfGyrationRestraint.h
 *  \brief Calculate score based on fit to SAXS profile.
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#include "RadiusOfGyrationRestraint.h"
#include "utility.h"

namespace foxs {

RadiusOfGyrationRestraint::RadiusOfGyrationRestraint(
                                                     const Molecule<Atom>& particles,
                                                     const Profile* exp_profile,
                                                     const double end_q_rg)
  : particles_(particles) {
  exp_rg_ = exp_profile->radius_of_gyration(end_q_rg);
}

//! Calculate the score and the derivatives for particles of the restraint.
/** \param[in] acc If true (not nullptr), partial first derivatives should be
    calculated.
    \return score associated with this restraint for the given state of
            the model.
*/
double RadiusOfGyrationRestraint::evaluate(
    Vector<algebra::Vector3D>* derivatives) const {

  // get centroid
  algebra::Vector3D centroid(0.0, 0.0, 0.0);
  Vector<algebra::Vector3D> coordinates(particles_.size());
  get_coordinates(particles_, coordinates);
  for (unsigned int i = 0; i < particles_.size(); i++) {
    centroid += coordinates[i];
  }
  centroid /= particles_.size();
  double radg = 0;
  for (unsigned int i = 0; i < particles_.size(); i++) {
    radg += get_squared_distance(coordinates[i], centroid);
  }
  radg /= particles_.size();
  radg = sqrt(radg);

  double score = (radg - exp_rg_) / exp_rg_;  // TODO: improve
  if (!derivatives) return score;

  double factor = 1.0 / (particles_.size() * radg);
  derivatives->assign(particles_.size(), algebra::Vector3D());
  for (unsigned int i = 0; i < particles_.size(); i++) {
    (*derivatives)[i] = (coordinates[i] - centroid) * factor;
  }
  return score;
}

}  // namespace foxs
