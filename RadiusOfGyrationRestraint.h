/**
 *  \file RadiusOfGyrationRestraint.h
 *  \brief Calculate score based on fit to SAXS profile.
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef FOXS_RADIUS_OF_GYRATION_RESTRAINT_H
#define FOXS_RADIUS_OF_GYRATION_RESTRAINT_H

#include "foxs_config.h"

#include "Profile.h"

namespace foxs {

//! Calculate score based on radius of gyration, taken from saxs profile
/** \ingroup exp_restraint

 */
class RadiusOfGyrationRestraint {
 public:
  //! Constructor
  /**
     \param[in] particles The particles participating in the fitting score
     \param[in] exp_profile  The experimental profile used in the fitting score
     \param[in] end_q_rg The range of profile used for approximation:
      i.e. q*rg < end_q_rg. Use 1.3 for globular proteins, 0.8 for elongated
  */
  RadiusOfGyrationRestraint(const Molecule<Atom>& particles,
                            const Profile* exp_profile,
                            const double end_q_rg = 1.3);

  double evaluate(Vector<algebra::Vector3D>* derivatives = nullptr) const;
  double unprotected_evaluate(
      Vector<algebra::Vector3D>* derivatives = nullptr) const {
    return evaluate(derivatives);
  }
  const Molecule<Atom>& get_inputs() const { return particles_; }

 private:
  Molecule<Atom> particles_;
  double exp_rg_;                // radius of gyration from experimental profile
};

}  // namespace foxs

#endif /* FOXS_RADIUS_OF_GYRATION_RESTRAINT_H */
