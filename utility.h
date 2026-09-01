/**
 *  \file utility.h
 *  \brief Functions to deal with very common saxs operations
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
*/

#ifndef FOXS_UTILITY_H
#define FOXS_UTILITY_H

#include "foxs_config.h"
#include "FormFactorTable.h"
#include "Profile.h"

#include "Atom.h"
#include "Molecule.h"

namespace foxs {

inline void get_coordinates(const Molecule<Atom>& particles,
                            Vector<algebra::Vector3D>& coordinates) {
  // copy everything in advance for fast access
  coordinates.resize(particles.size());
  for (unsigned int i = 0; i < particles.size(); i++) {
    coordinates[i] = algebra::Vector3D(particles[i][0], particles[i][1],
                                       particles[i][2]);
  }
}

inline void get_form_factors(const Molecule<Atom>& particles,
                             FormFactorTable* ff_table,
                             Vector<double>& form_factors,
                             FormFactorType ff_type) {
  form_factors.resize(particles.size());
  for (unsigned int i = 0; i < particles.size(); i++) {
    form_factors[i] = ff_table->get_form_factor(particles[i], ff_type);
  }
}

//! compute max distance
inline double compute_max_distance(const Molecule<Atom>& particles) {
  double max_dist2 = 0;
  Vector<algebra::Vector3D> coordinates(particles.size());
  get_coordinates(particles, coordinates);
  for (unsigned int i = 0; i < coordinates.size(); i++) {
    for (unsigned int j = i + 1; j < coordinates.size(); j++) {
      double dist2 = algebra::get_squared_distance(coordinates[i], coordinates[j]);
      if (dist2 > max_dist2) max_dist2 = dist2;
    }
  }
  return std::sqrt(max_dist2);
}

//! compute max distance between pairs of particles one from particles1
//! and the other from particles2
inline double compute_max_distance(const Molecule<Atom>& particles1,
                                  const Molecule<Atom>& particles2) {
  double max_dist2 = 0;
  Vector<algebra::Vector3D> coordinates1, coordinates2;
  get_coordinates(particles1, coordinates1);
  get_coordinates(particles2, coordinates2);

  for (unsigned int i = 0; i < coordinates1.size(); i++) {
    for (unsigned int j = i + 1; j < coordinates2.size(); j++) {
      double dist2 = algebra::get_squared_distance(coordinates1[i], coordinates2[j]);
      if (dist2 > max_dist2) max_dist2 = dist2;
    }
  }
  return std::sqrt(max_dist2);
}

//! compute radius_of_gyration
inline double radius_of_gyration(const Molecule<Atom>& particles) {
  algebra::Vector3D centroid(0.0, 0.0, 0.0);
  Vector<algebra::Vector3D> coordinates(particles.size());
  get_coordinates(particles, coordinates);
  for (unsigned int i = 0; i < particles.size(); i++) {
    centroid += coordinates[i];
  }
  centroid /= particles.size();
  double rg = 0;
  for (unsigned int i = 0; i < particles.size(); i++) {
    rg += algebra::get_squared_distance(coordinates[i], centroid);
  }
  rg /= particles.size();
  return std::sqrt(rg);
}

//! profile calculation for particles and a given set of options
Profile compute_profile(Molecule<Atom> particles,
                         double min_q = 0.0, double max_q = 0.5,
                         double delta_q = 0.001,
                         FormFactorTable* ft = get_default_form_factor_table(),
                         FormFactorType ff_type = HEAVY_ATOMS,
                         bool hydration_layer = true,
                         bool fit = true,
                         bool reciprocal = false,
                         bool ab_initio = false,
                         bool vacuum = false,
                         std::string beam_profile_file = "",
                         bool use_gpu = false);

//! Read PDB (or mmCIF) files
void read_pdb(const std::string& file,
              std::vector<std::string>& pdb_file_names,
              std::vector<Molecule<Atom>>& particles_vec,
              bool residue_level = false,
              bool heavy_atoms_only = true,
              int multi_model_pdb = 2,
              bool explicit_water = false);

//! Parse PDB and profile files
void read_files(const std::vector<std::string>& files,
                std::vector<std::string>& pdb_file_names,
                std::vector<std::string>& dat_files,
                std::vector<Molecule<Atom>>& particles_vec,
                Profiles& exp_profiles,
                bool residue_level = false,
                bool heavy_atoms_only = true,
                int multi_model_pdb = 2,
                bool explicit_water = false,
                float max_q = 0.0,
                int units = 1);

std::string trim_extension(const std::string file_name);


}  // namespace foxs

#endif /* FOXS_UTILITY_H */
