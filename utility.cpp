/**
 *  \file utility.cpp
 *  \brief Functions to deal with very common SAXS operations using GAMB.
 */

#include "utility.h"
#include "SolventAccessibleSurface.h"

#include <fstream>
#include <memory>
#include <sstream>

namespace foxs {

namespace {

class FoxsPDBSelector : public PDB::Selector {
 public:
  FoxsPDBSelector(bool residue_level, bool heavy_atoms_only,
                  bool explicit_water)
      : residue_level_(residue_level),
        heavy_atoms_only_(heavy_atoms_only),
        explicit_water_(explicit_water) {}

  bool operator()(const char* record) const override {
    const std::string line(record);
    const char alternate = PDB::atomAltLocIndicator(line);
    if (alternate != ' ' && alternate != 'A') return false;
    if (residue_level_) return PDB::CAlphaSelector()(record);

    const Atom atom(line);
    if (heavy_atoms_only_ && atom.isH()) return false;
    const std::string residue(atom.residueName());
    const bool water = residue == "HOH" || residue == "DOD" || residue == "WAT";
    return explicit_water_ || !water;
  }

 private:
  bool residue_level_;
  bool heavy_atoms_only_;
  bool explicit_water_;
};

}  // namespace

Profile compute_profile(Molecule<Atom> particles, double min_q,
                        double max_q, double delta_q, FormFactorTable* ft,
                        FormFactorType ff_type, bool hydration_layer, bool fit,
                        bool reciprocal, bool ab_initio, bool vacuum,
                        std::string beam_profile_file) {
  Profile profile(min_q, max_q, delta_q);
  if (reciprocal) profile.set_ff_table(ft);
  if (!beam_profile_file.empty()) profile.set_beam_profile(beam_profile_file);

  Vector<double> surface_area;
  SolventAccessibleSurface surface_calculator;
  double average_radius = 0.0;
  if (hydration_layer && !particles.empty()) {
    std::vector<double> radii;
    radii.reserve(particles.size());
    for (const Atom& atom : particles) {
      const double radius = ft->get_radius(atom, ff_type);
      radii.push_back(radius);
      average_radius += radius;
    }
    surface_area =
        surface_calculator.get_solvent_accessibility(particles, radii);
    profile.set_average_radius(average_radius / particles.size());
  }

  if (!fit) {
    if (ab_initio) {
      profile.calculate_profile_constant_form_factor(particles);
    } else if (vacuum) {
      profile.calculate_profile_partial(particles, surface_area, ff_type);
      profile.sum_partial_profiles(0.0, 0.0);
    } else {
      profile.calculate_profile(particles, ff_type, reciprocal);
    }
  } else if (reciprocal) {
    profile.calculate_profile_reciprocal_partial(particles, surface_area,
                                                 ff_type);
  } else {
    profile.calculate_profile_partial(particles, surface_area, ff_type);
  }
  return profile;
}

void read_pdb(const std::string& file,
              std::vector<std::string>& pdb_file_names,
              std::vector<Molecule<Atom>>& particles_vec,
              bool residue_level, bool heavy_atoms_only, int multi_model_pdb,
              bool explicit_water) {
  std::ifstream input(file.c_str());
  if (!input) return;

  FoxsPDBSelector selector(residue_level, heavy_atoms_only, explicit_water);
  std::vector<Molecule<Atom>> models;

  if (multi_model_pdb == 2) {
    while (input) {
      Molecule<Atom> model;
      model.readModelFromPDBfile(input, selector);
      if (!model.empty()) models.push_back(std::move(model));
    }
  } else {
    Molecule<Atom> model;
    if (multi_model_pdb == 3) {
      if (explicit_water)
        model.readAllPDBfile(input, selector);
      else
        model.readPDBfile(input, selector);
    } else {
      model.readModelFromPDBfile(input, selector);
    }
    if (!model.empty()) models.push_back(std::move(model));
  }

  const bool mmcif = file.size() >= 4 && file.substr(file.size() - 4) == ".cif";
  for (std::size_t index = 0; index < models.size(); ++index) {
    std::string name = file;
    if (models.size() > 1) {
      std::ostringstream numbered;
      numbered << trim_extension(file) << "_m" << index + 1
               << (mmcif ? ".cif" : ".pdb");
      name = numbered.str();
    }
    pdb_file_names.push_back(name);
    particles_vec.push_back(std::move(models[index]));
  }
}

void read_files(const std::vector<std::string>& files,
                std::vector<std::string>& pdb_file_names,
                std::vector<std::string>& dat_files,
                std::vector<Molecule<Atom>>& particles_vec,
                Profiles& exp_profiles, bool residue_level,
                bool heavy_atoms_only, int multi_model_pdb,
                bool explicit_water, float max_q, int units) {
  for (const std::string& file : files) {
    std::ifstream input(file.c_str());
    if (!input) {
      std::cerr << "Can't open file " << file << std::endl;
      continue;
    }

    const std::size_t old_structure_count = particles_vec.size();
    read_pdb(file, pdb_file_names, particles_vec, residue_level,
             heavy_atoms_only, multi_model_pdb, explicit_water);
    if (particles_vec.size() != old_structure_count) continue;

    std::shared_ptr<Profile> profile =
        std::make_shared<Profile>(file, false, max_q, units);
    if (profile->empty()) {
      std::cerr << "can't parse input file " << file << std::endl;
    } else {
      dat_files.push_back(file);
      exp_profiles.push_back(std::move(profile));
    }
  }
}

std::string trim_extension(const std::string file_name) {
  if (file_name.size() >= 4 && file_name[file_name.size() - 4] == '.') {
    return file_name.substr(0, file_name.size() - 4);
  } else if (file_name.size() >= 5 && file_name[file_name.size() - 5] == '.') {
    return file_name.substr(0, file_name.size() - 5);
  }
  return file_name;
}

}  // namespace foxs
