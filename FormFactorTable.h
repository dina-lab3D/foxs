/**
 *  \file FormFactorTable.h   \brief A class for computation of
 * atomic and residue level form factors for SAXS calculations
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef FOXS_FORM_FACTOR_TABLE_H
#define FOXS_FORM_FACTOR_TABLE_H

#include "foxs_config.h"

#include "Atom.h"

#include <iostream>
#include <vector>


namespace foxs {

//! type of the form factors for profile calculations
/*
 ALL_ATOMS - all atoms including hydrogens
 HEAVY_ATOMS - no hydrogens, all other atoms included
 CA_ATOMS - residue level, residue represented by CA
 RESIDUES - residue level, represented by residue bead
*/
enum FormFactorType {
  ALL_ATOMS,
  HEAVY_ATOMS,
  CA_ATOMS,
  RESIDUES
};

/**
   class that deals with form factor computation
   two form factors are supported:
   (i) zero form factors for faster approximated calculations
   (ii) full form factors for slower accurate calculations

   Each form factor can be divided into two parts: vacuum and dummy.
   dummy is an approximated excluded volume (solvent) form factor.
   The approximation is done using Fraser, MacRae and Suzuki (1978) model.
*/
class FormFactorTable {
 public:
  //! default constructor
  FormFactorTable();

  //! constructor with form factor table file (required for full form factors)
  FormFactorTable(const std::string& table_name, double min_q, double max_q,
                  double delta_q);

  // 1. Zero form factors

  //! get f(0), ie q=0 for real space profile calculation
  double get_form_factor(const Atom& atom,
                         FormFactorType ff_type = HEAVY_ATOMS) const;

  //! f(0) in vacuum
  double get_vacuum_form_factor(const Atom& atom,
                                FormFactorType ff_type = HEAVY_ATOMS) const;

  //! f(0) for solvent
  double get_dummy_form_factor(const Atom& atom,
                               FormFactorType ff_type = HEAVY_ATOMS) const;

  //! f(0) for water
  double get_water_form_factor() const { return zero_form_factors_[OH2]; }

  //! f(0) for water in vacuum
  double get_vacuum_water_form_factor() const {
    return vacuum_zero_form_factors_[OH2];
  }

  //! f(0) for water (solvent)
  double get_dummy_water_form_factor() const {
    return dummy_zero_form_factors_[OH2];
  }

  // 2. Full form factors

  //! full form factor for reciprocal space profile calculation
  const std::vector<double>& get_form_factors(const Atom& atom,
                                 FormFactorType ff_type = HEAVY_ATOMS) const;

  //! for reciprocal space profile calculation
  const std::vector<double>& get_vacuum_form_factors(const Atom& atom,
                                 FormFactorType ff_type = HEAVY_ATOMS) const;

  //! for reciprocal space profile calculation
  const std::vector<double>& get_dummy_form_factors(const Atom& atom,
                                 FormFactorType ff_type = HEAVY_ATOMS) const;

  //! full water form factor
  const std::vector<double>& get_water_form_factors() const {
    return form_factors_[OH2];
  }

  //! full water vacuum form factor
  const std::vector<double>& get_water_vacuum_form_factors() const {
    return vacuum_form_factors_[OH2];
  }

  //! full water dummy form factor
  const std::vector<double>& get_water_dummy_form_factors() const {
    return dummy_form_factors_[OH2];
  }

  //! radius
  double get_radius(const Atom& atom,
                    FormFactorType ff_type = HEAVY_ATOMS) const;

  //! volume
  double get_volume(const Atom& atom,
                    FormFactorType ff_type = HEAVY_ATOMS) const;

  //! print tables
  void show(std::ostream& out = std::cout, std::string prefix = "") const;

  // electron density of solvent - default=0.334 e/A^3 (H2O)
  static double rho_;

 private:
  // atom types for heavy atoms according to the number of hydrogens
  // connected to them
  // ALL_ATOM_SIZE is number of types needed for all atom representation
  // this indexing is used in form_factors arrays
  enum FormFactorAtomType {
    H,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne,  // periodic table, lines 1-2 (10)
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Cl,
    Ar,  // line 3 (8)
    K,
    Ca,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Se,
    Br,  // line 4 (11)
    Ag,
    I,
    Ir,
    Pt,
    Au,
    Hg,
    ALL_ATOM_SIZE = 35,
    CH = 35,
    CH2 = 36,
    CH3 = 37,
    NH = 38,
    NH2 = 39,
    NH3 = 40,
    OH = 41,
    OH2 = 42,
    SH = 43,
    HEAVY_ATOM_SIZE = 44,
    UNK = 45
  };

  // map between atom element and FormFactorAtomType
  static std::map<std::string, FormFactorAtomType> element_ff_type_map_;

  struct FormFactor {
    FormFactor() {}
    FormFactor(double ff, double vacuum_ff, double dummy_ff)
        : ff_(ff), vacuum_ff_(vacuum_ff), dummy_ff_(dummy_ff) {}
    double ff_, vacuum_ff_, dummy_ff_;
  };

  // map between residue type and residue level form factors
  static std::map<std::string, FormFactor> residue_type_form_factor_map_;

  // form factors for q=0, the order as in the FormFactorAtomType enum
  static double zero_form_factors_[];

  static double vacuum_zero_form_factors_[];
  // those represent excluded volume
  static double dummy_zero_form_factors_[];

  // class for storing form factors solvation table
  class AtomFactorCoefficients {
   public:
    std::string atom_type_;
    double a_[5];
    double b_[5];
    double c_;
    double excl_vol_;
  };

#ifndef SWIG
  // read entry
  friend std::istream& operator>>(
      std::istream& s, AtomFactorCoefficients& atom_factor_coefficients);

  // write entry
  friend std::ostream& operator<<(
      std::ostream& s, const AtomFactorCoefficients& atom_factor_coefficients);
#endif

 private:
  int read_form_factor_table(const std::string& table_name);

  void init_element_form_factor_map();

  void init_residue_type_form_factor_map();

  void compute_form_factors_all_atoms();

  void compute_form_factors_heavy_atoms();

  double get_form_factor(const std::string& residue_type) const;

  double get_vacuum_form_factor(const std::string& residue_type) const;

  double get_dummy_form_factor(const std::string& residue_type) const;

  FormFactorAtomType get_form_factor_atom_type(
      const std::string& element) const;

  FormFactorAtomType get_form_factor_atom_type(const Atom& atom,
                                               FormFactorType ff_type) const;

  FormFactorAtomType get_carbon_atom_type(const std::string& atom_type,
                                          const std::string& residue_type)
      const;

  FormFactorAtomType get_nitrogen_atom_type(
      const std::string& atom_type,
      const std::string& residue_type) const;

  FormFactorAtomType get_oxygen_atom_type(const std::string& atom_type,
                                          const std::string& residue_type)
      const;

  FormFactorAtomType get_sulfur_atom_type(const std::string& atom_type,
                                          const std::string& residue_type)
      const;

 private:
  // read from lib file
  std::vector<AtomFactorCoefficients> form_factors_coefficients_;

  // table of full form factors for 14 atom types
  std::vector<std::vector<double>> form_factors_;

  // vacuum full form factors for 14 atom types
  std::vector<std::vector<double>> vacuum_form_factors_;

  // dummy full form factors for 14 atom types
  std::vector<std::vector<double>> dummy_form_factors_;

  // min/max q and sampling resolution for form factor computation
  double min_q_, max_q_, delta_q_;

};

/** Get the default table.*/
FormFactorTable* get_default_form_factor_table();

}  // namespace foxs

#endif /* FOXS_FORM_FACTOR_TABLE_H */
