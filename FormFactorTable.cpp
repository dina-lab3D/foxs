/**
 *  \file FormFactorTable.h   \brief A class for computation of
 * atomic and residue level form factors for SAXS calculations
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#include "FormFactorTable.h"
#include <fstream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace foxs {

namespace {

std::string normalize_name(std::string value) {
  value.erase(std::remove_if(value.begin(), value.end(),
                             [](unsigned char c) { return std::isspace(c); }),
              value.end());
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::toupper(c));
                 });
  std::replace(value.begin(), value.end(), '*', '\'');
  return value;
}

bool is_one_of(const std::string& value,
               std::initializer_list<const char*> choices) {
  for (const char* choice : choices) {
    if (value == choice) return true;
  }
  return false;
}

bool is_adenine(const std::string& residue) {
  return is_one_of(residue, {"A", "ADE", "DA", "DADE"});
}
bool is_cytosine(const std::string& residue) {
  return is_one_of(residue, {"C", "CYT", "DC", "DCYT"});
}
bool is_uracil(const std::string& residue) {
  return is_one_of(residue, {"U", "URA", "DU", "DURA"});
}
bool is_thymine(const std::string& residue) {
  return is_one_of(residue, {"T", "THY", "DT", "DTHY"});
}
bool is_guanine(const std::string& residue) {
  return is_one_of(residue, {"G", "GUA", "DG", "DGUA"});
}

}  // namespace

std::map<std::string, FormFactorTable::FormFactorAtomType>
    FormFactorTable::element_ff_type_map_;

std::map<std::string, FormFactorTable::FormFactor>
    FormFactorTable::residue_type_form_factor_map_;

double FormFactorTable::zero_form_factors_[] = {
    -0.720147, -0.720228,
    //   H       He - periodic table line 1
    1.591,     2.591,     3.591,   0.50824,  6.16294, 4.94998, 7.591,   6.993,
    // Li     Be      B     C       N        O       F      Ne - line 2
    7.9864,    8.9805,    9.984,   10.984,   13.0855, 9.36656, 13.984,  16.591,
    //  Na      Mg        Al       Si        P        S       Cl    Ar - line 3
    15.984,    14.9965,   20.984,  21.984,   20.9946, 23.984,
    // K       Ca2+       Cr      Mn      Fe2+      Co - line 4
    24.984,   25.984,     24.9936, 30.9825,  31.984,  43.984, 49.16,
    // Ni     Cu          Zn2+      Se       Br       Ag      I
    70.35676,  71.35676,  72.324,  73.35676,
    // Ir         Pt      Au      Hg
    -0.211907, -0.932054, -1.6522, 5.44279,  4.72265, 4.0025,  4.22983, 3.50968,
    8.64641
    //  CH        CH2        CH3     NH       NH2       NH3     OH       OH2
    // SH
};

double FormFactorTable::vacuum_zero_form_factors_[] = {
    //   H       He - periodic table line 1
    0.999953, 0.999872,
    // Li  Be    B     C       N       O       F     Ne - line 2
    2.99,  3.99, 4.99, 5.9992, 6.9946, 7.9994, 8.99, 9.999,
    //  Na     Mg     Al     Si      P        S        Cl     Ar - line 3
    10.9924, 11.9865, 12.99, 13.99,  14.9993, 15.9998, 16.99, 17.99,
    // K    Ca2+     Cr     Mn     Fe2+     Co - line 4
    18.99,  18.0025, 23.99, 24.99, 24.0006, 26.99,
    // Ni   Cu      Zn2+     Se     Br - line 4 cont.
    27.99,  28.99,  27.9996, 33.99, 34.99,
    // Ag    I       Ir     Pt      Au     Hg - some elements from lines 5, 6
    46.99, 52.99,   76.99,   77.99, 78.9572,  79.99,
    // CH      CH2     CH3     NH       NH2       NH3     OH      OH2      SH
    6.99915,  7.99911,  8.99906, 7.99455, 8.99451, 9.99446, 8.99935, 9.9993,
    16.9998
};

double FormFactorTable::dummy_zero_form_factors_[] = {
    1.7201,  1.7201,  1.399,   1.399,   1.399,   5.49096, 0.83166, 3.04942,
    1.399,   3.006,
    //  H     He     Li?    Be?    B?       C        N        O      F?     Ne
    3.006,   3.006,   3.006,   3.006,   1.91382, 6.63324, 3.006,   1.399,
    // Na     Mg    Al?    Si?      P        S      Cl?    Ar?
    3.006,   3.006,   3.006,   3.006,   3.006,   3.006,
    // K?   Ca2+    Cr?    Mn?   Fe2+   Co?
    3.006,   3.006,   3.006,   3.006,   3.006,
    // Ni?   Cu?   Zn2+    Se     Br?
    3.006,   3.83,    6.63324, 6.63324, 6.63324, 6.63324,
    // Ag?   I?       Ir?      Pt?       Au      Hg
    7.21106, 8.93116, 10.6513, 2.55176, 4.27186, 5.99196, 4.76952, 6.48962,
    8.35334
    //  CH       CH2      CH3     NH       NH2       NH3     OH       OH2   SH
};

// electron density of solvent - default=0.334 e/A^3 (H2O)
double FormFactorTable::rho_ = 0.334;

std::istream& operator>>(
    std::istream& s,
    FormFactorTable::AtomFactorCoefficients& atom_factor_coefficients) {
  // read atom type
  s >> atom_factor_coefficients.atom_type_;
  // read ai values
  for (unsigned int i = 0; i < 5; i++) {
    s >> atom_factor_coefficients.a_[i];
  }
  s >> atom_factor_coefficients.c_;  // c value
  // read bi values
  for (unsigned int i = 0; i < 5; i++) {
    s >> atom_factor_coefficients.b_[i];
  }
  return s >> atom_factor_coefficients.excl_vol_;  // excluded volume
}

std::ostream& operator<<(
    std::ostream& s,
    const FormFactorTable::AtomFactorCoefficients& atom_factor_coefficients) {
  s << atom_factor_coefficients.atom_type_ << ' ';
  for (unsigned int i = 0; i < 5; i++) {
    s << atom_factor_coefficients.a_[i] << ' ';
  }
  s << atom_factor_coefficients.c_ << ' ';
  for (unsigned int i = 0; i < 5; i++) {
    s << atom_factor_coefficients.b_[i] << ' ';
  }
  return s << atom_factor_coefficients.excl_vol_ << std::endl;
}

FormFactorTable::FormFactorTable() {
  init_element_form_factor_map();
  init_residue_type_form_factor_map();
}

FormFactorTable::FormFactorTable(const std::string& table_name, double min_q,
                                 double max_q, double delta_q)
    : min_q_(min_q), max_q_(max_q), delta_q_(delta_q) {
  init_element_form_factor_map();
  init_residue_type_form_factor_map();

  // read form factor coefficients from file
  int ffnum = read_form_factor_table(table_name);

  if (ffnum > 0) {  // form factors found in form factor file
    // init zero_form_factors so that they are computed from  the file
    for (int i = 0; i < HEAVY_ATOM_SIZE; i++) {
      zero_form_factors_[i] = 0.0;
      vacuum_zero_form_factors_[i] = 0.0;
      dummy_zero_form_factors_[i] = 0.0;
    }
    // init all the tables
    unsigned int number_of_q_entries =
        algebra::get_rounded((max_q_ - min_q_) / delta_q_) + 1;
    std::vector<double> form_factor_template(number_of_q_entries, 0.0);
    form_factors_ =
        std::vector<std::vector<double>>(HEAVY_ATOM_SIZE, form_factor_template);
    vacuum_form_factors_ =
        std::vector<std::vector<double>>(HEAVY_ATOM_SIZE, form_factor_template);
    dummy_form_factors_ =
        std::vector<std::vector<double>>(HEAVY_ATOM_SIZE, form_factor_template);

    // compute all the form factors
    compute_form_factors_all_atoms();
    compute_form_factors_heavy_atoms();
  }
}

void FormFactorTable::init_element_form_factor_map() {
  element_ff_type_map_ = {
      {"H", H},   {"D", H},   {"HE", He}, {"LI", Li}, {"BE", Be},
      {"B", B},   {"C", C},   {"N", N},   {"O", O},   {"F", F},
      {"NE", Ne}, {"NA", Na}, {"MG", Mg}, {"AL", Al}, {"SI", Si},
      {"P", P},   {"S", S},   {"CL", Cl}, {"AR", Ar}, {"K", K},
      {"CA", Ca}, {"CR", Cr}, {"MN", Mn}, {"FE", Fe}, {"CO", Co},
      {"NI", Ni}, {"CU", Cu}, {"ZN", Zn}, {"SE", Se}, {"BR", Br},
      {"AG", Ag}, {"I", I},   {"IR", Ir}, {"PT", Pt}, {"AU", Au},
      {"HG", Hg}};
}

void FormFactorTable::init_residue_type_form_factor_map() {
  residue_type_form_factor_map_ = {
      {"ALA", FormFactor(9.037, 37.991, 28.954)},
      {"ARG", FormFactor(23.289, 84.972, 61.683)},
      {"ASP", FormFactor(20.165, 58.989, 38.824)},
      {"ASN", FormFactor(19.938, 59.985, 40.047)},
      {"CYS", FormFactor(18.403, 53.991, 35.588)},
      {"GLN", FormFactor(19.006, 67.984, 48.978)},
      {"GLU", FormFactor(19.233, 66.989, 47.755)},
      {"GLY", FormFactor(10.689, 28.992, 18.303)},
      {"HIS", FormFactor(21.235, 78.977, 57.742)},
      {"ILE", FormFactor(6.241, 61.989, 55.748)},
      {"LEU", FormFactor(6.241, 61.989, 55.748)},
      {"LYS", FormFactor(10.963, 70.983, 60.020)},
      {"MET", FormFactor(16.539, 69.989, 53.450)},
      {"PHE", FormFactor(9.206, 77.986, 68.7806)},
      {"PRO", FormFactor(8.613, 51.9897, 43.377)},
      {"SER", FormFactor(13.987, 45.991, 32.004)},
      {"THR", FormFactor(13.055, 53.99, 40.935)},
      {"TYR", FormFactor(14.156, 85.986, 71.83)},
      {"TRP", FormFactor(14.945, 98.979, 84.034)},
      {"VAL", FormFactor(7.173, 53.9896, 46.817)},
      {"POP", FormFactor(45.616, 365.99, 320.41)},
      {"UNK", FormFactor(9.037, 37.991, 28.954)}};
}

int FormFactorTable::read_form_factor_table(const std::string& table_name) {
  std::ifstream s(table_name.c_str());
  if (!s) {
    throw std::runtime_error("Can't find form factor table file " + table_name);
  }

  // init coefficients table
  form_factors_coefficients_ =
      std::vector<AtomFactorCoefficients>(ALL_ATOM_SIZE);

  // skip the comment lines
  char c;
  const int MAX_LENGTH = 1000;
  char line[MAX_LENGTH];
  while (s.get(c)) {
    if (c == '#') {  // if comment line, read the whole line and move on
      s.getline(line, MAX_LENGTH);
    } else {  // return the first character
      s.putback(c);
      break;
    }
  }

  // read the data files
  AtomFactorCoefficients coeff;
  int counter = 0;
  while (!s.eof()) {
    s >> coeff;
    // find FormFactorAtomType
    FormFactorAtomType ff_type = get_form_factor_atom_type(coeff.atom_type_);
    if (ff_type != UNK) {
      form_factors_coefficients_[ff_type] = coeff;
      counter++;
      std::clog << "read_form_factor_table: Atom type found: "
                << coeff.atom_type_ << std::endl;
    } else {
      std::clog << "Atom type is not supported " << coeff.atom_type_
                << std::endl;
    }
  }
  std::clog << counter << " form factors were read from file " << std::endl;
  return counter;
}

void FormFactorTable::show(std::ostream& out, std::string prefix) const {
  for (unsigned int i = 0; i < HEAVY_ATOM_SIZE; i++) {
    out << prefix << " FFATOMTYPE " << i << " zero_ff " << zero_form_factors_[i]
        << " vacuum_ff " << vacuum_zero_form_factors_[i] << " dummy_ff "
        << dummy_zero_form_factors_[i] << std::endl;
  }
}

/*
f(q) = f_atomic(q) - f_solvent(q)
f_atomic(q) = c + SUM [ a_i * EXP( - b_i * (q/4pi)^2 )]
                 i=1,5
f_solvent(q) = rho * v_i * EXP( (- v_i^(2/3) / (4pi)) * q^2 )
*/
void FormFactorTable::compute_form_factors_all_atoms() {
  int number_of_q_entries = (int)std::ceil((max_q_ - min_q_) / delta_q_);

  // iterate over different atom types
  for (unsigned int i = 0; i < ALL_ATOM_SIZE; i++) {
    // form factors for all the q range
    // volr_coeff = - v_i^(2/3) / 4PI
    double volr_coeff =
        -std::pow(form_factors_coefficients_[i].excl_vol_, (2.0 / 3.0)) /
        (16 * PI);

    // iterate over q
    for (int iq = 0; iq < number_of_q_entries; iq++) {
      double q = min_q_ + (double)iq * delta_q_;
      double s = q / (4 * PI);

      // c
      vacuum_form_factors_[i][iq] = form_factors_coefficients_[i].c_;

      // SUM [a_i * EXP( - b_i * (q/4pi)^2 )] Waasmaier and Kirfel (1995)
      for (unsigned int j = 0; j < 5; j++) {
        vacuum_form_factors_[i][iq] +=
            form_factors_coefficients_[i].a_[j] *
            std::exp(-form_factors_coefficients_[i].b_[j] * s * s);
      }
      // subtract solvation: rho * v_i * EXP( (- v_i^(2/3) / (4pi)) * q^2  )
      dummy_form_factors_[i][iq] = rho_ *
                                   form_factors_coefficients_[i].excl_vol_ *
                                   std::exp(volr_coeff * q * q);

      form_factors_[i][iq] =
          vacuum_form_factors_[i][iq] - dummy_form_factors_[i][iq];
    }

    // zero form factors
    zero_form_factors_[i] = form_factors_coefficients_[i].c_;
    for (unsigned int j = 0; j < 5; j++) {
      zero_form_factors_[i] += form_factors_coefficients_[i].a_[j];
    }
    vacuum_zero_form_factors_[i] = zero_form_factors_[i];
    dummy_zero_form_factors_[i] =
        rho_ * form_factors_coefficients_[i].excl_vol_;
    // subtract solvation
    zero_form_factors_[i] -= rho_ * form_factors_coefficients_[i].excl_vol_;
  }
}

void FormFactorTable::compute_form_factors_heavy_atoms() {
  int number_of_q_entries = (int)std::ceil((max_q_ - min_q_) / delta_q_);
  FormFactorAtomType element_type = UNK;
  unsigned int h_num = 0;  // bonded hydrogens number

  for (unsigned int i = ALL_ATOM_SIZE; i < HEAVY_ATOM_SIZE; i++) {
    switch (i) {
      case CH:
        element_type = C;
        h_num = 1;
        break;
      case CH2:
        element_type = C;
        h_num = 2;
        break;
      case CH3:
        element_type = C;
        h_num = 3;
        break;
      case NH:
        element_type = N;
        h_num = 1;
        break;
      case NH2:
        element_type = N;
        h_num = 2;
        break;
      case NH3:
        element_type = N;
        h_num = 3;
        break;
      case OH:
        element_type = O;
        h_num = 1;
        break;
      case OH2:
        element_type = O;
        h_num = 2;
        break;
      case SH:
        element_type = S;
        h_num = 1;
        break;
      default:
        break;
    }

    // full form factors
    for (int iq = 0; iq < number_of_q_entries; iq++) {
      // ff(i) = ff(element) + h_num*ff(hydrogen)
      form_factors_[i][iq] =
          form_factors_[element_type][iq] + h_num * form_factors_[H][iq];
      vacuum_form_factors_[i][iq] = vacuum_form_factors_[element_type][iq] +
                                    h_num * vacuum_form_factors_[H][iq];
      dummy_form_factors_[i][iq] = dummy_form_factors_[element_type][iq] +
                                   h_num * dummy_form_factors_[H][iq];
    }

    // zero form factors
    zero_form_factors_[i] =
        zero_form_factors_[element_type] + h_num * zero_form_factors_[H];
    vacuum_zero_form_factors_[i] = vacuum_zero_form_factors_[element_type] +
                                   h_num * vacuum_zero_form_factors_[H];
    dummy_zero_form_factors_[i] = dummy_zero_form_factors_[element_type] +
                                  h_num * dummy_zero_form_factors_[H];
  }
}


FormFactorTable::FormFactorAtomType FormFactorTable::get_carbon_atom_type(
    const std::string& atom_type, const std::string& residue_type) const {
  const std::string a = normalize_name(atom_type);
  const std::string r = normalize_name(residue_type);
  if (a == "CH") return CH;
  if (a == "CH2") return CH2;
  if (a == "CH3") return CH3;
  if (a == "C") return C;
  if (a == "CA") return r == "GLY" ? CH2 : CH;
  if (a == "CB") {
    if (is_one_of(r, {"ILE", "THR", "VAL"})) return CH;
    return r == "ALA" ? CH3 : CH2;
  }
  if (a == "CG") {
    if (is_one_of(r, {"ASN", "ASP", "HIS", "PHE", "TRP", "TYR"})) return C;
    return r == "LEU" ? CH : CH2;
  }
  if (a == "CG1") return r == "ILE" ? CH2 : CH3;
  if (a == "CG2") return CH3;
  if (a == "CD") return is_one_of(r, {"GLU", "GLN"}) ? C : CH2;
  if (a == "CD1") {
    if (is_one_of(r, {"LEU", "ILE"})) return CH3;
    return is_one_of(r, {"PHE", "TRP", "TYR"}) ? CH : C;
  }
  if (a == "CD2") {
    if (r == "LEU") return CH3;
    return is_one_of(r, {"PHE", "HIS", "TYR"}) ? CH : C;
  }
  if (a == "CE") {
    if (r == "LYS") return CH2;
    if (r == "MET") return CH3;
    return C;
  }
  if (a == "CE1") return is_one_of(r, {"PHE", "HIS", "TYR"}) ? CH : C;
  if (a == "CE2") return is_one_of(r, {"PHE", "TYR"}) ? CH : C;
  if (a == "CZ") return r == "PHE" ? CH : C;
  if (is_one_of(a, {"CZ2", "CZ3", "CE3"})) return r == "TRP" ? CH : C;

  if (a == "C5'") return CH2;
  if (is_one_of(a, {"C1'", "C2'", "C3'", "C4'"})) return CH;
  if (a == "C2") return is_adenine(r) ? CH : C;
  if (a == "C4") return C;
  if (a == "C5") return (is_cytosine(r) || is_uracil(r)) ? CH : C;
  if (a == "C6") {
    return (is_cytosine(r) || is_uracil(r) || is_thymine(r)) ? CH : C;
  }
  if (a == "C7") return CH3;
  if (a == "C8") return CH;
  std::cerr << "Carbon atom not found, using default C form factor for "
            << a << " " << r << std::endl;
  return C;
}

FormFactorTable::FormFactorAtomType FormFactorTable::get_nitrogen_atom_type(
    const std::string& atom_type, const std::string& residue_type) const {
  const std::string a = normalize_name(atom_type);
  const std::string r = normalize_name(residue_type);
  if (a == "N") return r == "PRO" ? N : NH;
  if (a == "ND1") return r == "HIS" ? NH : N;
  if (a == "ND2") return r == "ASN" ? NH2 : N;
  if (is_one_of(a, {"NH1", "NH2"})) return r == "ARG" ? NH2 : N;
  if (a == "NE") return r == "ARG" ? NH : N;
  if (a == "NE1") return r == "TRP" ? NH : N;
  if (a == "NE2") return r == "GLN" ? NH2 : N;
  if (a == "NZ") return r == "LYS" ? NH3 : N;
  if (a == "N1") return is_guanine(r) ? NH : N;
  if (is_one_of(a, {"N2", "N4", "N6"})) return NH2;
  if (a == "N3") return is_uracil(r) ? NH : N;
  if (is_one_of(a, {"N7", "N9"})) return N;
  std::cerr << "Nitrogen atom not found, using default N form factor for "
            << a << " " << r << std::endl;
  return N;
}

FormFactorTable::FormFactorAtomType FormFactorTable::get_oxygen_atom_type(
    const std::string& atom_type, const std::string& residue_type) const {
  const std::string a = normalize_name(atom_type);
  const std::string r = normalize_name(residue_type);
  if (is_one_of(r, {"HOH", "DOD", "WAT"})) return OH2;
  if ((a == "OG" && r == "SER") || (a == "OG1" && r == "THR") ||
      (a == "OH" && r == "TYR") || a == "O2'") {
    return OH;
  }
  return O;
}

FormFactorTable::FormFactorAtomType FormFactorTable::get_sulfur_atom_type(
    const std::string& atom_type, const std::string& residue_type) const {
  return normalize_name(atom_type) == "SG" &&
                 normalize_name(residue_type) == "CYS"
             ? SH
             : S;
}

FormFactorTable::FormFactorAtomType FormFactorTable::get_form_factor_atom_type(
    const Atom& atom, FormFactorType ff_type) const {
  const std::string raw_name = atom.type();
  const std::string atom_name = normalize_name(raw_name);
  const std::string residue_name = normalize_name(atom.residueName());

  std::string element;
  if (atom_name.empty()) {
    element = "";
  } else if (atom_name[0] == 'D' || atom_name[0] == 'H' ||
             std::isdigit(static_cast<unsigned char>(atom_name[0]))) {
    element = "H";
  } else if (!raw_name.empty() && raw_name[0] != ' ' &&
             atom_name.size() >= 2 &&
             element_ff_type_map_.count(atom_name.substr(0, 2)) > 0) {
    element = atom_name.substr(0, 2);
  } else {
    element = atom_name.substr(0, 1);
  }

  FormFactorAtomType ret_type = get_form_factor_atom_type(element);
  if (ff_type == HEAVY_ATOMS) {
    if (ret_type == C) ret_type = get_carbon_atom_type(atom_name, residue_name);
    if (ret_type == N) ret_type = get_nitrogen_atom_type(atom_name, residue_name);
    if (ret_type == O) ret_type = get_oxygen_atom_type(atom_name, residue_name);
    if (ret_type == S) ret_type = get_sulfur_atom_type(atom_name, residue_name);
  }

  if (ret_type >= HEAVY_ATOM_SIZE) {
    std::cerr << "Can't find form factor for atom " << atom_name
              << " using default value of nitrogen" << std::endl;
    ret_type = N;
  }
  return ret_type;
}

FormFactorTable::FormFactorAtomType FormFactorTable::get_form_factor_atom_type(
    const std::string& element) const {
  const auto found = element_ff_type_map_.find(normalize_name(element));
  return found == element_ff_type_map_.end() ? UNK : found->second;
}

double FormFactorTable::get_form_factor(
    const std::string& residue_type) const {
  const auto found =
      residue_type_form_factor_map_.find(normalize_name(residue_type));
  if (found != residue_type_form_factor_map_.end()) return found->second.ff_;
  std::cerr << "Can't find form factor for residue " << residue_type
            << " using default value of ALA" << std::endl;
  return residue_type_form_factor_map_.find("UNK")->second.ff_;
}

double FormFactorTable::get_vacuum_form_factor(
    const std::string& residue_type) const {
  const auto found =
      residue_type_form_factor_map_.find(normalize_name(residue_type));
  if (found != residue_type_form_factor_map_.end())
    return found->second.vacuum_ff_;
  std::cerr << "Can't find form factor for residue " << residue_type
            << " using default value of ALA" << std::endl;
  return residue_type_form_factor_map_.find("UNK")->second.vacuum_ff_;
}

double FormFactorTable::get_dummy_form_factor(
    const std::string& residue_type) const {
  const auto found =
      residue_type_form_factor_map_.find(normalize_name(residue_type));
  if (found != residue_type_form_factor_map_.end())
    return found->second.dummy_ff_;
  std::cerr << "Can't find form factor for residue " << residue_type
            << " using default value of ALA" << std::endl;
  return residue_type_form_factor_map_.find("UNK")->second.dummy_ff_;
}

double FormFactorTable::get_form_factor(const Atom& atom,
                                        FormFactorType ff_type) const {
  if (ff_type == CA_ATOMS || ff_type == RESIDUES)
    return get_form_factor(atom.residueName());
  return zero_form_factors_[get_form_factor_atom_type(atom, ff_type)];
}

double FormFactorTable::get_vacuum_form_factor(
    const Atom& atom, FormFactorType ff_type) const {
  if (ff_type == CA_ATOMS || ff_type == RESIDUES)
    return get_vacuum_form_factor(atom.residueName());
  return vacuum_zero_form_factors_[get_form_factor_atom_type(atom, ff_type)];
}

double FormFactorTable::get_dummy_form_factor(
    const Atom& atom, FormFactorType ff_type) const {
  if (ff_type == CA_ATOMS || ff_type == RESIDUES)
    return get_dummy_form_factor(atom.residueName());
  return dummy_zero_form_factors_[get_form_factor_atom_type(atom, ff_type)];
}

double FormFactorTable::get_radius(const Atom& atom,
                                   FormFactorType ff_type) const {
  static const double c = 3.0 / (4 * PI * rho_);
  return std::cbrt(c * get_dummy_form_factor(atom, ff_type));
}

double FormFactorTable::get_volume(const Atom& atom,
                                   FormFactorType ff_type) const {
  return get_dummy_form_factor(atom, ff_type) / rho_;
}

const std::vector<double>& FormFactorTable::get_form_factors(
    const Atom& atom, FormFactorType ff_type) const {
  return form_factors_[get_form_factor_atom_type(atom, ff_type)];
}

const std::vector<double>& FormFactorTable::get_vacuum_form_factors(
    const Atom& atom, FormFactorType ff_type) const {
  return vacuum_form_factors_[get_form_factor_atom_type(atom, ff_type)];
}

const std::vector<double>& FormFactorTable::get_dummy_form_factors(
    const Atom& atom, FormFactorType ff_type) const {
  return dummy_form_factors_[get_form_factor_atom_type(atom, ff_type)];
}

FormFactorTable* get_default_form_factor_table() {
  static FormFactorTable ff;
  return &ff;
}

}  // namespace foxs
