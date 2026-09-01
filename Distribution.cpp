/*
 *  \file Distribution.cpp \brief computes
 *
 *  distribution classes implementation
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#include "Distribution.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

namespace foxs {

RadialDistributionFunction::RadialDistributionFunction(double bin_size)
    : Distribution<double>(bin_size) {}

RadialDistributionFunction::RadialDistributionFunction(
    const std::string& file_name)
    : Distribution<double>(pr_resolution) {
  read_pr_file(file_name);
}

void RadialDistributionFunction::scale(double c) {
  for (unsigned int i = 0; i < size(); i++) (*this)[i] *= c;
}

void RadialDistributionFunction::add(
    const RadialDistributionFunction& other_rd) {
  for (unsigned int i = 0; i < other_rd.size(); i++) {
    add_to_distribution(other_rd.get_distance_from_index(i), other_rd[i]);
  }
}

double RadialDistributionFunction::R_factor_score(
    const RadialDistributionFunction& model_pr,
    const std::string& file_name) const {
  double sum1 = 0.0, sum2 = 0.0;
  unsigned int distribution_size = std::min(size(), model_pr.size());

  for (unsigned int i = 0; i < distribution_size; i++) {
    sum1 += std::abs((*this)[i] - model_pr[i]);
    sum2 += std::abs((*this)[i]);
  }

  if (file_name.length() > 0) write_fit_file(model_pr, 1.0, file_name);

  return sum1 / sum2;
}

double RadialDistributionFunction::fit(
    const RadialDistributionFunction& model_pr,
    const std::string& file_name) const {
  double max_value = 0.0;
  unsigned int max_index = 0;
  for (unsigned int i = 0; i < size(); i++) {
    if ((*this)[i] > max_value) {
      max_value = (*this)[i];
      max_index = i;
    }
  }
  double c = max_value / model_pr[max_index];

  double sum1 = 0.0, sum2 = 0.0;
  unsigned int distribution_size = std::min(size(), model_pr.size());
  for (unsigned int i = 0; i < distribution_size; i++) {
    sum1 += std::abs((*this)[i] - c * model_pr[i]);
    sum2 += std::abs((*this)[i]);
  }

  if (file_name.length() > 0) write_fit_file(model_pr, c, file_name);

  return sum1 / sum2;
}

// double RadialDistributionFunction::
// chi_score(const RadialDistributionFunction& model_pr) const
// {
//   double chi_square = 0.0;
//   unsigned int distribution_size = std::min(size(), model_pr.size());

//   // compute chi
//   for(unsigned int i = 0; i < distribution_size; i++) {
//     chi_square += square(model_pr[i] - (*this)[i]);
//   }
//   chi_square /= distribution_size;
//   return sqrt(chi_square);
// }

void RadialDistributionFunction::write_fit_file(
    const RadialDistributionFunction& model_pr, double c,
    const std::string& file_name) const {
  std::ofstream out_file(file_name.c_str());
  if (!out_file) {
    throw std::runtime_error("Can't open file " + file_name);
  }

  unsigned int distribution_size = std::min(size(), model_pr.size());
  for (unsigned int i = 0; i < distribution_size; i++) {
    out_file << get_distance_from_index(i) << " " << (*this)[i] << " "
             << c * model_pr[i] << std::endl;
  }
  out_file.close();
}

void RadialDistributionFunction::show(std::ostream& out) const {
  const std::string TITLE_LINE = "Distance distribution";
  out << TITLE_LINE << std::endl;
  for (unsigned int i = 0; i < size(); i++) {
    out << get_distance_from_index(i) << " " << (*this)[i] << std::endl;
  }
}

void RadialDistributionFunction::normalize() {
  // calculate area
  double sum = 0.0;
  for (unsigned int i = 0; i < size(); i++) sum += (*this)[i];

  // normalize
  for (unsigned int i = 0; i < size(); i++) (*this)[i] /= sum;
}

void RadialDistributionFunction::read_pr_file(const std::string& file_name) {
  const std::string TITLE_LINE = "Distance distribution";
  // std::cerr << "start reading pr file " << file_name << std::endl;
  std::ifstream in_file(file_name.c_str());
  if (!in_file) {
    throw std::runtime_error("Can't open file " + file_name);
  }

  double count = 0.0;
  std::string line;
  bool in_distribution = false;
  bool bin_size_set = false;
  while (!in_file.eof()) {
    getline(in_file, line);
    boost::trim(line);  // remove all spaces
    // std::cerr << line << std::endl;
    if (line.substr(0, TITLE_LINE.length()) == TITLE_LINE) {
      in_distribution = true;
      continue;
    }
    if (!in_distribution || line.length() == 0 || line[0] == '\0' ||
        !isdigit(line[0]))
      continue;
    // read distribution line
    std::vector<std::string> split_results;
    boost::split(split_results, line, boost::is_any_of("\t "),
                 boost::token_compress_on);
    if (split_results.size() < 2) continue;
    double r = atof(split_results[0].c_str());
    double pr = atof(split_results[1].c_str());
    if (!bin_size_set && r > 0.0) {
      init(r);
      bin_size_set = true;
      // std::cerr << "read_pr_file: bin_size set to " << r << std::endl;
    }
    add_to_distribution(r, pr);
    count += pr;
  }

  std::clog << "read_pr_file: " << file_name << " size=" << size()
            << " area=" << count << std::endl;
}

}  // namespace foxs
