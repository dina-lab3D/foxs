/**
 * \file SolventAccessibleSurface \brief
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */

#include "SolventAccessibleSurface.h"
#include <limits>
#include <cmath>

namespace foxs {

Vector<double> SolventAccessibleSurface::get_solvent_accessibility(
    const std::vector<Sphere>& points, double probe_radius, double density) {
  std::vector<double> radii;
  radii.reserve(points.size());
  for (const Sphere& point : points) radii.push_back(point.radius);
  create_sphere_dots(points, density);
  Vector<double> result;
  result.reserve(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    const double ratio = (points[i].radius + probe_radius) / points[i].radius;
    const algebra::Vector3Ds& dots = get_sphere_dots(points[i].radius);
    std::size_t accessible = 0;
    for (const algebra::Vector3D& dot : dots) {
      const algebra::Vector3D probe_center = points[i].center + ratio * dot;
      bool collides = false;
      for (std::size_t j = 0; j < points.size(); ++j) {
        if (i != j && is_intersecting(probe_center, points[j].center,
                                      probe_radius, points[j].radius)) {
          collides = true;
          break;
        }
      }
      if (!collides) ++accessible;
    }
    result.push_back(dots.empty() ? 0.0
                                  : static_cast<double>(accessible) / dots.size());
  }
  return result;
}

Vector<double> SolventAccessibleSurface::get_solvent_accessibility(
    const Molecule<Atom>& atoms, const std::vector<double>& radii,
    double probe_radius, double density) {
  if (atoms.size() != radii.size())
    throw std::invalid_argument("Each GAMB atom must have a matching SAXS radius");
  Vector<double> result;
  result.reserve(atoms.size());
  create_sphere_dots(radii, density);

  algebra::Vector3Ds coordinates;
  coordinates.reserve(atoms.size());
  for (const Atom& atom : atoms) {
    coordinates.push_back(algebra::Vector3D(atom[0], atom[1], atom[2]));
  }

  for (unsigned int i = 0; i < atoms.size(); ++i) {
    const double atom_radius = radii[i];
    const double ratio = (atom_radius + probe_radius) / atom_radius;
    const algebra::Vector3Ds& sphere_points = get_sphere_dots(atom_radius);
    unsigned int accessible = 0;

    for (const algebra::Vector3D& sphere_point : sphere_points) {
      const algebra::Vector3D probe_center =
          coordinates[i] + ratio * sphere_point;
      bool collides = false;
      for (unsigned int j = 0; j < atoms.size(); ++j) {
        if (i == j) continue;
        if (is_intersecting(probe_center, coordinates[j], probe_radius,
                            radii[j])) {
          collides = true;
          break;
        }
      }
      if (!collides) ++accessible;
    }
    result.push_back(sphere_points.empty()
                         ? 0.0
                         : static_cast<double>(accessible) /
                               sphere_points.size());
  }
  return result;
}

algebra::Vector3Ds SolventAccessibleSurface::create_sphere_dots(double radius,
                                                                double density) {
  algebra::Vector3Ds res;
  double num_equat = 2 * PI * radius * sqrt(density);
  double vert_count = 0.5 * num_equat;

  for (int i = 0; i < vert_count; i++) {
    double phi = (PI * i) / vert_count;
    double z = cos(phi);
    double xy = sin(phi);
    double horz_count = xy * num_equat;
    for (int j = 0; j < horz_count - 1; j++) {
      double teta = (2 * PI * j) / horz_count;
      double x = xy * cos(teta);
      double y = xy * sin(teta);
      res.push_back(algebra::Vector3D(radius * x, radius * y, radius * z));
    }
  }
  return res;
}

void SolventAccessibleSurface::create_sphere_dots(const std::vector<Sphere>& ps,
                                                  double density) {

  if (radii2type_.size() > 0 && density_ != density) {
    radii2type_.clear();
    sphere_dots_.clear();
    density_ = density;
  }
  for (unsigned int i = 0; i < ps.size(); i++) {
    double r = ps[i].radius;
    std::unordered_map<double, int>::const_iterator it = radii2type_.find(r);
    if (it == radii2type_.end()) {
      int type = radii2type_.size();
      radii2type_[r] = type;
      algebra::Vector3Ds dots = create_sphere_dots(r, density);
      sphere_dots_.push_back(dots);
    }
  }
}

void SolventAccessibleSurface::create_sphere_dots(
    const std::vector<double>& radii, double density) {
  if (!radii2type_.empty() && density_ != density) {
    radii2type_.clear();
    sphere_dots_.clear();
  }
  density_ = density;
  for (double radius : radii) {
    if (radii2type_.find(radius) == radii2type_.end()) {
      const int type = static_cast<int>(radii2type_.size());
      radii2type_[radius] = type;
      sphere_dots_.push_back(create_sphere_dots(radius, density));
    }
  }
}

}  // namespace foxs
