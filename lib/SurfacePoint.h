#ifndef _SurfacePoint_h
#define _SurfacePoint_h

#include <iosfwd>
#include <array>

#include "Vector3.h"
#include "Matrix3.h"
#include "RigidTrans3.h"

/*
CLASS
  SurfacePoint

  A descendant of the Vector3 class. Serves as a base class for surface
  points. Besides their positions holds geometrical information such as nornal
  to the surface and type of surface point (caps, pits and belts).

KEYWORDS
  vector, position, algebra, linear, rigid, matrix, distance, match, vector3,
  normal

AUTHOR
  Zipi Fligelman (zipo@math.tau.ac.il)
  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

GOALS
  Defines a class to hold information about surface points. Besides the
  position inherited from the Vector3 class this class holds the normal to
  surface at the given point and the type of surface point (caps, pits or
  belts). More surface properties may be added including chemical ones.

  The class also holds a parameter which designates which fragment or part of
  the molecule the surface point belongs to. This is useful for hinge-
  flexible molecules.

USAGE
*/

class SurfacePoint : public Vector3
{
public:

  enum SurfaceType { Undef, Caps , Belts , Pits };

  //// Creates a Surface Point with 0-position and a general empty normal
  SurfacePoint() = default;

  //// Creates a point with a given position newPos and a given normal and
  // a given type
  explicit SurfacePoint(const Vector3& newPos, const Vector3& norm,
			float surface = 0.0f,
                        int atomIndex1 = -1,
                        int atomIndex2 = -1,
                        int atomIndex3 = -1);

  //// Returns normal of surface point.
  const Vector3& normal() const noexcept { return norm; }

  //// Returns area of surface represented by surface point and normal
  float surfaceArea() const noexcept { return area; }

  //// Set the atom indices of the surface point. This method allows the class
  // user to set the indices of the atoms from which the surface point was
  // calculated. The user may update them at will.
  void setAtomIndices(int atomIndex1 = -1,
                      int atomIndex2 = -1,
                      int atomIndex3 = -1);

  //// Returns atom indices that participated in creating the surface point.
  // Up to 3 atoms may define a surface point
  // depeneding on the surface point type (Caps = 1, Belts = 2, Pits = 3)
  // indx may take value of 0..2. If less then 3 atoms defined the surface
  // point then atomIndex will return a -1.
  int atomIndex(const unsigned int indx) const noexcept { return atoms[indx]; }

  //// Returns true if all atom indices are greater than or equal to the given
  // atom index. Disregards undefined atom indices (i.e. for caps and belts).
  bool atomIndicesGTE(const unsigned int atomIndex) const noexcept;

  //// Returns true if all atom indices are less than the given
  // atom index. Disregards undefined atom indices (i.e. for caps and belts).
  bool atomIndicesLT(const unsigned int atomIndex) const noexcept;

  //// Returns surface type (caps, pits or belts)
  SurfaceType surfaceType() const noexcept;

  //// Changing a point orientation with a new normal
  void orient(const Vector3 &v) noexcept { norm = v/v.norm(); }

  //// Defines symmetric adjacency relation. This relation stems from a
  // natural triangulation of the protein surface. Every
  // cap point is adjacent to all belts and pits surrounding it, every belt
  // is adjacent to its two pits and two caps neighbors and hence every pit
  // is adjacent to its three caps and three belts in its immediate
  // neighborhood.
  bool adjacentTo(const SurfacePoint& sp) const;

  //// Shift SurfacePoint position and normal using a linear tranfromation
  SurfacePoint& operator*=(const Matrix3&) noexcept;

  //// shift SurfacePoint position and normal using a rigid transformation
  SurfacePoint& operator*=(const RigidTrans3&) noexcept;

  friend std::ostream& operator<<(std::ostream& s, const SurfacePoint& sp);


protected:
  Vector3 norm{};

private:
  // atoms creating the patch filled according to SurfaceType with 1/2/3 atoms
  std::array<int, 3> atoms{{-1, -1, -1}};

  //The area on which the surface point is based
  float area{0.0f};
};

#endif
