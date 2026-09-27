#include "SurfacePoint.h"

#include <algorithm>
#include <functional>
#include <ostream>

typedef SurfacePoint::SurfaceType SurfaceType;

SurfacePoint::SurfacePoint(const Vector3& newPos, const Vector3& newNorm,
			   float surface,
                           int atomIndex1,
                           int atomIndex2,
                           int atomIndex3)
  : Vector3(newPos), norm(newNorm), area(surface)
{
  setAtomIndices(atomIndex1, atomIndex2, atomIndex3);
}

void SurfacePoint::setAtomIndices(int atomIndex1,
                                  int atomIndex2,
                                  int atomIndex3)
{
  atoms[0] = atomIndex1;
  atoms[1] = atomIndex2;
  atoms[2] = atomIndex3;
  std::sort(atoms.begin(), atoms.end(), std::greater<int>());
}

bool SurfacePoint::atomIndicesGTE(const unsigned int atomIndex) const noexcept
{
  bool result = (atoms[0] >= (int)atomIndex);
  if (result && (atoms[1] >= 0)) {
    result = (atoms[1] >= (int)atomIndex);
    if (result && (atoms[2] >= 0))
      result = (atoms[2] >= (int)atomIndex);
  }
  return result;
}

bool SurfacePoint::atomIndicesLT(const unsigned int atomIndex) const noexcept
{
  return ((atoms[0] < (int)atomIndex) && (atoms[1] < (int)atomIndex) &&
          (atoms[2] < (int)atomIndex));
}

SurfaceType SurfacePoint::surfaceType() const noexcept
{
  if (atoms[0] == -1)
    return Undef;
  if (atoms[1] == -1)
    return Caps;
  if (atoms[2] == -1)
    return Belts;
  return Pits;
}

bool SurfacePoint::adjacentTo(const SurfacePoint& sp) const
{
  bool result = false;
  switch (surfaceType()) {
  case Caps:
    result =  (atoms[0] == sp.atoms[0] ||
               atoms[0] == sp.atoms[1] ||
               atoms[0] == sp.atoms[2]);
    break;
  case Belts:
    switch (sp.surfaceType()) {
    case Caps:
      result = (atoms[0] == sp.atoms[0] || atoms[1] == sp.atoms[0]);
      break;
    case Belts:
      result = (atoms[0] == sp.atoms[0] && atoms[1] == sp.atoms[1]);
      break;
    case Pits:
      result = ((atoms[0] == sp.atoms[0] && atoms[1] == sp.atoms[1]) ||
                (atoms[0] == sp.atoms[1] && atoms[1] == sp.atoms[2]) ||
                (atoms[0] == sp.atoms[0] && atoms[1] == sp.atoms[2]));
      break;
    case Undef:
      break;
    }
    break;
  case Pits:
    switch (sp.surfaceType()) {
    case Caps:
      result = (atoms[0] == sp.atoms[0] ||
                atoms[1] == sp.atoms[0] ||
                atoms[2] == sp.atoms[0]);
      break;
    case Belts:
      result = ((atoms[0] == sp.atoms[0] && atoms[1] == sp.atoms[1]) ||
                (atoms[1] == sp.atoms[0] && atoms[2] == sp.atoms[1]) ||
                (atoms[0] == sp.atoms[0] && atoms[2] == sp.atoms[1]));
      break;
    case Pits:
      result = (atoms[0] == sp.atoms[0] &&
                atoms[1] == sp.atoms[1] &&
                atoms[2] == sp.atoms[2]);
      break;
    case Undef:
      break;
    }
  case Undef:
    break;
  }
  return result;
}

SurfacePoint& SurfacePoint::operator*=(const Matrix3& lt) noexcept {
  (Vector3&)(*this)=(lt*position());
  orient(lt*normal());
  return *this;
}

SurfacePoint& SurfacePoint::operator*=(const RigidTrans3& rt) noexcept {
  (Vector3&)(*this)=(rt*position());
  orient(rt.rotation()*normal());
  return *this;
}

std::ostream& operator<<(std::ostream& s, const SurfacePoint& sp)
{
  s.setf(std::ios::fixed, std::ios::floatfield);
  s.setf(std::ios::right, std::ios::adjustfield);
  s.precision(3);

  s.width(5);
  s << sp.atoms[0]+1;
  s.width(5);
  s << sp.atoms[1]+1;
  s.width(5);
  s << sp.atoms[2]+1;
  s.width(8);
  s << sp.x();
  s.width(8);
  s << sp.y();
  s.width(8);
  s << sp.z();
  s.width(8);
  s << sp.area;
  s.width(7);
  s << sp.norm.x();
  s.width(7);
  s << sp.norm.y();
  s.width(7);
  s << sp.norm.z();
  s << std::endl;
  return s;
}
