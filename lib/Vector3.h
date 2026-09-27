#ifndef _Vector3_h
#define _Vector3_h

#include <array>
#include <iosfwd>
#include <cmath>

/*
CLASS
  Vector3

  Defines a 3D vector class with operators.

KEYWORDS
  linear, algebra, vector, angle, distance, norm, real

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)
  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

GOALS
  Since most of the work on protein structure is done in 3D space, a fast,
  well designed 3D vector class was in place. Knowing the dimension is 3
  allows for a fast, loopless implemetation of Vector3 and Matrix3 to be
  written.

USAGE
  Many 3D vector operations are supplied including vector addition and
  subtraction, dot-product mult., multiplication by scalar, distance and
  angular operators

  A new Vector3 may be constructed using 3 floats or an array of 3 floats.
  EXAMPLE
    Vector3 v(1.0, 2.0, 3.0);
    float x[3];
    Vector3 v(x);
  END

  Coordinates may be accessed using the x(int) method which returns the
  vector's given coordinate (1..3) or using the [] operator (0..2).
  EXAMPLE
    if (v[0] == v.x(1))  cout << "OK\n";
  END

  Addition, subtraction, multiplication and division by scalar, dot-product
  distance, angle.
  EXAMPLE
    w = u+v;                           // addition
    w-= v;                             // in-place subtraction w==u
    if ((w|u) == 0) cout << "OK\n";    // distance operator |
    if (w^u == 0) cout << "OK\n";      // angle operator ^
    w=2*w;                             // multiplicatio by scalar
    float f=w*v;                       // dot-product;
  END
*/
class Vector3 {

public:
  //// Defines the type with which coordinates are defined. Changing this to
  // double will casue the Vector3 and Matrix3 classes to work with double
  // precision coordinates. float was chosen for considerations of speed.
  // double precision is rarely needed since structural data is inexact.
  using real = float;

  // GROUP: Constructors.

  //// Default
  constexpr Vector3() noexcept = default;
  //constexpr Vector3(const Vector3&) noexcept = default;
  //constexpr Vector3& operator=(const Vector3&) noexcept = default;

  //// Initialize using 3 coordinates x,y,z.
  constexpr Vector3(real nx, real ny, real nz) noexcept : c{{nx, ny, nz}} {}

  //// initialize using 3 coordinates in array.
  explicit constexpr Vector3(const real x[3]) noexcept : c{{x[0], x[1], x[2]}} {}

  //// initialize using 3 double precision coordinates in array.
  explicit Vector3(const double x[3]) noexcept
      : c{{static_cast<real>(x[0]), static_cast<real>(x[1]), static_cast<real>(x[2])}} {}

  // GROUP: Inspection.

  //// Returns position (redundant function)
  const Vector3& position() const { return *this; }

  //// return the x, y or z coordinates.
  constexpr real x() const noexcept { return c[0]; }
  constexpr real y() const noexcept { return c[1]; }
  constexpr real z() const noexcept { return c[2]; }

  //// Returns x=x(1), y=x(2) or z=x(3) coordinate. If other values are given
  // result are unexpected.
  real x(const unsigned short coord) const { return c[coord-1]; }

  //// Returns true if the vector is the zero vector
  constexpr bool isZero() const noexcept { return ((c[0] == 0) && (c[1] == 0) && (c[2] == 0)); }

  // Group Update
  void updateX(real x) noexcept { c[0] = x; }
  void updateY(real y) noexcept { c[1] = y; }
  void updateZ(real z) noexcept { c[2] = z; }
  void update(const Vector3& v) noexcept { *this = v; }

  //// Return Vector3's distance from origin.
  real norm() const noexcept { return std::sqrt(norm2()); }

  //// Return Vector3's squared distance from origin
  real norm2() const noexcept { return c[0]*c[0]+c[1]*c[1]+c[2]*c[2]; }

  //// Returns Vector3's distance from Vector3 p.
  real dist(const Vector3& p) const noexcept { return std::sqrt(dist2(p)); }

  //// Returns Vector3's distance squared from Vector3 p.
  real dist2(const Vector3& p) const noexcept {
    real d = c[0]-p.c[0];
    real s = d*d;
    d = c[1]-p.c[1];
    s += d*d;
    d = c[2]-p.c[2];
    return s+d*d;
  }

  //// Access: acceptable indices 0, 1, 2. Other indices may cause run-time errors.
  constexpr const real& operator[](const unsigned short i) const noexcept { return c[i]; }

  //// Access: acceptable indices 0, 1, 2. Other indices may cause run-time errors.
  real& operator[](const unsigned short i) noexcept { return c[i]; }

  ////
  // Returns a perpendicular vector
  // (A zero vector is returned in case of a zero vector). <br>
  // Author: Oranit Dror (oranit@tau.ac.il)
  Vector3 getPerpendicularVector() const {
    if (x() != 0) {
      return Vector3((-y()-z())/x(),1,1);
    } else if (y() != 0) {
      return Vector3(1,(-x()-z())/y(),1);
    } else if (z() != 0) {
      return Vector3(1,1,(-x()-y())/z());
    } else {
      return Vector3(0,0,0);
    }
  }

  ////
  // Returns a unit vector in the same direction
  // (A zero vector is returned in case of a zero vector). <br>
  // Author: Oranit Dror (oranit@tau.ac.il)
  Vector3 getUnitVector() const {
    real n2 = norm2();
    if (n2 == 0) return *this;
    return *this / std::sqrt(n2);
  }

  // GROUP: Operators.

  ////
  // Returns true if the two vectors have the same coordinates
  // TODO: is it correct to check this way for float equality?
  constexpr bool operator==(const Vector3& v) const noexcept {
    return ((x() == v.x()) && (y() == v.y()) && (z() == v.z()));
  }

  //// Returns Vector3 addition of two Vector3s.
  friend constexpr Vector3 operator+(const Vector3& p1, const Vector3& p2) noexcept {
    return Vector3(p1.c[0]+p2.c[0], p1.c[1]+p2.c[1], p1.c[2]+p2.c[2]);
  }

  //// Returns Vector3 subtraction of two Vector3s.
  friend constexpr Vector3 operator-(const Vector3& p1, const Vector3& p2) noexcept {
    return Vector3(p1.c[0]-p2.c[0], p1.c[1]-p2.c[1], p1.c[2]-p2.c[2]);
  }

  //// Returns Vector3 negative.
  friend constexpr Vector3 operator-(const Vector3& p) noexcept {
    return Vector3(-p.c[0], -p.c[1], -p.c[2]);
  }

  //// Returns dot product of two Vector3s.
  friend constexpr real operator*(const Vector3& p1, const Vector3& p2) noexcept {
    return p1.c[0]*p2.c[0]+p1.c[1]*p2.c[1]+p1.c[2]*p2.c[2];
  }

  //// Returns vector multiplied by scalar m.
  friend constexpr Vector3 operator*(const Vector3& p, const real& m) noexcept {
    return Vector3(p.c[0]*m, p.c[1]*m, p.c[2]*m);
  }

  //// Returns vector multiplied by scalar m.
  friend constexpr Vector3 operator*(const real& m, const Vector3& p) noexcept {
    return Vector3(p.c[0]*m, p.c[1]*m, p.c[2]*m);
  }

  //// Returns vector divided by scalar m.
  friend constexpr Vector3 operator/(const Vector3& p, const real& m) noexcept {
    return Vector3(p.c[0]/m, p.c[1]/m, p.c[2]/m);
  }

  //// Returns vector that is vertical to both vectors.
  friend constexpr Vector3 operator&(const Vector3& p1, const Vector3& p2) noexcept {
    return Vector3(p1.c[1]*p2.c[2]-p1.c[2]*p2.c[1],
		   p1.c[2]*p2.c[0]-p1.c[0]*p2.c[2],
		   p1.c[0]*p2.c[1]-p1.c[1]*p2.c[0]);
  }

  //// Returns distance between 2 Vector3s (similar to dist).
  friend real operator|(const Vector3& p1, const Vector3& p2) noexcept {
    return p1.dist(p2);
  }

  //// Returns angle between p1 and p2 in radians.
  friend real operator^(const Vector3& p1, const Vector3& p2) noexcept {
    real angle=p1*p2/(p1.norm()*p2.norm());
    if(angle<-1)
      angle=-1;
    else
      if(angle>1)
	angle=1;
    return acos(angle);
  }

  //// Vector3 addition of another Vector3's coordinates.
  Vector3& operator+=(const Vector3& p) noexcept {
    c[0]+=p.c[0];  c[1]+=p.c[1];  c[2]+=p.c[2];
    return *this;
  }

  //// Vector3 subtraction of another Vector3's coordinates.
  Vector3& operator-=(const Vector3& p) noexcept {
    c[0]-=p.c[0];  c[1]-=p.c[1];  c[2]-=p.c[2];
    return *this;
  }

  //// adds a scalar
  Vector3& operator+=(real m) noexcept {
    c[0]+=m; c[1]+=m; c[2]+=m;
    return *this;
  }

  //// substracts a scalar
  Vector3& operator-=(real m) noexcept {
    c[0]-=m; c[1]-=m; c[2]-=m;
    return *this;
  }

  //// Multiplies coordinates by scalar m.
  Vector3& operator*=(real m) noexcept {
    c[0]*=m;  c[1]*=m;  c[2]*=m;
    return *this;
  }

  //// Divides coordinates by scalar m.
  Vector3& operator/=(real m) noexcept {
    c[0]/=m;  c[1]/=m;  c[2]/=m;
    return *this;
  }

  //// Outputs coordinates delimited by single space.
  friend std::ostream& operator<<(std::ostream& s, const Vector3& v);

  //// Inputs coordinates delimited by single space.
  friend std::istream& operator>>(std::istream& s, Vector3& v);

private:
  std::array<real, 3> c{};
};

#endif
