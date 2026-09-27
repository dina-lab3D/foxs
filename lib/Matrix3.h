#ifndef _Matrix3_h
#define _Matrix3_h

#include <cmath>
#include <iosfwd>

#include "Vector3.h"

/*
CLASS
  Matrix3

  Defines a 3x3 matrix with a load of constructors, operators, and methods.

KEYWORDS
  linear, algebra, matrix, real, vector, rotation, determinant, adjoint,
  inverse, angle, vector

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)
  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

GOALS
  Matrix3 defines a 3x3 matrix. The class was written so that vector-matrix
  operations in 3D will be performed with highest efficiency. The Matrix3 class
  is implemented without any loops since the dimension is small and known.

USAGE
  Matrix3 hosts a set of constructors, operators and inspection methods. Below
  are just several examples of how one may use the class.
  EXAMPLE
    Matrix3 A(Vector3(1, 2, 3),   // construct matric out of 3 vector3
              Vector3(1, 0, 1),
	      Vector3(4, 5, 6));
    Matrix3 B=A;
    B.transpose();
    C=A+B;                        // + operator. C is symmetric.
    Vector3 v = C.eigenvalues();  // eigenvalues.
    Vector3 u = C*v;              // multiplication by vector.
    Matrix3 D(v);                 // D has v in the diagonal. 0's elsewhere.
    Matrix3 E(v,u);               // E=v*u'
  END
  A special note should be made regarding the rotational matrix constructors
  and methods. No special rotational matrix class was defined for reasons of
  conciseness and efficiency. Matrix3 hosts constructors for rotational
  matrices and inspectors for inspecting angular parameters. These inspection
  methods are probably meaningless when invoked on non-rotational matrices, but
  will not cause run-time errors and their likes.
  EXAMPLE
    Matrix3 A(pi, pi, pi);   // rotational parameters for x-y-z. We get I.
    Matrix3 B(0,  0,  0 );   // again we get I.
    if (A.rotX() == 0) cout << "ok\n";
  END
*/
class Matrix3 {

public:
  //// real is defined using the real defined in Vector3 class. By default
  // this is of type float
  using real = Vector3::real;

  // GROUP: Constructors

  //// Default constructor: 0 matrix.
  constexpr Matrix3() noexcept = default;

  //// Constructs a Matrix3 with d's in the diagonal and 0's else where.
  constexpr explicit Matrix3(real d) noexcept
    : vrow{ Vector3(d,0,0), Vector3(0,d,0), Vector3(0,0,d) } {}

  //// Creates a rotational Matrix3. given rotational values around the
  // x axis (xr), around the y axis (yr) and around the z axis (zr) a
  // linear transformation that rotates space in the order x-y-z is created.
  // Rotation in this order is the default representation for the gamb++
  // library rotational transformations.
  explicit Matrix3(real xr, real yr, real zr) noexcept;

  //// Creates a rotational Matrix3. given rotational values around the
  // z axis (psi), around the x axis (theta) and around the z axis (phi), a
  // linear transformation that rotates space in the order z-x-z is created.
  // Rotation in this order is NOT the default representation for the gamb++
  // library rotational transformations. The int zxz specifies that the
  // rotation order is ZXZ and not default XYZ.
  explicit Matrix3(real psi, real theta, real phi, const int &zxz) noexcept;

  //// Initializes a Matrix3 object by stating all the matrix elements
  // explicitly
  explicit Matrix3(real a11, real a12, real a13,
                             real a21, real a22, real a23,
                             real a31, real a32, real a33) noexcept;

  //// Creates a Matrix3 out of a real type array with first 3 elements
  // for the 1st row, 2nd 3 elements for second row and last 3 elements for
  // the 3rd row.
  constexpr explicit Matrix3(const real trn[9]) noexcept
    : vrow{ Vector3(trn), Vector3(trn + 3), Vector3(trn + 6) } {}

  //// Creates a Matrix3 out of a double type array with first 3 elements
  // for the 1st row, 2nd 3 elements for second row and last 3 elements for
  // the 3rd row.
  explicit Matrix3(const double trn[9]) noexcept
    : vrow{ Vector3(trn), Vector3(trn + 3), Vector3(trn + 6) } {}

  //// Use a 3x3 real 2D array to initialize matrix.
  constexpr explicit Matrix3(const real trn[3][3]) noexcept
    : vrow{ Vector3(trn[0]), Vector3(trn[1]), Vector3(trn[2]) } {}

  //// Use a 3x3 double 2D array to initialize matrix.
  explicit Matrix3(const double trn[3][3]) noexcept
    : vrow{ Vector3(trn[0]), Vector3(trn[1]), Vector3(trn[2]) } {}

  //// Constructs a Matrix3 using the 3 given row vectors. For construction
  // with column vectors use this constructor and then the transpose method.
  constexpr explicit Matrix3(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3) noexcept
    : vrow{v1, v2, v3} {}

  //// Constructs a Matrix3 with Vector3 diag's values in the matrix diagonal
  // and 0's elsewhere.
  constexpr explicit Matrix3(const Vector3 &diag) noexcept
    : vrow{ Vector3(diag[0], 0, 0),
            Vector3(0, diag[1], 0),
            Vector3(0, 0, diag[2]) } {}

  //// Creates the matrix v1*v2' from the two given Vector3s
  constexpr explicit Matrix3(const Vector3 &v1, const Vector3 &v2) noexcept
  : vrow{ v1[0] * v2, v1[1] * v2, v1[2] * v2 } {}

  // GROUP: Inspection methods and properties.

  //// Returns the requested row of the matrix. r must be within 0..2. If
  // not run-time errors may be caused. Notice that a const reference is
  // returned so that this method will perform better than the row method
  // and should be used when const rows are good enough.
  const Vector3& operator[](const unsigned short r) const {
    return vrow[r];
  }

  //// Returns a vector matching the r'th row of the matrix. r must be in
  // the range 1..3 or run time errors may occur.
  Vector3 row(const unsigned short r) const {
    return Vector3(vrow[r-1]);
  }

  //// Returns a vector matching the col'th column of the matrix. col must
  // be in the range of 1..3 or run time errors may occur.
  Vector3 column(const unsigned short col) const {
    return Vector3(vrow[0].x(col), vrow[1].x(col), vrow[2].x(col));
  }

  //// Returns a Vector3 with eigenvalues of Matrix3. If eigenvalues are not
  // real (i.e. complex) the function will return the 0 Vector3. So if the
  // Matrix3 is not a 0-matrix and the function returned a 0-vector the
  // Matrix's eigenvalues are not real.
  Vector3 eigenvalues() const;

  //// Returns the numerical value at position row, col in the 3X3 matrix
  // row, col should be within the range of 1..3 or run time errors may occur.
  constexpr real element(const unsigned short r, const unsigned short col) const noexcept {
    return vrow[r-1][col-1];
  }

  //// Returns the rotation around the x-axis for a rotational matrix assuming
  // rotations are performed in the order of x-y-z.
  // Note that results are meaningless if Matrix3 is not a rotational matrix.
  real rotX() const {
    return atan2(element(3,2), element(3,3));
  }

  //// Returns the rotation around the y-axis for a rotational matrix assuming
  // rotations are performed in the order of x-y-z.
  // Note that results are meaningless if Matrix3 is not a rotational matrix.
  real rotY() const {
    real e11 = element(1,1);
    real e21 = element(2,1);
    return atan2(element(3,1), std::sqrt(e21*e21 + e11*e11));
  }

  //// Returns the rotation around the z-axis for a rotational matrix assuming
  // rotations are performed in the order of x-y-z.
  // Note that results are meaningless if Matrix3 is not a rotational matrix.
  real rotZ() const {
    return atan2(element(2,1), element(1,1));
  }

  //// Returns 3 rotation angles assuming
  // rotations are performed in the order of x-y-z.
  // Note that results are meaningless if Matrix3 is not a rotational matrix.
  Vector3 rot() const {
    return Vector3(rotX(), rotY(), rotZ());
  }

  //// Returns 3 rotation angles assuming
  // rotations are performed in the order of z-x-z.
  // Note that results are meaningless if Matrix3 is not a rotational matrix.
  Vector3 rotZXZ() const {
    real e13 = element(1,3);
    real e23 = element(2,3);
    real psi = atan2(e13, e23);
    real theta = atan2(std::sqrt(e13*e13+e23*e23), element(3,3));
    real phi = atan2(element(3,1), -element(3,2));
    return Vector3(psi, theta, phi);
  }

  //// Returns the determinant of the matrix
  real determinant() const;

  //// Returns the trace of the matrix
  real trace() const {
    return element(1,1) + element(2,2) + element(3,3);
  }

  //// Returns the adjoint matrix.
  Matrix3 adjoint() const;

  // GROUP: Operators

  //// Transposes the given matrix.
  void  transpose(){
    real temp12, temp13, temp23;

    temp12 = element(1,2);
    temp13 = element(1,3);
    temp23 = element(2,3);
    vrow[0] = Vector3(element(1,1), element(2,1), element(3,1));
    vrow[1] = Vector3(temp12,       element(2,2), element(3,2));
    vrow[2] = Vector3(temp13,       temp23,       element(3,3));
  }

  //// Returns the matrix transpose.
  Matrix3 transposeMatrix() const{
    return Matrix3(Vector3(vrow[0][0], vrow[1][0], vrow[2][0]),
		   Vector3(vrow[0][1], vrow[1][1], vrow[2][1]),
		   Vector3(vrow[0][2], vrow[1][2], vrow[2][2]));
  }

  //// In place addition of given matrice's members.
  Matrix3& operator+=(const Matrix3&);

  //// In-place subtraction of give matrice's members.
  Matrix3& operator-=(const Matrix3&);

  //// Addition of two matrices.
  friend Matrix3 operator+(const Matrix3& t1, const Matrix3& t2) {
    return Matrix3(t1.vrow[0]+t2.vrow[0], t1.vrow[1]+t2.vrow[1],
		   t1.vrow[2]+t2.vrow[2]);
  }

  //// Substraction of two matrices
  friend Matrix3 operator-(const Matrix3& t1, const Matrix3& t2) {
    return Matrix3(t1.vrow[0]-t2.vrow[0], t1.vrow[1]-t2.vrow[1],
		   t1.vrow[2]-t2.vrow[2]);
  }

  //// In place multiplication with given matrix.
  Matrix3& operator*=(const Matrix3 &);

  //// Matrix multiplication of two matrices.
  friend Matrix3 operator*(const Matrix3& t1, const Matrix3& t2) {
    Vector3 c1 = t2.column(1);
    Vector3 c2 = t2.column(2);
    Vector3 c3 = t2.column(3);
    return Matrix3(Vector3(t1.vrow[0]*c1, t1.vrow[0]*c2, t1.vrow[0]*c3),
                 Vector3(t1.vrow[1]*c1, t1.vrow[1]*c2, t1.vrow[1]*c3),
		   Vector3(t1.vrow[2]*c1, t1.vrow[2]*c2, t1.vrow[2]*c3));
  }

  //// In-place matrix mutiplication by real scalar.
  Matrix3& operator*=(real m);

  //// In-place matrix division by real scalar.
  Matrix3& operator/=(real m);

  //// Matrix mutiplication by real scalar.
  friend Matrix3 operator*(const Matrix3& t, real m) {
    return Matrix3(m*t.vrow[0], m*t.vrow[1], m*t.vrow[2]);
  }

  //// Matrix inverse. If matrix inverse is un-defined a 0-matrix is returned
  friend Matrix3 operator!(const Matrix3&);

  //// Matrix by vector multiplication.
  friend Vector3 operator*(const Matrix3& t, const Vector3& v) {
    return Vector3(t.vrow[0]*v, t.vrow[1]*v, t.vrow[2]*v);
  }

  //// Default Matrix output. 3 row vectors are output using Vector3's <<
  // operator with new-lines between rows.
  friend std::ostream& operator<<(std::ostream& s, const Matrix3 &t);

private:
  Vector3 vrow[3]{};     /* 3 columns of the matrix are held in 3 vectors */
};

#endif
