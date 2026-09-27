#ifndef _numerics_h
#define _numerics_h

#include <cmath>

constexpr float MIN_FLOAT = std::numeric_limits<float>::lowest();
constexpr float MAX_FLOAT = std::numeric_limits<float>::max();

constexpr float pi     = 3.14159265358979323846f;
constexpr float sqrt_2 = 1.4142135623730950488f;
constexpr float sqrt_3 = 1.7320508075688772935f;
constexpr float sqrt_6 = 2.4494897427831780982f;

// calcualtes c where c = sqrt(a^2 + b^2).
float pythag(float a, float b);

// recieves a 3x3 matrix a and returns its DVDecompisions in a,w and v such
// that original = a*w*vt
bool svd3d(double a[3][3], double w[3], double v[3][3]);

#endif
