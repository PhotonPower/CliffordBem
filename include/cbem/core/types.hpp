#pragma once
// Grundtypen des Clifford-BEM-Rechenkerns.
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

namespace cbem {

using real = double;
using cplx = std::complex<double>;
constexpr real pi = 3.14159265358979323846;

struct Vec3 {
    real x{0}, y{0}, z{0};
    constexpr Vec3() = default;
    constexpr Vec3(real a, real b, real c) : x(a), y(b), z(c) {}
    real& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    real operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(real s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(real s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
};
inline Vec3 operator*(real s, const Vec3& v) { return v * s; }
inline real dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline real norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// Komplexer 3-Vektor (Vektoranteil von Kernintegralen)
using CVec3 = std::array<cplx, 3>;

// Die vier Kernkomponenten eines Dreieckspaars: [s, v_x, v_y, v_z]
// s = Skalarteil, v = Vektorteil von int int Phi_k(x - y) dS_y dS_x.
using KernelComp = std::array<cplx, 4>;

}  // namespace cbem
