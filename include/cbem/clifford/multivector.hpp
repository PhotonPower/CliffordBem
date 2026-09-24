#pragma once
// Komplexifizierte Clifford-Algebra Cl_3(C).
// Blades als Bitmasken: 1=0, e1=1, e2=2, e12=3, e3=4, e13=5, e23=6, e123=7.
// Die imaginaere Einheit i (komplexe Koeffizienten) ist vom Pseudoskalar I = e123 verschieden.
#include <array>
#include "cbem/core/types.hpp"

namespace cbem {

namespace detail {
constexpr int popcount(unsigned v) { int c = 0; while (v) { c += static_cast<int>(v & 1u); v >>= 1u; } return c; }
constexpr int reorder_sign(unsigned a, unsigned b) {
    a >>= 1u; int s = 0;
    while (a) { s += popcount(a & b); a >>= 1u; }
    return (s & 1) ? -1 : 1;
}
}  // namespace detail

constexpr int blade_grade(int b) { return detail::popcount(static_cast<unsigned>(b)); }
constexpr int blade_sign(int a, int b) { return detail::reorder_sign(static_cast<unsigned>(a), static_cast<unsigned>(b)); }

// Linksmultiplikation als 8x8-Matrix, zeilenweise: M[r*8 + c]
using Mat8 = std::array<cplx, 64>;

struct Multivector {
    std::array<cplx, 8> c{};
    static Multivector blade(int b, cplx v = 1.0) { Multivector m; m.c[b] = v; return m; }
    static Multivector vector(const Vec3& v) { Multivector m; m.c[1] = v.x; m.c[2] = v.y; m.c[4] = v.z; return m; }
    static Multivector vector(const CVec3& v) { Multivector m; m.c[1] = v[0]; m.c[2] = v[1]; m.c[4] = v[2]; return m; }
    Multivector operator+(const Multivector& o) const { Multivector r; for (int i = 0; i < 8; ++i) r.c[i] = c[i] + o.c[i]; return r; }
    Multivector operator-(const Multivector& o) const { Multivector r; for (int i = 0; i < 8; ++i) r.c[i] = c[i] - o.c[i]; return r; }
    Multivector operator*(cplx s) const { Multivector r; for (int i = 0; i < 8; ++i) r.c[i] = c[i] * s; return r; }
    Multivector operator*(const Multivector& o) const {
        Multivector r;
        for (int a = 0; a < 8; ++a) {
            if (c[a] == cplx(0)) continue;
            for (int b = 0; b < 8; ++b) r.c[a ^ b] += static_cast<real>(blade_sign(a, b)) * c[a] * o.c[b];
        }
        return r;
    }
    Multivector reverse() const {
        static constexpr int s[4] = {1, 1, -1, -1};
        Multivector r; for (int b = 0; b < 8; ++b) r.c[b] = static_cast<real>(s[blade_grade(b)]) * c[b]; return r;
    }
    Multivector involute() const {
        Multivector r; for (int b = 0; b < 8; ++b) r.c[b] = (blade_grade(b) % 2 ? -1.0 : 1.0) * c[b]; return r;
    }
    Mat8 left_matrix() const {
        Mat8 M{};
        for (int b = 0; b < 8; ++b) {
            Multivector col = (*this) * blade(b);
            for (int r = 0; r < 8; ++r) M[r * 8 + b] = col.c[r];
        }
        return M;
    }
};

// y = M x
inline void apply(const Mat8& M, const cplx* x, cplx* y) {
    for (int r = 0; r < 8; ++r) { cplx s = 0; for (int c = 0; c < 8; ++c) s += M[r * 8 + c] * x[c]; y[r] = s; }
}

}  // namespace cbem
