#include "cbem/sources/fields.hpp"
#include <cmath>
#include "cbem/geometry/quadrature.hpp"

namespace cbem {

std::vector<cplx> project_plane_wave(const TriangleMesh& m, cplx k, cplx eps, const Vec3& d, const Vec3& p, int sub) {
    MeshQuadrature q(m, QuadRule::subdivided(sub));
    const cplx se = std::sqrt(eps); const Vec3 dp = cross(d, p);
    const Multivector A = Multivector::vector(p) * se + Multivector::blade(7) * Multivector::vector(dp) * se;
    std::vector<cplx> h(8 * m.size(), cplx(0));
    for (std::size_t t = 0; t < m.size(); ++t) {
        for (int a = 0; a < q.q; ++a) {
            cplx ph = std::exp(cplx(0, 1) * k * dot(d, q.points(t)[a])) * q.weights(t)[a];
            for (int b = 0; b < 8; ++b) h[8 * t + b] += ph * A.c[b];
        }
        for (int b = 0; b < 8; ++b) h[8 * t + b] /= std::sqrt(m.area[t]);
    }
    return h;
}

Multivector far_field(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, const Vec3& xh, int sub) {
    MeshQuadrature q(m, QuadRule::subdivided(sub));
    const Multivector onepx = Multivector::blade(0) + Multivector::vector(xh);
    Multivector F;
    for (std::size_t t = 0; t < m.size(); ++t) {
        Multivector h; for (int b = 0; b < 8; ++b) h.c[b] = hs[8 * t + b] / std::sqrt(m.area[t]);
        Multivector G = onepx * (Multivector::vector(m.normal[t]) * h);
        cplx s = 0; for (int a = 0; a < q.q; ++a) s += std::exp(cplx(0, -1) * k * dot(xh, q.points(t)[a])) * q.weights(t)[a];
        F = F + G * s;
    }
    return F * (cplx(0, -1) * k / (4 * pi));
}

real extinction_cross_section(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, cplx eps, const Vec3& d, const Vec3& p) {
    Multivector F = far_field(m, hs, k, d);
    const cplx se = std::sqrt(eps);
    cplx pe = (p.x * F.c[1] + p.y * F.c[2] + p.z * F.c[4]) / se;
    return std::real(4 * pi / k * std::imag(pe));
}

}  // namespace cbem
