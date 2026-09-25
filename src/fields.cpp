#include "cbem/sources/fields.hpp"
#include <cmath>
#include "cbem/geometry/quadrature.hpp"

namespace cbem {

std::vector<cplx> project_plane_wave(const TriangleMesh& m, cplx k, cplx eps, const Vec3& d, const Vec3& p, int sub) {
    return project_plane_wave(m, k, eps, d, CVec3{p.x, p.y, p.z}, sub);
}

CVec3 circular_polarization(const Vec3& d, int s) {
    Vec3 a = std::abs(d.z) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
    Vec3 v = cross(d, a); v = v / norm(v); Vec3 u = cross(v, d);        // (u, v, d) rechtshaendig
    const cplx is(0, s);
    return {u.x + is * v.x, u.y + is * v.y, u.z + is * v.z};
}

std::vector<cplx> project_plane_wave(const TriangleMesh& m, cplx k, cplx eps, const Vec3& d, const CVec3& p, int sub) {
    MeshQuadrature q(m, QuadRule::subdivided(sub));
    const cplx se = std::sqrt(eps);
    const CVec3 dp = {d.y * p[2] - d.z * p[1], d.z * p[0] - d.x * p[2], d.x * p[1] - d.y * p[0]};
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
    return extinction_cross_section(m, hs, k, eps, d, CVec3{p.x, p.y, p.z});
}

real extinction_cross_section(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
    Multivector F = far_field(m, hs, k, d);
    const cplx se = std::sqrt(eps);
    cplx pe = (std::conj(p[0]) * F.c[1] + std::conj(p[1]) * F.c[2] + std::conj(p[2]) * F.c[4]) / se;
    real pn = std::norm(p[0]) + std::norm(p[1]) + std::norm(p[2]);
    return std::real(4 * pi / k * std::imag(pe)) / pn;
}

}  // namespace cbem
