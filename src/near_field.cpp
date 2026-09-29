#include "cbem/sources/near_field.hpp"
#include <cmath>
#include "cbem/geometry/quadrature.hpp"
#include "cbem/kernel/dirac_kernel.hpp"
#include "cbem/kernel/triangle_integrals.hpp"

namespace cbem {

namespace {
// Vorzeichen der Darstellung des Aussenfeldes durch das Cauchy-Integral: +1, konsistent mit far_field (Fernfeldgrenze)
// und mit der Randwertgrenze F_s -> h_s; gegen die Mie-Loesung in tests/test_near_field.cpp
constexpr real kSign = +1.0;
Multivector projector(int s) { return Multivector::blade(0, 0.5) + Multivector::blade(7, cplx(0, 0.5 * s)); }
}  // namespace

std::vector<Multivector> scattered_field(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, const std::vector<Vec3>& pts) {
    const std::size_t N = m.size();
    MeshQuadrature q(m, QuadRule::dunavant7());
    std::vector<Multivector> G(N);                                       // n_tau u_tau je Dreieck
    std::vector<real> rad(N);
    for (std::size_t t = 0; t < N; ++t) {
        Multivector u; for (int c = 0; c < 8; ++c) u.c[c] = hs[8 * t + c] / std::sqrt(m.area[t]);
        G[t] = Multivector::vector(m.normal[t]) * u;
        const auto tri = m.vertices(t); rad[t] = 0;
        for (const Vec3& v : tri) rad[t] = std::max(rad[t], norm(v - m.centroid[t]));
    }
    const cplx ik = cplx(0, 1) * k;
    std::vector<Multivector> out(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const Vec3& x = pts[i];
        Multivector F;
        for (std::size_t t = 0; t < N; ++t) {
            cplx S = 0; CVec3 V{};
            const real dc = norm(x - m.centroid[t]);
            const Vec3* qp = q.points(t); const real* qw = q.weights(t);
            if (dc > 4 * rad[t]) {                                        // fern: 7-Punkt-Regel mit dem vollen Kern
                for (int a = 0; a < q.q; ++a) {
                    const Vec3 z = x - qp[a]; const KernelValue kv = dirac_kernel_full(z, k);
                    S += qw[a] * kv.s; const cplx c = qw[a] * kv.vcoef; V[0] += c * z.x; V[1] += c * z.y; V[2] += c * z.z;
                }
            } else {                                                      // nah: Phi_0 und 1/r analytisch, Rest numerisch
                Vec3 Ig; real Ii; triangle_integrals(x, m.vertices(t), m.normal[t], Ig, Ii);
                S += -ik * Ii / (4 * pi); V[0] += Ig.x / (4 * pi); V[1] += Ig.y / (4 * pi); V[2] += Ig.z / (4 * pi);
                for (int a = 0; a < q.q; ++a) {
                    const Vec3 z = x - qp[a]; const KernelValue kv = dirac_kernel_remainder(z, k);
                    S += qw[a] * kv.s; const cplx c = qw[a] * kv.vcoef; V[0] += c * z.x; V[1] += c * z.y; V[2] += c * z.z;
                }
            }
            F = F + (Multivector::blade(0, S) + Multivector::vector(V)) * G[t];
        }
        out[i] = F * kSign;
    }
    return out;
}

std::vector<NearFieldPoint> exterior_near_field(const TriangleMesh& outer0, const std::vector<cplx>& h, const Medium& m, real omega,
                                                const Vec3& d0, const CVec3& p, const std::vector<Vec3>& pts) {
    TriangleMesh tmp;                                                   // Geometrie (auch Elementgroessen) bei Bedarf berechnen
    const TriangleMesh& outer = (outer0.hmax.size() == outer0.size() && outer0.centroid.size() == outer0.size()) ? outer0 : (tmp = outer0, tmp.compute_geometry(), tmp);
    const Vec3 d = d0 / norm(d0);
    const PlaneWaveIncidence inc = plane_wave_incidence(m, omega, d, p);
    const std::vector<cplx> b = project_plane_wave(outer, inc.k, m.eps, d, p);
    std::vector<cplx> hs(h.size()); for (std::size_t i = 0; i < h.size(); ++i) hs[i] = h[i] - b[i];
    std::vector<Multivector> Fs;
    if (std::abs(m.chi) > 0) {                                          // je Helizitaet mit k_pm
        Fs.assign(pts.size(), Multivector{});
        for (int s : {+1, -1}) {
            const auto part = scattered_field(outer, helicity_part(hs, s), m.k(omega, s), pts);
            const Multivector P = projector(s);
            for (std::size_t i = 0; i < pts.size(); ++i) Fs[i] = Fs[i] + P * part[i];
        }
    } else Fs = scattered_field(outer, hs, inc.k, pts);
    const cplx se = std::sqrt(m.eps), sm = std::sqrt(m.mu);
    const CVec3 dxp{d.y * p[2] - d.z * p[1], d.z * p[0] - d.x * p[2], d.x * p[1] - d.y * p[0]};
    const real p2 = std::norm(p[0]) + std::norm(p[1]) + std::norm(p[2]);
    const real C0 = std::abs(std::real(se / sm)) * p2;                  // |Im(conj(E0).H0)| der zirkularen Welle, H0 = sqrt(eps/mu) d x E0
    std::vector<NearFieldPoint> out(pts.size());
    real hm = 0; for (real hh : outer.hmax) hm += hh; hm /= std::max<std::size_t>(1, outer.hmax.size());
    const std::vector<real> dist = distance_to_surface(outer, pts, 0.02 * hm);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        NearFieldPoint& r = out[i];
        r.inside = winding_number(outer, pts[i]) > 0.5;
        r.too_close = dist[i] < 0.02 * hm;
        const cplx ph = std::exp(cplx(0, 1) * inc.k * dot(d, pts[i]));
        const int VEC[3] = {1, 2, 4}, BIV[3] = {6, 5, 3}; const real BS[3] = {1, -1, 1};
        for (int a = 0; a < 3; ++a) {
            r.E[a] = Fs[i].c[VEC[a]] / se + p[a] * ph;
            r.H[a] = Fs[i].c[BIV[a]] * BS[a] / sm + (se / sm) * dxp[a] * ph;
        }
        real e2 = 0; cplx eh = 0;
        for (int a = 0; a < 3; ++a) { e2 += std::norm(r.E[a]); eh += std::conj(r.E[a]) * r.H[a]; }
        r.enhancement = e2 / p2;
        r.chirality = std::imag(eh) / C0;
    }
    return out;
}

}  // namespace cbem
