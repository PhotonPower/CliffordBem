#include "cbem/sources/curved_near_field.hpp"

#include <algorithm>
#include <stdexcept>

#include "cbem/assembly/curved_entries.hpp"
#include "cbem/geometry/quadrature.hpp"
#include "cbem/problems/curved_problem.hpp"
#include "cbem/sources/chiral_incidence.hpp"

namespace cbem {

std::vector<Vec3> projection_points_curved(const QuadraticMesh& m, int sub) { return CurvedQuadrature(m, QuadRule::subdivided(sub)).x; }

std::vector<cplx> project_samples_curved(const QuadraticMesh& m, const std::vector<CVec3>& Es, const std::vector<CVec3>& Hs, const Medium& med,
                                         int sub) {
    CurvedQuadrature Q(m, QuadRule::subdivided(sub));
    if (Es.size() != Q.x.size() || Hs.size() != Q.x.size()) throw std::invalid_argument("project_samples_curved: je Quadraturpunkt E und H");
    const auto S = curved_psi_matrices(m);
    const cplx se = std::sqrt(med.eps), sm = std::sqrt(med.mu);
    std::vector<cplx> h(24 * m.size(), cplx(0));
    CBEM_OMP(omp parallel for schedule(dynamic, 16))
    for (long t = 0; t < static_cast<long>(m.size()); ++t) {
        for (int q = 0; q < Q.q; ++q) {
            const CVec3& E = Es[t * Q.q + q]; const CVec3& H = Hs[t * Q.q + q];
            const Multivector F = Multivector::vector(E) * se + Multivector::blade(7) * Multivector::vector(H) * sm;
            const real w = Q.weights(t)[q];
            for (int a = 0; a < 3; ++a) {
                real psi = 0; for (int k = 0; k < 3; ++k) psi += S[t][a * 3 + k] * Q.lam[q][k];
                for (int b = 0; b < 8; ++b) h[8 * (3 * t + a) + b] += w * psi * F.c[b];
            }
        }
    }
    return h;
}

std::vector<cplx> project_incident_curved(const QuadraticMesh& m, const IncidentField& inc, const Medium& med, int sub) {
    const std::vector<Vec3> x = projection_points_curved(m, sub);
    std::vector<CVec3> E(x.size()), H(x.size());
    CBEM_OMP(omp parallel for schedule(static))
    for (long i = 0; i < static_cast<long>(x.size()); ++i) inc.eval(x[i], E[i], H[i]);
    return project_samples_curved(m, E, H, med, sub);
}

std::vector<Multivector> scattered_field_curved(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, const std::vector<Vec3>& pts) {
    if (hs.size() != 24 * m.size()) throw std::invalid_argument("scattered_field_curved: Spur mit 24 Koeffizienten je Element erwartet");
    EntryParams ep; ep.cache_near = false;
    const CurvedKernelEntries E(m, k, ep);
    // Z_(t a, c) = e_c u_(t a) mit den Blades der sieben Kernkomponenten (wie CurvedHMatrix)
    const std::size_t NI = 3 * m.size();
    std::vector<Multivector> Z(NI * kCurvedComps);
    for (std::size_t I = 0; I < NI; ++I) {
        Multivector u; for (int q = 0; q < 8; ++q) u.c[q] = hs[8 * I + q];
        for (int c = 0; c < kCurvedComps; ++c) Z[I * kCurvedComps + c] = Multivector::blade(kCurvedBlade[c]) * u;
    }
    std::vector<Multivector> out(pts.size());
    CBEM_OMP(omp parallel for schedule(dynamic, 8))
    for (long i = 0; i < static_cast<long>(pts.size()); ++i) {
        Multivector F;
        for (std::size_t t = 0; t < m.size(); ++t) {
            const auto G = E.point_integrals(pts[i], t);
            for (int a = 0; a < 3; ++a)
                for (int c = 0; c < kCurvedComps; ++c) {
                    const Multivector& z = Z[(3 * t + a) * kCurvedComps + c]; const cplx g = G[a][c];
                    for (int q = 0; q < 8; ++q) F.c[q] += g * z.c[q];
                }
        }
        out[i] = F;                                                     // Vorzeichen +1 wie scattered_field (near_field.cpp)
    }
    return out;
}

std::vector<NearFieldPoint> exterior_near_field_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b,
                                                       const Medium& m, real omega, const IncidentField& incf, const std::vector<Vec3>& pts) {
    if (h.size() != b.size()) throw std::invalid_argument("exterior_near_field_curved: h und b verschieden lang");
    std::vector<cplx> hs(h.size()); for (std::size_t i = 0; i < h.size(); ++i) hs[i] = h[i] - b[i];
    std::vector<Multivector> Fs;
    if (std::abs(m.chi) > 0) {                                          // chirales Aussenmedium (v0.59): je Helizitaet mit k_pm
        Fs.assign(pts.size(), Multivector{});
        for (int s : {+1, -1}) {
            const auto part = scattered_field_curved(outer, helicity_part(hs, s), m.k(omega, s), pts);
            const Multivector P = Multivector::blade(0, 0.5) + Multivector::blade(7, cplx(0, 0.5 * s));   // P_s = (1 + s iI)/2
            for (std::size_t i = 0; i < pts.size(); ++i) Fs[i] = Fs[i] + P * part[i];
        }
    } else Fs = scattered_field_curved(outer, hs, m.k(omega), pts);
    const cplx se = std::sqrt(m.eps), sm = std::sqrt(m.mu);
    const real p2 = incf.reference_E2(), C0 = incf.reference_C();
    const TriangleMesh& flat = outer.flat;
    real hm = 0; for (real hh : flat.hmax) hm += hh; hm /= std::max<std::size_t>(1, flat.hmax.size());
    real bmax = 0; for (real bb : outer.bulge) bmax = std::max(bmax, bb);
    const real close = std::max(0.02 * hm, 4.0 / 3.0 * bmax);         // Woelbung: groesster Abstand der Flaeche von der Sehne
    const std::vector<real> dist = distance_to_surface(flat, pts, close);
    std::vector<NearFieldPoint> out(pts.size());
    CBEM_OMP(omp parallel for schedule(dynamic, 16))
    for (long i = 0; i < static_cast<long>(pts.size()); ++i) {
        NearFieldPoint& r = out[i];
        r.inside = winding_number(flat, pts[i]) > 0.5;
        r.too_close = dist[i] < close;
        CVec3 Ei, Hi; incf.eval(pts[i], Ei, Hi);
        const int VEC[3] = {1, 2, 4}, BIV[3] = {6, 5, 3}; const real BS[3] = {1, -1, 1};
        for (int a = 0; a < 3; ++a) {
            r.E[a] = Fs[i].c[VEC[a]] / se + Ei[a];
            r.H[a] = Fs[i].c[BIV[a]] * BS[a] / sm + Hi[a];
        }
        real e2 = 0; cplx eh = 0;
        for (int a = 0; a < 3; ++a) { e2 += std::norm(r.E[a]); eh += std::conj(r.E[a]) * r.H[a]; }
        r.enhancement = e2 / p2;
        r.chirality = std::imag(eh) / C0;
    }
    return out;
}

std::vector<NearFieldPoint> exterior_near_field_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega,
                                                       const Vec3& d0, const CVec3& p, const std::vector<Vec3>& pts) {
    const Vec3 d = d0 / norm(d0);
    const std::vector<cplx> b = project_plane_wave_curved(outer, plane_wave_incidence(m, omega, d, p).k, m.eps, d, p);   // wie im Loeser
    return exterior_near_field_curved(outer, h, b, m, omega, PlaneWaveField(m, omega, d, p), pts);
}

NearFieldEval make_near_field_eval_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b, const Medium& m,
                                          real omega, std::shared_ptr<const IncidentField> inc) {
    return [&outer, &h, b, m, omega, inc](const std::vector<Vec3>& pts) { return exterior_near_field_curved(outer, h, b, m, omega, *inc, pts); };
}

NearFieldEval make_plane_wave_eval_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d,
                                          const CVec3& p) {
    return [&outer, &h, m, omega, d, p](const std::vector<Vec3>& pts) { return exterior_near_field_curved(outer, h, m, omega, d, p, pts); };
}

}  // namespace cbem
