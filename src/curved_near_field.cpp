#include "cbem/sources/curved_near_field.hpp"

#include <algorithm>
#include <chrono>
#include <functional>
#include <stdexcept>

#include "cbem/assembly/curved_entries.hpp"
#include "cbem/geometry/quadrature.hpp"
#include "cbem/hmatrix/aca.hpp"
#include "cbem/hmatrix/cluster_tree.hpp"
#include "cbem/kernel/dirac_kernel.hpp"
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

std::vector<Multivector> scattered_field_curved_hmatrix(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, const std::vector<Vec3>& pts,
                                                        real eps, HStats* stats) {
    if (hs.size() != 24 * m.size()) throw std::invalid_argument("scattered_field_curved_hmatrix: Spur mit 24 Koeffizienten je Element erwartet");
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t N = m.size(), M = pts.size();
    std::vector<Multivector> out(M);
    if (stats) *stats = HStats{};
    if (M == 0 || N == 0) return out;
    EntryParams ep; ep.cache_near = false;
    const CurvedKernelEntries E(m, k, ep);
    const CurvedQuadrature Q(m, QuadRule::dunavant7());                // dieselbe Regel wie die Fernregel von point_integrals
    const int nq = Q.q;
    // Dichten: Z21_(t a, c) = e_c u_(t a) fuer die dichten Bloecke (wie scattered_field_curved); je Quadraturpunkt
    // W_q = n_q sum_a w_q psi_a(y_q) u_(t a) und Z4_(q, 0) = W_q, Z4_(q, c) = e_c W_q fuer den Kern s + v z
    std::vector<Multivector> Z21(3 * N * kCurvedComps), Z4(N * nq * 4);
    CBEM_OMP(omp parallel for schedule(static))
    for (long t = 0; t < static_cast<long>(N); ++t) {
        Multivector u[3];
        for (int a = 0; a < 3; ++a) {
            for (int q = 0; q < 8; ++q) u[a].c[q] = hs[8 * (3 * t + a) + q];
            for (int c = 0; c < kCurvedComps; ++c) Z21[(3 * t + a) * kCurvedComps + c] = Multivector::blade(kCurvedBlade[c]) * u[a];
        }
        const auto& S = E.S(t);
        for (int p = 0; p < nq; ++p) {
            Multivector U;
            for (int a = 0; a < 3; ++a) {
                real psi = 0; for (int j = 0; j < 3; ++j) psi += S[a * 3 + j] * Q.lam[p][j];
                const real g = Q.weights(t)[p] * psi;
                for (int q = 0; q < 8; ++q) U.c[q] += g * u[a].c[q];
            }
            const Multivector W = Multivector::vector(Q.normals(t)[p]) * U;
            Multivector* z = &Z4[(t * nq + p) * 4];
            z[0] = W; z[1] = Multivector::blade(1) * W; z[2] = Multivector::blade(2) * W; z[3] = Multivector::blade(4) * W;
        }
    }
    // Clusterbaeume und Blockzerlegung
    const std::size_t leaf = 64; const real eta = 2.0, sep = 3.0, far = 4.0, max_kdiam = 20.0;
    const ClusterTree rows(pts, std::vector<real>(M, 0.0), leaf), cols(m.flat, leaf);
    std::vector<real> rmax(cols.nodes.size(), 0.0);                    // groesster Elementradius (wie point_integrals) je Cluster
    for (std::size_t b = 0; b < cols.nodes.size(); ++b)
        for (std::size_t i = cols.nodes[b].begin; i < cols.nodes[b].end; ++i) {
            const std::size_t t = cols.perm[i]; const auto v = m.flat.vertices(t); const Vec3& c = m.flat.centroid[t];
            rmax[b] = std::max(rmax[b], std::max(norm(v[0] - c), std::max(norm(v[1] - c), norm(v[2] - c))) + m.bulge[t]);
        }
    struct Block { int r, c; bool lowrank; };
    std::vector<Block> blocks;
    std::function<void(int, int)> part = [&](int t, int s) {
        const ClusterNode& a = rows.nodes[t]; const ClusterNode& b = cols.nodes[s];
        const real dist = box_distance(a, b);
        if (std::min(a.diam, b.diam) <= eta * dist && dist > sep * b.hmax && dist > far * rmax[s] && std::abs(k) * std::max(a.diam, b.diam) <= max_kdiam) {
            blocks.push_back({t, s, true}); return;
        }
        if (a.leaf() && b.leaf()) { blocks.push_back({t, s, false}); return; }
        if (b.leaf() || (!a.leaf() && a.diam >= b.diam)) { for (int x : a.child) part(x, s); }
        else { for (int y : b.child) part(t, y); }
    };
    part(0, 0);
    // Aufbau und sofortige Anwendung je Block; Summen je Thread
    std::size_t nlr = 0, ndense = 0, rsum = 0, rmaxk = 0;
    CBEM_OMP(omp parallel)
    {
        std::vector<Multivector> Fl(M);
        std::vector<Multivector> tmp;
        std::size_t lr_n = 0, d_n = 0, r_s = 0, r_m = 0;
        CBEM_OMP(omp for schedule(dynamic) nowait)
        for (long ib = 0; ib < static_cast<long>(blocks.size()); ++ib) {
            const Block& bl = blocks[ib];
            const std::vector<std::size_t> R = rows.indices(rows.nodes[bl.r]), C = cols.indices(cols.nodes[bl.c]);
            const std::size_t mr = R.size(), nc = C.size();
            if (!bl.lowrank) {
                ++d_n;
                for (std::size_t i = 0; i < mr; ++i) {
                    Multivector& f = Fl[R[i]];
                    for (std::size_t j = 0; j < nc; ++j) {
                        const auto G = E.point_integrals(pts[R[i]], C[j]);
                        for (int a = 0; a < 3; ++a)
                            for (int c = 0; c < kCurvedComps; ++c) {
                                const Multivector& z = Z21[(3 * C[j] + a) * kCurvedComps + c]; const cplx g = G[a][c];
                                for (int q = 0; q < 8; ++q) f.c[q] += g * z.c[q];
                            }
                    }
                }
                continue;
            }
            // Spalten (Komponente c, Quadraturpunkt J = j nq + p): Kern s, v z_x, v z_y, v z_z
            const std::size_t np = nc * nq;
            std::vector<std::vector<cplx>> cache(np);
            RowFn row = [&](std::size_t i, cplx* o) {
                const Vec3& x = pts[R[i]];
                for (std::size_t j = 0; j < nc; ++j) {
                    const Vec3* y = Q.points(C[j]);
                    for (int p = 0; p < nq; ++p) {
                        const Vec3 z = x - y[p]; const KernelValue kv = dirac_kernel_fast(z, k); const std::size_t J = j * nq + p;
                        o[J] = kv.s; o[np + J] = kv.vcoef * z.x; o[2 * np + J] = kv.vcoef * z.y; o[3 * np + J] = kv.vcoef * z.z;
                    }
                }
            };
            ColFn col = [&](std::size_t Jc, cplx* o) {
                const std::size_t c = Jc / np, J = Jc % np;
                std::vector<cplx>& cc = cache[J];
                if (cc.empty()) {                                       // alle vier Komponenten des Quadraturpunkts auf einmal
                    cc.resize(4 * mr); const Vec3 y = Q.points(C[J / nq])[J % nq];
                    for (std::size_t i = 0; i < mr; ++i) {
                        const Vec3 z = pts[R[i]] - y; const KernelValue kv = dirac_kernel_fast(z, k);
                        cc[i] = kv.s; cc[mr + i] = kv.vcoef * z.x; cc[2 * mr + i] = kv.vcoef * z.y; cc[3 * mr + i] = kv.vcoef * z.z;
                    }
                }
                std::copy(cc.begin() + c * mr, cc.begin() + (c + 1) * mr, o);
            };
            const LowRank f = aca_plus(row, col, mr, 4 * np, eps);      // ohne Nachkompression (einmalige Anwendung)
            const std::size_t r = f.rank();
            ++lr_n; r_s += r; r_m = std::max(r_m, r);
            tmp.assign(r, Multivector{});
            for (std::size_t kk = 0; kk < r; ++kk) {                   // tmp = V^T Z4
                const cplx* v = f.V.col(kk); Multivector& tk = tmp[kk];
                for (int c = 0; c < 4; ++c)
                    for (std::size_t J = 0; J < np; ++J) {
                        const cplx vv = v[c * np + J]; const Multivector& z = Z4[(C[J / nq] * nq + J % nq) * 4 + c];
                        for (int q = 0; q < 8; ++q) tk.c[q] += vv * z.c[q];
                    }
            }
            for (std::size_t i = 0; i < mr; ++i) {                     // F += U tmp
                Multivector& fi = Fl[R[i]];
                for (std::size_t kk = 0; kk < r; ++kk) { const cplx uu = f.U(i, kk); for (int q = 0; q < 8; ++q) fi.c[q] += uu * tmp[kk].c[q]; }
            }
        }
        CBEM_OMP(omp critical)
        {
            for (std::size_t i = 0; i < M; ++i) for (int q = 0; q < 8; ++q) out[i].c[q] += Fl[i].c[q];
            nlr += lr_n; ndense += d_n; rsum += r_s; rmaxk = std::max(rmaxk, r_m);
        }
    }
    if (stats) {
        stats->n_dense = ndense; stats->n_lowrank = nlr; stats->max_rank = rmaxk;
        stats->mean_rank = nlr ? double(rsum) / nlr : 0.0;
        stats->seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    return out;                                                         // Vorzeichen +1 wie scattered_field_curved
}

std::vector<NearFieldPoint> exterior_near_field_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b,
                                                       const Medium& m, real omega, const IncidentField& incf, const std::vector<Vec3>& pts,
                                                       const NearFieldOptions& opt) {
    if (h.size() != b.size()) throw std::invalid_argument("exterior_near_field_curved: h und b verschieden lang");
    std::vector<cplx> hs(h.size()); for (std::size_t i = 0; i < h.size(); ++i) hs[i] = h[i] - b[i];
    const bool useH = pts.size() >= opt.hmatrix_min_points;
    auto field = [&](const std::vector<cplx>& hh, cplx kk) {
        return useH ? scattered_field_curved_hmatrix(outer, hh, kk, pts, opt.eps) : scattered_field_curved(outer, hh, kk, pts);
    };
    std::vector<Multivector> Fs;
    if (std::abs(m.chi) > 0) {                                          // chirales Aussenmedium (v0.59): je Helizitaet mit k_pm
        Fs.assign(pts.size(), Multivector{});
        for (int s : {+1, -1}) {
            const auto part = field(helicity_part(hs, s), m.k(omega, s));
            const Multivector P = Multivector::blade(0, 0.5) + Multivector::blade(7, cplx(0, 0.5 * s));   // P_s = (1 + s iI)/2
            for (std::size_t i = 0; i < pts.size(); ++i) Fs[i] = Fs[i] + P * part[i];
        }
    } else Fs = field(hs, m.k(omega));
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
                                                       const Vec3& d0, const CVec3& p, const std::vector<Vec3>& pts, const NearFieldOptions& opt) {
    const Vec3 d = d0 / norm(d0);
    const std::vector<cplx> b = project_plane_wave_curved(outer, plane_wave_incidence(m, omega, d, p).k, m.eps, d, p);   // wie im Loeser
    return exterior_near_field_curved(outer, h, b, m, omega, PlaneWaveField(m, omega, d, p), pts, opt);
}

NearFieldEval make_near_field_eval_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b, const Medium& m,
                                          real omega, std::shared_ptr<const IncidentField> inc, const NearFieldOptions& opt) {
    return [&outer, &h, b, m, omega, inc, opt](const std::vector<Vec3>& pts) { return exterior_near_field_curved(outer, h, b, m, omega, *inc, pts, opt); };
}

NearFieldEval make_plane_wave_eval_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d,
                                          const CVec3& p, const NearFieldOptions& opt) {
    return [&outer, &h, m, omega, d, p, opt](const std::vector<Vec3>& pts) { return exterior_near_field_curved(outer, h, m, omega, d, p, pts, opt); };
}

}  // namespace cbem
