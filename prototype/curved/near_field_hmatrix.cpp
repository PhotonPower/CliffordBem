// H-Matrix fuer Nahfeldkarten auf gekruemmten Elementen (wie NearFieldOperator im ebenen Pfad): Zeilen = Auswertepunkte,
// Spalten = Elemente; unzulaessige Bloecke direkt mit point_integrals (21 Werte je Punkt und Element), zulaessige per ACA+
// in zwei Varianten:
//   E (Elementebene): Eintraege point_integrals (3 Basisfunktionen x 7 Komponenten je Element), Spalten 21 |E|;
//   P (Punktebene):   Eintraege des Dirac-Kerns (s, v z) an den 7 Quadraturpunkten je Element, Spalten 4 x 7 |E|; die Dichte
//                     wird je Quadraturpunkt zu W_q = n_q sum_a w_q psi_a(y_q) u_a zusammengefasst (F = sum_q (s + v z) W_q,
//                     genau die 7-Punkt-Regel von point_integrals fuer ferne Paare).
// Zulaessig: min(diam) <= eta dist, dist > 3 h_max, dist > 4 r_max (Fernregel von point_integrals), |k| diam <= 20.
// Fehler gegen scattered_field_curved (direkte Summation), bezogen auf max |F_s|. Goldkugel in Wasser, Karte in y = 0.
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -DCBEM_USE_OPENMP -fcx-fortran-rules -Iinclude prototype/curved/near_field_hmatrix.cpp build/libcbem.a -o near_field_hmatrix
//   ./near_field_hmatrix [n] [Punkte]   (Ikosaederkugel mit 20 n^2 Elementen)
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "cbem/assembly/curved_entries.hpp"
#include "cbem/geometry/quadrature.hpp"
#include "cbem/hmatrix/aca.hpp"
#include "cbem/hmatrix/cluster_tree.hpp"
#include "cbem/kernel/dirac_kernel.hpp"
#include "cbem/problems/curved_problem.hpp"
#include "cbem/sources/curved_near_field.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

// Kopie von aca_plus (src/aca.cpp) mit Betragsquadraten statt std::abs (hypot) bei Pivotsuche und Referenzpruefung
static LowRank aca_plus_norm(const RowFn& row, const ColFn& col, std::size_t m, std::size_t n, real eps, std::size_t rmax = 300) {
    std::vector<std::vector<cplx>> Us, Vs;
    std::vector<char> urow(m, 0), ucol(n, 0);
    const std::size_t kmax = std::min({rmax, m, n});
    // Rest einer Zeile bzw. Spalte gegen die bisherigen Kreuze
    auto res_row = [&](std::size_t i, std::vector<cplx>& r) { r.resize(n); row(i, r.data()); for (std::size_t k = 0; k < Us.size(); ++k) { const cplx a = Us[k][i]; for (std::size_t c = 0; c < n; ++c) r[c] -= a * Vs[k][c]; } };
    auto res_col = [&](std::size_t j, std::vector<cplx>& u) { u.resize(m); col(j, u.data()); for (std::size_t k = 0; k < Us.size(); ++k) { const cplx a = Vs[k][j]; for (std::size_t q = 0; q < m; ++q) u[q] -= Us[k][q] * a; } };
    auto nrm2 = [](const std::vector<cplx>& x) { real s = 0; for (auto v : x) s += std::norm(v); return s; };
    // Referenzen: gleichmaessig verteilte, noch nicht verwendete Indizes (deterministisch), deren Rest nicht verschwindet
    std::size_t nref_r = 0, nref_c = 0;
    auto next_ref = [](std::size_t len, std::size_t& counter, const std::vector<char>& used, std::size_t avoid) {
        for (std::size_t t = 0; t < len; ++t) {
            const std::size_t idx = (counter * 7919 + len / 2) % len; ++counter;
            if (!used[idx] && idx != avoid) return idx;
        }
        for (std::size_t q = 0; q < len; ++q) if (!used[q] && q != avoid) return q;
        return len;
    };
    std::size_t ir = next_ref(m, nref_r, urow, m), jr = next_ref(n, nref_c, ucol, n);
    std::vector<cplx> rref, cref, u, v;
    if (ir < m) res_row(ir, rref);
    if (jr < n) res_col(jr, cref);
    // Masstab fuer "verschwunden": relativ zum groessten bisher gesehenen Eintrag (Rundungsrauschen nach exakter Erfassung ist
    // nicht null; ein absoluter Grenzwert liesse die Referenz auf dem Rauschen pivotieren)
    real scale = 0;
    auto maxabs = [](const std::vector<cplx>& x) { real a = 0; for (auto v0 : x) a = std::max(a, std::norm(v0)); return std::sqrt(a); };
    scale = std::max(maxabs(rref), maxabs(cref));
    auto vanished = [&](const std::vector<cplx>& x) { return maxabs(x) <= 1e-12 * scale; };
    real S2 = 0;
    for (std::size_t it = 0; it < kmax; ++it) {
        // Referenzen ersetzen, deren Rest verschwunden ist oder die schon Pivot waren
        for (int tries = 0; ir < m && (urow[ir] || vanished(rref)); ++tries) {
            if (tries == 8) { ir = m; break; }                                  // keine Referenzzeile mit Rest mehr
            ir = next_ref(m, nref_r, urow, ir); if (ir < m) { res_row(ir, rref); scale = std::max(scale, maxabs(rref)); }
        }
        for (int tries = 0; jr < n && (ucol[jr] || vanished(cref)); ++tries) {
            if (tries == 8) { jr = n; break; }
            jr = next_ref(n, nref_c, ucol, jr); if (jr < n) { res_col(jr, cref); scale = std::max(scale, maxabs(cref)); }
        }
        if (ir == m && jr == n) break;                                          // beide Referenzen ohne Rest: Block erfasst
        std::size_t js = n, is = m; real a = -1, b = -1;
        if (ir < m) for (std::size_t c = 0; c < n; ++c) if (!ucol[c] && std::norm(rref[c]) > a) { a = std::norm(rref[c]); js = c; }
        if (jr < n) for (std::size_t q = 0; q < m; ++q) if (!urow[q] && std::norm(cref[q]) > b) { b = std::norm(cref[q]); is = q; }
        if (a > 0) a = std::sqrt(a);
        if (b > 0) b = std::sqrt(b);
        if (std::max(a, b) <= 1e-12 * scale) break;
        std::size_t i, j; cplx piv;
        if (a >= b) {                                                           // Spalte aus der Referenzzeile
            j = js; res_col(j, u);
            i = m; real bu = -1; for (std::size_t q = 0; q < m; ++q) if (!urow[q] && std::norm(u[q]) > bu) { bu = std::norm(u[q]); i = q; }
            if (bu > 0) bu = std::sqrt(bu);
            if (i == m || bu <= 1e-12 * scale) { ucol[j] = 1; continue; }
            res_row(i, v);
        } else {                                                                // Zeile aus der Referenzspalte
            i = is; res_row(i, v);
            j = n; real bv = -1; for (std::size_t c = 0; c < n; ++c) if (!ucol[c] && std::norm(v[c]) > bv) { bv = std::norm(v[c]); j = c; }
            if (bv > 0) bv = std::sqrt(bv);
            if (j == n || bv <= 1e-12 * scale) { urow[i] = 1; continue; }
            res_col(j, u);
        }
        piv = v[j];
        for (auto& x : v) x /= piv;                                             // u v^T mit v(j) = 1
        const real nu = std::sqrt(nrm2(u)), nv = std::sqrt(nrm2(v));
        real cross = 0;
        for (std::size_t k = 0; k < Us.size(); ++k) {
            cplx pa = 0, pb = 0;
            for (std::size_t q = 0; q < m; ++q) pa += std::conj(Us[k][q]) * u[q];
            for (std::size_t c = 0; c < n; ++c) pb += std::conj(Vs[k][c]) * v[c];
            cross += 2 * std::real(pa * pb);
        }
        S2 += nu * nu * nv * nv + cross;
        Us.push_back(u); Vs.push_back(v); urow[i] = 1; ucol[j] = 1;
        if (ir < m) { const cplx c0 = u[ir]; for (std::size_t c = 0; c < n; ++c) rref[c] -= c0 * v[c]; }   // Referenzreste nachfuehren
        if (jr < n) { const cplx c0 = v[jr]; for (std::size_t q = 0; q < m; ++q) cref[q] -= u[q] * c0; }
        const real S = std::sqrt(std::abs(S2));
        const real er = ir < m ? std::sqrt(nrm2(rref) * m) : 0, ec = jr < n ? std::sqrt(nrm2(cref) * n) : 0;   // auf den Block hochgerechnet
        if (nu * nv <= eps * S && er <= eps * S && ec <= eps * S) break;
    }
    LowRank lr; lr.U = Matrix(m, Us.size()); lr.V = Matrix(n, Vs.size());
    for (std::size_t k = 0; k < Us.size(); ++k) { std::copy(Us[k].begin(), Us[k].end(), lr.U.col(k)); std::copy(Vs[k].begin(), Vs[k].end(), lr.V.col(k)); }
    return lr;
}

struct Result { double build, apply, err, rank, mb, dense_frac, kfrac, taca, trec, tsetup, tdense; };

static Result run(const QuadraticMesh& q, cplx k, const std::vector<cplx>& hs, const std::vector<Vec3>& pts, const std::vector<Multivector>& ref,
                  bool point_level, real eps, real eta, std::size_t leaf, bool fast, bool recomp) {
    auto t0 = clk::now();
    const CurvedKernelEntries E(q, k, [] { EntryParams ep; ep.cache_near = false; return ep; }());
    const CurvedQuadrature Q(q, QuadRule::dunavant7());
    const int nq = Q.q;
    const auto S = curved_psi_matrices(q);
    const std::size_t N = q.size(), M = pts.size();
    const ClusterTree rows(pts, std::vector<real>(M, 0.0), leaf), cols(q.flat, leaf);
    // groesster Elementradius je Spaltencluster (Fernregel von point_integrals: |x - c| > 4 rad)
    std::vector<real> rad(N), rmax(cols.nodes.size(), 0.0);
    for (std::size_t t = 0; t < N; ++t) {
        const auto v = q.flat.vertices(t); const Vec3& c = q.flat.centroid[t];
        rad[t] = std::max(norm(v[0] - c), std::max(norm(v[1] - c), norm(v[2] - c))) + q.bulge[t];
    }
    for (std::size_t b = 0; b < cols.nodes.size(); ++b)
        for (std::size_t i = cols.nodes[b].begin; i < cols.nodes[b].end; ++i) rmax[b] = std::max(rmax[b], rad[cols.perm[i]]);
    std::vector<std::pair<int, int>> adm, inadm;
    std::function<void(int, int)> part = [&](int t, int s) {
        const ClusterNode& a = rows.nodes[t]; const ClusterNode& b = cols.nodes[s];
        const real dist = box_distance(a, b);
        const bool ok = std::min(a.diam, b.diam) <= eta * dist && dist > 3.0 * b.hmax && dist > 4.0 * rmax[s] && std::abs(k) * std::max(a.diam, b.diam) <= 20.0;
        if (ok) { adm.emplace_back(t, s); return; }
        if (a.leaf() && b.leaf()) { inadm.emplace_back(t, s); return; }
        if (b.leaf() || (!a.leaf() && a.diam >= b.diam)) { for (int x : a.child) part(x, s); } else { for (int y : b.child) part(t, y); }
    };
    part(0, 0);
    const double tsetup = since(t0); const auto td = clk::now();
    struct Dense { std::vector<std::size_t> R, C; std::vector<std::array<CurvedComp, 3>> G; };
    struct LR { std::vector<std::size_t> R, C; LowRank f; };
    std::vector<Dense> dense(inadm.size()); std::vector<LR> lr(adm.size());
    CBEM_OMP(omp parallel for schedule(dynamic))
    for (long b = 0; b < static_cast<long>(inadm.size()); ++b) {
        Dense& D = dense[b];
        D.R = rows.indices(rows.nodes[inadm[b].first]); D.C = cols.indices(cols.nodes[inadm[b].second]);
        D.G.resize(D.R.size() * D.C.size());
        for (std::size_t a = 0; a < D.R.size(); ++a) for (std::size_t c = 0; c < D.C.size(); ++c) D.G[a * D.C.size() + c] = E.point_integrals(pts[D.R[a]], D.C[c]);
    }
    const double tdense = since(td); const auto tl = clk::now();
    std::vector<double> tblk(adm.size());
    double kev = 0, taca = 0, trec = 0;                                  // Kernauswertungen (Paare x 7 bzw. Punktpaare), Zeiten
    CBEM_OMP(omp parallel for schedule(dynamic) reduction(+ : kev, taca, trec))
    for (long b = 0; b < static_cast<long>(adm.size()); ++b) {
        LR& B = lr[b];
        double ev = 0; const auto ta = clk::now();
        B.R = rows.indices(rows.nodes[adm[b].first]); B.C = cols.indices(cols.nodes[adm[b].second]);
        const std::size_t m = B.R.size(), n = B.C.size();
        if (!point_level) {
            std::vector<std::vector<cplx>> cc(n);
            RowFn row = [&](std::size_t i, cplx* out) {
                ev += 7.0 * n;
                for (std::size_t j = 0; j < n; ++j) { const auto G = E.point_integrals(pts[B.R[i]], B.C[j]);
                    for (int a = 0; a < 3; ++a) for (int c = 0; c < kCurvedComps; ++c) out[(a * kCurvedComps + c) * n + j] = G[a][c]; }
            };
            ColFn col = [&](std::size_t J, cplx* out) {
                const std::size_t ac = J / n, j = J % n;
                if (cc[j].empty()) {
                    cc[j].resize(21 * m); ev += 7.0 * m;
                    for (std::size_t i = 0; i < m; ++i) { const auto G = E.point_integrals(pts[B.R[i]], B.C[j]);
                        for (int a = 0; a < 3; ++a) for (int c = 0; c < kCurvedComps; ++c) cc[j][(a * kCurvedComps + c) * m + i] = G[a][c]; }
                }
                std::copy(cc[j].begin() + ac * m, cc[j].begin() + (ac + 1) * m, out);
            };
            B.f = fast ? aca_plus_norm(row, col, m, 21 * n, eps) : aca_plus(row, col, m, 21 * n, eps);
        } else {
            const std::size_t np = n * nq;
            std::vector<std::vector<cplx>> cc(np);
            RowFn row = [&](std::size_t i, cplx* out) {
                const Vec3& x = pts[B.R[i]]; ev += double(np);
                for (std::size_t j = 0; j < n; ++j) { const Vec3* y = Q.points(B.C[j]);
                    for (int p = 0; p < nq; ++p) { const Vec3 z = x - y[p]; const KernelValue kv = dirac_kernel_fast(z, k); const std::size_t J = j * nq + p;
                        out[J] = kv.s; out[np + J] = kv.vcoef * z.x; out[2 * np + J] = kv.vcoef * z.y; out[3 * np + J] = kv.vcoef * z.z; } }
            };
            ColFn col = [&](std::size_t J, cplx* out) {
                const std::size_t c = J / np, jp = J % np;
                if (cc[jp].empty()) {
                    cc[jp].resize(4 * m); ev += double(m); const Vec3 y = Q.points(B.C[jp / nq])[jp % nq];
                    for (std::size_t i = 0; i < m; ++i) { const Vec3 z = pts[B.R[i]] - y; const KernelValue kv = dirac_kernel_fast(z, k);
                        cc[jp][i] = kv.s; cc[jp][m + i] = kv.vcoef * z.x; cc[jp][2 * m + i] = kv.vcoef * z.y; cc[jp][3 * m + i] = kv.vcoef * z.z; }
                }
                std::copy(cc[jp].begin() + c * m, cc[jp].begin() + (c + 1) * m, out);
            };
            B.f = fast ? aca_plus_norm(row, col, m, 4 * np, eps) : aca_plus(row, col, m, 4 * np, eps);
        }
        const auto tr = clk::now();
        if (recomp) recompress(B.f, eps);
        tblk[b] = since(ta); kev += ev; taca += std::chrono::duration<double>(tr - ta).count(); trec += since(tr);
    }
    {
        const double tlr = since(tl);
        std::size_t bmax = 0; double tsum = 0; for (std::size_t b = 0; b < adm.size(); ++b) { tsum += tblk[b]; if (tblk[b] > tblk[bmax]) bmax = b; }
        std::printf("      %zu zulaessige Bloecke, Wandzeit %.2f s, Summe %.2f s, groesster %.2f s (%zu x %zu, Rang %zu)%c", adm.size(), tlr, tsum,
                    tblk[bmax], lr[bmax].R.size(), lr[bmax].C.size(), lr[bmax].f.rank(), 10);
    }
    Result res{};
    res.build = since(t0);
    res.kfrac = kev / (7.0 * M * N); res.taca = taca; res.trec = trec; res.tsetup = tsetup; res.tdense = tdense;
    // Anwendung
    t0 = clk::now();
    std::vector<Multivector> Z21(3 * N * kCurvedComps), Z4(N * nq * 4);
    for (std::size_t t = 0; t < N; ++t) {
        Multivector u[3];
        for (int a = 0; a < 3; ++a) { for (int c = 0; c < 8; ++c) u[a].c[c] = hs[8 * (3 * t + a) + c];
            for (int c = 0; c < kCurvedComps; ++c) Z21[(3 * t + a) * kCurvedComps + c] = Multivector::blade(kCurvedBlade[c]) * u[a]; }
        for (int p = 0; p < nq; ++p) {
            Multivector U;
            for (int a = 0; a < 3; ++a) { real psi = 0; for (int kk = 0; kk < 3; ++kk) psi += S[t][a * 3 + kk] * Q.lam[p][kk];
                const real g = Q.weights(t)[p] * psi; for (int c = 0; c < 8; ++c) U.c[c] += g * u[a].c[c]; }
            const Multivector W = Multivector::vector(Q.normals(t)[p]) * U;
            const std::size_t J = t * nq + p;
            Z4[J * 4] = W;
            for (int c = 0; c < 3; ++c) { Vec3 e(0, 0, 0); e[c] = 1; Z4[J * 4 + 1 + c] = Multivector::vector(e) * W; }
        }
    }
    std::vector<Multivector> F(M);
    const long nd = static_cast<long>(dense.size()), nb = nd + static_cast<long>(lr.size());
    CBEM_OMP(omp parallel)
    {
        std::vector<Multivector> Fl(M);
        std::vector<Multivector> tmp;
        CBEM_OMP(omp for schedule(dynamic) nowait)
        for (long b = 0; b < nb; ++b) {
            if (b < nd) {
                const Dense& D = dense[b]; const std::size_t n = D.C.size();
                for (std::size_t a = 0; a < D.R.size(); ++a) { Multivector& f = Fl[D.R[a]];
                    for (std::size_t j = 0; j < n; ++j) { const auto& G = D.G[a * n + j];
                        for (int s = 0; s < 3; ++s) for (int c = 0; c < kCurvedComps; ++c) { const Multivector& z = Z21[(3 * D.C[j] + s) * kCurvedComps + c];
                            for (int r = 0; r < 8; ++r) f.c[r] += G[s][c] * z.c[r]; } } }
            } else {
                const LR& B = lr[b - nd]; const std::size_t m = B.R.size(), n = B.C.size(), r = B.f.rank();
                tmp.assign(r, Multivector{});
                for (std::size_t kk = 0; kk < r; ++kk) {
                    const cplx* v = B.f.V.col(kk); Multivector& t = tmp[kk];
                    if (!point_level) {
                        for (int ac = 0; ac < 21; ++ac) for (std::size_t j = 0; j < n; ++j) { const cplx vv = v[ac * n + j];
                            const Multivector& z = Z21[(3 * B.C[j] + ac / kCurvedComps) * kCurvedComps + ac % kCurvedComps]; for (int s = 0; s < 8; ++s) t.c[s] += vv * z.c[s]; }
                    } else {
                        const std::size_t np = n * nq;
                        for (int c = 0; c < 4; ++c) for (std::size_t J = 0; J < np; ++J) { const cplx vv = v[c * np + J];
                            const Multivector& z = Z4[(B.C[J / nq] * nq + J % nq) * 4 + c]; for (int s = 0; s < 8; ++s) t.c[s] += vv * z.c[s]; }
                    }
                }
                for (std::size_t a = 0; a < m; ++a) { Multivector& f = Fl[B.R[a]];
                    for (std::size_t kk = 0; kk < r; ++kk) { const cplx uu = B.f.U(a, kk); for (int s = 0; s < 8; ++s) f.c[s] += uu * tmp[kk].c[s]; } }
            }
        }
        CBEM_OMP(omp critical)
        for (std::size_t i = 0; i < M; ++i) for (int s = 0; s < 8; ++s) F[i].c[s] += Fl[i].c[s];
    }
    res.apply = since(t0);
    real emax = 0, fmax = 0;
    for (std::size_t i = 0; i < M; ++i) { real e = 0, f = 0; for (int s = 0; s < 8; ++s) { e += std::norm(F[i].c[s] - ref[i].c[s]); f += std::norm(ref[i].c[s]); }
        emax = std::max(emax, std::sqrt(e)); fmax = std::max(fmax, std::sqrt(f)); }
    res.err = emax / fmax;
    std::size_t rs = 0, ent = 0, dpairs = 0;
    for (auto& B : lr) { rs += B.f.rank(); ent += B.f.storage(); }
    for (auto& D : dense) { ent += 21 * D.G.size(); dpairs += D.G.size(); }
    res.rank = lr.empty() ? 0 : double(rs) / lr.size();
    res.mb = 16.0 * ent / 1e6;
    res.dense_frac = double(dpairs) / (double(M) * N);
    return res;
}

int main(int argc, char** argv) {
    const int n = argc > 1 ? std::atoi(argv[1]) : 8;
    const int Mreq = argc > 2 ? std::atoi(argv[2]) : 10000;
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, water{1.7689, 1.0, 0.0};
    const real om = 0.5; const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    SolveOptions so; so.tol = 1e-8;
    CurvedScatteringProblem C({quadratic_icosphere(n)}, {gold}, om, water);
    const auto rc = C.solve_plane_wave(d, px, so);
    const cplx k = water.k(om);
    const auto b = project_plane_wave_curved(C.mesh(), plane_wave_incidence(water, om, d, px).k, water.eps, d, px);
    std::vector<cplx> hs(b.size()); for (std::size_t i = 0; i < b.size(); ++i) hs[i] = rc.h[i] - b[i];
    std::vector<Vec3> pts;
    const int g = static_cast<int>(std::ceil(std::sqrt(16.0 * Mreq / (16 - pi))));
    for (int i = 0; i < g && static_cast<int>(pts.size()) < Mreq; ++i)
        for (int j = 0; j < g && static_cast<int>(pts.size()) < Mreq; ++j) {
            const Vec3 x(-2 + 4.0 * (i + 0.5) / g, 0, -2 + 4.0 * (j + 0.5) / g);
            if (norm(x) > 1.02) pts.push_back(x);
        }
    auto t0 = clk::now();
    const auto ref = scattered_field_curved(C.mesh(), hs, k, pts);
    std::printf("%zu Elemente, %zu Punkte: direkt %.2f s\n", C.mesh().size(), pts.size(), since(t0));
    struct Cfg { bool pl; real eta, eps; std::size_t leaf; bool fast, recomp; };
    for (const Cfg& c : {Cfg{false, 2, 1e-4, 32, false, true}, Cfg{false, 2, 1e-4, 32, true, true}, Cfg{true, 2, 1e-4, 64, false, true},
                         Cfg{true, 2, 1e-4, 64, true, true}, Cfg{true, 2, 1e-4, 64, true, false}, Cfg{false, 2, 1e-4, 32, true, false},
                         Cfg{true, 2, 1e-6, 64, true, true}, Cfg{true, 2, 1e-6, 64, true, false}}) {
        const Result r = run(C.mesh(), k, hs, pts, ref, c.pl, c.eps, c.eta, c.leaf, c.fast, c.recomp);
        std::printf("  %s%s%s eta %.0f eps %.0e Blatt %zu: Aufbau %.2f s, Anwendung %.3f s, Fehler %.1e, Rang %.1f, %.0f MB, dicht %.1f %%" "%c",
                    c.pl ? "P" : "E", c.fast ? " norm" : " abs", c.recomp ? "" : " ohne Nachkompr.", c.eta, c.eps, c.leaf, r.build, r.apply, r.err, r.rank, r.mb, 100 * r.dense_frac, 10);
        std::printf("      Kernauswertungen %.1f %% der direkten, ACA %.1f s, Nachkompression %.1f s (Summe ueber Threads); Vorbereitung %.2f s, dicht %.2f s%c", 100 * r.kfrac, r.taca, r.trec, r.tsetup, r.tdense, 10);
        std::fflush(stdout);
    }
}
