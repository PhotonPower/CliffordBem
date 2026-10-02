// Gekruemmte (quadratische) Elemente mit unstetig linearen Dichten (v0.47, Stufe 2b):
// (1) Geometrie: Mitten auf den Sehnen = ebenes Netz; Volumen der quadratischen Kugel O(h^4) (wie Stufe 1b);
// (2) Eintraege: ebene Gegenprobe gegen LinearKernelEntries (fern, Ecke, Kante, Selbstterm auf Maschinengenauigkeit, nahe
//     getrennte Paare auf die Genauigkeit der Nahquadraturen);
// (3) Plemelj auf der gekruemmten Flaeche: E h = h fuer innere, E h = -h fuer aeussere Loesungen mit O(h^2) (ein Fehler im
//     Selbstterm, etwa ohne den schwach singulaeren Anteil K_s, gibt O(h) und einen 40-fach groesseren Fehler);
// (4) Streuproblem: ebene Gegenprobe gegen LinearScatteringProblem; Kugel aus Glas und Gold gegen Mie (Fehler bei 320
//     Elementen unter 0,02 %, heute 5,7 % bzw. 7,3 %); chirale Spiegelsymmetrie;
// (5) Nahquadratur mit Singularitaetssubtraktion (v0.49; Gauss-Blattregeln v0.52) gegen das doppelt adaptive Verfahren;
// (6) Fernbloecke direkt in der psi-Basis (v0.50) gegen die Umrechnung der lambda-Bloecke;
// (7) anisotrope Sauter-Schwab-Ordnungen (v0.53) gegen eine hohe isotrope Ordnung.
#include <cstdio>

#include "check.hpp"
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/kernel/dirac_kernel.hpp"
#include "cbem/problems/curved_problem.hpp"
#include "cbem/problems/linear_problem.hpp"
#include "cbem/sources/dipole.hpp"
#include "cbem/sources/fields.hpp"

using namespace cbem;

static void test_geometry() {
    const TriangleMesh f = make_icosphere(4);
    const QuadraticMesh qf = make_quadratic(f, [](const Vec3& p) { return p; });
    double da = 0; for (std::size_t t = 0; t < f.size(); ++t) da = std::max(da, std::abs(qf.area[t] / f.area[t] - 1));
    const double v4 = signed_volume(quadratic_icosphere(4)) / (4 * pi / 3) - 1, v8 = signed_volume(quadratic_icosphere(8)) / (4 * pi / 3) - 1;
    std::printf("(1) eben: Flaechen %.1e; quadratische Kugel: Volumenfehler %.3e (320), %.3e (1280), Faktor %.1f\n", da, v4, v8, v4 / v8);
    CHECK(da < 1e-13, "ebene quadratische Elemente: Flaechen weichen ab");
    CHECK(std::abs(v4 + 1.458e-4) < 2e-6, "Volumenfehler bei 320 Elementen %.4e statt -1.458e-4 (Stufe 1b)", v4);
    CHECK(v4 / v8 > 14 && v4 / v8 < 17, "Volumen konvergiert nicht mit O(h^4): Faktor %.2f", v4 / v8);
}

static CurvedComp to7(const KernelComp& K, const Vec3& n) {
    return CurvedComp{K[1] * n.x + K[2] * n.y + K[3] * n.z, K[0] * n.x, K[0] * n.y, K[0] * n.z,
                      K[1] * n.y - K[2] * n.x, K[1] * n.z - K[3] * n.x, K[2] * n.z - K[3] * n.y};
}

static void test_entries_flat() {
    const TriangleMesh f = make_icosphere(3);
    const QuadraticMesh qf = make_quadratic(f, [](const Vec3& p) { return p; });
    EntryParams ep; ep.cache_near = false;
    // Referenz mit verschaerfter aeusserer Regel: mit der Voreinstellung (0,5) hat die analytische Nahquadratur auf diesem groben
    // Netz selbst 1,7e-5 Fehler; verschaerft bleibt der Boden ihres Restkerns (7 Gauss-Punkte) von etwa 6e-6. Die gekruemmte
    // Quadratur ist dagegen konvergiert (0,3/0,15 und 0,15/0,08 geben dasselbe).
    EntryParams epl = ep; epl.adapt_ratio = 0.1;
    LinearKernelEntries L(f, cplx(1.3, 0.05), epl);
    CurvedNearParams iso; iso.ss_orders = {};   // dieselbe isotrope Sauter-Schwab-Regel wie LinearKernelEntries
    CurvedKernelEntries C(qf, cplx(1.3, 0.05), ep, iso);
    double worst[5] = {0, 0, 0, 0, 0};
    for (std::size_t i = 0; i < f.size(); i += 3)
        for (std::size_t j = 0; j < f.size(); ++j) {
            const int typ = !L.is_near(i, j) ? 0 : 1 + static_cast<int>(L.adjacency(i, j));
            if (typ == 0 && (j % 11)) continue;
            const LinearBlock B = L.lambda_exact(i, j); const CurvedBlock K = C.lambda_exact(i, j);
            double nb = 0, d = 0;
            for (int q = 0; q < 9; ++q) {
                const CurvedComp r = to7(B[q], f.normal[j]);
                for (int c = 0; c < kCurvedComps; ++c) { nb += std::norm(r[c]); d += std::norm(r[c] - K[q][c]); }
            }
            worst[typ] = std::max(worst[typ], std::sqrt(d / nb));
        }
    std::printf("(2) ebene Gegenprobe der Eintraege: fern %.1e, nah getrennt %.1e, Ecke %.1e, Kante %.1e, Selbstterm %.1e\n",
                worst[0], worst[1], worst[2], worst[3], worst[4]);
    for (int t : {0, 2, 3, 4}) CHECK(worst[t] < 1e-13, "Eintraege (Typ %d) weichen ab: %.2e", t, worst[t]);
    CHECK(worst[1] < 1e-5, "nahe getrennte Paare: %.2e", worst[1]);   // beide Nahquadraturen ~2-4e-6
}

template <class F> static std::vector<cplx> trace_curved(const QuadraticMesh& m, F field, cplx eps) {
    CurvedQuadrature Q(m, QuadRule::subdivided(3));
    const auto S = curved_psi_matrices(m);
    std::vector<cplx> h(24 * m.size(), cplx(0));
    for (std::size_t t = 0; t < m.size(); ++t)
        for (int q = 0; q < Q.q; ++q) {
            CVec3 E, H; field(Q.points(t)[q], E, H);
            const Multivector Fm = Multivector::vector(E) * std::sqrt(eps) + Multivector::blade(7) * Multivector::vector(H);
            for (int a = 0; a < 3; ++a) {
                real psi = 0; for (int k = 0; k < 3; ++k) psi += S[t][a * 3 + k] * Q.lam[q][k];
                for (int b = 0; b < 8; ++b) h[8 * (3 * t + a) + b] += Q.weights(t)[q] * psi * Fm.c[b];
            }
        }
    return h;
}

static double rel(const std::vector<cplx>& y, const std::vector<cplx>& h, double s) {
    double n = 0, d = 0; for (std::size_t i = 0; i < h.size(); ++i) { n += std::norm(h[i]); d += std::norm(y[i] - s * h[i]); }
    return std::sqrt(d / n);
}

static void test_plemelj() {
    const real om = 0.5; const Medium W{1.7689}; const cplx k = W.k(om);
    auto plane = [&](const Vec3& x, CVec3& E, CVec3& H) { const cplx ph = std::exp(cplx(0, 1) * k * x.z); E = {ph, 0, 0}; H = {0, std::sqrt(W.eps) * ph, 0}; };
    auto dip = [&](const Vec3& x, CVec3& E, CVec3& H) { dipole_field(W, om, Vec3(0.1, -0.05, 0.08), CVec3{1, 0.5, cplx(0, 0.3)}, x, E, H); };
    HMatrixParams hp; hp.eps = 1e-8;
    double ei[2], ee[2]; int idx = 0;
    for (int n : {4, 6}) {
        const QuadraticMesh q = quadratic_icosphere(n);
        CurvedKernelEntries Ce(q, k); CurvedHMatrix Ch(Ce, hp); CurvedCauchyOperator CE(q, Ch);
        std::vector<cplx> y;
        const auto hi = trace_curved(q, plane, W.eps); CE.apply(hi, y); ei[idx] = rel(y, hi, +1);
        const auto he = trace_curved(q, dip, W.eps); CE.apply(he, y); ee[idx] = rel(y, he, -1);
        ++idx;
    }
    // 320 -> 720 Elemente: O(h^2) = Faktor 2,25, O(h) waere 1,5
    std::printf("(3) Plemelj gekruemmt: innen %.2e -> %.2e (Faktor %.2f), aussen %.2e -> %.2e (Faktor %.2f)\n", ei[0], ei[1], ei[0] / ei[1],
                ee[0], ee[1], ee[0] / ee[1]);
    CHECK(ei[0] < 2e-3 && ee[0] < 1.2e-2, "Plemelj-Fehler zu gross: %.2e / %.2e", ei[0], ee[0]);
    CHECK(ei[0] / ei[1] > 2.0 && ee[0] / ee[1] > 2.0, "Plemelj konvergiert nicht mit O(h^2): %.2f / %.2f", ei[0] / ei[1], ee[0] / ee[1]);
}

static void test_scattering() {
    SolveOptions so; so.tol = 1e-9;
    // ebene Gegenprobe gegen LinearScatteringProblem
    const TriangleMesh f = make_icosphere(3);
    LinearScatteringProblem L({f}, {Medium{2.25}}, 1.0);
    CurvedScatteringProblem Cf({make_quadratic(f, [](const Vec3& p) { return p; })}, {Medium{2.25}}, 1.0);
    const real sl = L.solve_plane_wave(Vec3(0, 0, 1), CVec3{1, 0, 0}, so).sigma_ext, sc = Cf.solve_plane_wave(Vec3(0, 0, 1), CVec3{1, 0, 0}, so).sigma_ext;
    std::printf("(4) ebene Gegenprobe des Streuproblems (180 Elemente): relativ %.1e\n", std::abs(sc / sl - 1));
    CHECK(std::abs(sc / sl - 1) < 2e-5, "gekruemmter Pfad auf ebenen Elementen weicht ab: %.2e", std::abs(sc / sl - 1));
    // Kugel gegen Mie: heute -5,68 % (Glas) bzw. +7,31 % (Gold) bei 320 Elementen
    struct Case { const char* name; cplx eps; real om, mie; };
    for (const Case& c : {Case{"Glas", 2.25, 1.0, 0.2150978}, Case{"Gold", cplx(-11, 1.2), 0.5, 0.5900182}}) {
        CurvedScatteringProblem P({quadratic_icosphere(4)}, {Medium{c.eps}}, c.om);
        const PlaneWaveResult r = P.solve_plane_wave(Vec3(0, 0, 1), CVec3{1, 0, 0}, so);
        const real e = r.sigma_ext / pi / c.mie - 1;
        std::printf("    %s, 320 gekruemmte Elemente: %+.4f %% gegen Mie (%d It.)\n", c.name, 100 * e, r.iterations);
        CHECK(std::abs(e) < 2e-4, "%s: Fehler %.5f gegen Mie zu gross", c.name, e);
        CHECK(std::abs(r.sigma_ext - 4 * pi * std::real(r.forward) / (c.om * c.om)) < 1e-9 * r.sigma_ext, "optisches Theorem");
    }
    // chirale Kugel: Spiegelsymmetrie (ACA eps = 1e-6)
    HMatrixParams hp6; hp6.eps = 1e-6;
    const QuadraticMesh q3 = quadratic_icosphere(3);
    CurvedScatteringProblem Cp({q3}, {Medium{2.25, 1.0, 0.1}}, 1.0, Medium{}, hp6), Cm({q3}, {Medium{2.25, 1.0, -0.1}}, 1.0, Medium{}, hp6);
    const real spp = Cp.solve_plane_wave(Vec3(0, 0, 1), circular_polarization(Vec3(0, 0, 1), +1), so).sigma_ext;
    const real spm = Cp.solve_plane_wave(Vec3(0, 0, 1), circular_polarization(Vec3(0, 0, 1), -1), so).sigma_ext;
    const real smm = Cm.solve_plane_wave(Vec3(0, 0, 1), circular_polarization(Vec3(0, 0, 1), -1), so).sigma_ext;
    std::printf("    chiral: s+ %.6f, s- %.6f, s-(chi=-0,1) %.6f\n", spp, spm, smm);
    CHECK(std::abs(spp - smm) < 1e-6 * spp, "Spiegelsymmetrie verletzt: %.2e", std::abs(spp - smm) / spp);
    CHECK(std::abs(spp - spm) > 1e-3 * spp, "kein Zirkulardichroismus");
}

static void test_subtraction() {
    const QuadraticMesh q = quadratic_icosphere(3);
    EntryParams ep; ep.cache_near = false;
    CurvedNearParams ref; ref.subtract = false; ref.outer_ratio = 0.15; ref.inner_ratio = 0.08; ref.inner_depth = 20;
    CurvedNearParams fast; fast.subtract_outer_ratio = 1.0; fast.correction_ratio = 1.0; fast.outer_rule = 4; fast.correction_rule = 4;
    CurvedNearParams d7; d7.subtract_outer_ratio = 0.3; d7.correction_ratio = 0.3; d7.outer_rule = 0; d7.correction_rule = 0;   // v0.49
    const cplx k(1.3, 0.05);
    CurvedKernelEntries R(q, k, ep, ref), D(q, k, ep), F(q, k, ep, fast), O(q, k, ep, d7);
    double wd = 0, wf = 0, wo = 0;
    for (std::size_t i = 0; i < q.size(); i += 7)
        for (std::size_t j = 0; j < q.size(); ++j) {
            if (!R.is_near(i, j) || R.adjacency(i, j) != Adjacency::None) continue;
            const CurvedBlock A = R.lambda_near(i, j), B = D.lambda_near(i, j), C = F.lambda_near(i, j), G = O.lambda_near(i, j);
            double na = 0, db = 0, dc = 0, dg = 0;
            for (int p = 0; p < 9; ++p)
                for (int c = 0; c < kCurvedComps; ++c) {
                    na += std::norm(A[p][c]); db += std::norm(B[p][c] - A[p][c]); dc += std::norm(C[p][c] - A[p][c]); dg += std::norm(G[p][c] - A[p][c]);
                }
            wd = std::max(wd, std::sqrt(db / na)); wf = std::max(wf, std::sqrt(dc / na)); wo = std::max(wo, std::sqrt(dg / na));
        }
    std::printf("(5) Nahquadratur mit Subtraktion gegen doppelt adaptiv (streng): Voreinstellung (Gauss 5, 1,5) %.1e, schnell (Gauss 4, 1,0) "
                "%.1e, Dunavant 7 mit 0,3 (v0.49) %.1e\n", wd, wf, wo);
    CHECK(wd < 3e-6, "Subtraktion (Voreinstellung) weicht ab: %.2e", wd);
    CHECK(wf < 6e-5, "Subtraktion (schnell) weicht ab: %.2e", wf);
    CHECK(wo < 1e-5, "Subtraktion (Dunavant 7, v0.49) weicht ab: %.2e", wo);
}

static void test_far_blocks() {
    // block_far (Gewichte w psi, zweistufig, dirac_kernel_fast; v0.50) gegen die Umrechnung S L S^T der lambda-Bloecke
    const QuadraticMesh q = quadratic_icosphere(4);
    EntryParams ep; ep.cache_near = false;
    CurvedKernelEntries C(q, cplx(0.09, 1.66), ep);
    double w = 0;
    for (std::size_t i = 0; i < q.size(); i += 5)
        for (std::size_t j = 0; j < q.size(); j += 3) {
            if (C.is_near(i, j)) continue;
            const CurvedBlock L = C.lambda_far(i, j), F = C.block_far(i, j);
            const auto& Si = C.S(i); const auto& Sj = C.S(j);
            double nr = 0, d = 0;
            for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) for (int c = 0; c < kCurvedComps; ++c) {
                cplx r = 0;
                for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) r += Si[k * 3 + a] * L[a * 3 + b][c] * Sj[l * 3 + b];
                nr += std::norm(r); d += std::norm(r - F[k * 3 + l][c]);
            }
            w = std::max(w, std::sqrt(d / nr));
        }
    double wk = 0;
    for (int t = 0; t < 2000; ++t) {
        const Vec3 z(0.001 + 0.004 * t, 0.3 * std::sin(0.1 * t), -0.2);
        for (cplx k : {cplx(0.5), cplx(0.09, 1.66), cplx(2.0, 0.3)}) {
            const KernelValue a = dirac_kernel_full(z, k), b = dirac_kernel_fast(z, k);
            wk = std::max(wk, std::max(std::abs(a.s - b.s) / std::abs(a.s), std::abs(a.vcoef - b.vcoef) / std::abs(a.vcoef)));
        }
    }
    std::printf("(6) Fernbloecke in der psi-Basis gegen S L S^T: %.1e; dirac_kernel_fast gegen dirac_kernel_full: %.1e\n", w, wk);
    CHECK(w < 1e-13, "block_far weicht von der Umrechnung der lambda-Bloecke ab: %.2e", w);
    CHECK(wk < 1e-14, "dirac_kernel_fast weicht ab: %.2e", wk);
}

static void test_sauter_schwab_orders() {
    // anisotrope Sauter-Schwab-Ordnungen (Voreinstellung v0.53) und isotrop 5 (bis v0.52) gegen isotrop 12
    const QuadraticMesh q = quadratic_icosphere(3);
    const cplx k(0.09, 1.66);
    EntryParams ep; ep.cache_near = false;
    EntryParams e12 = ep; e12.ss_order = 12;
    CurvedNearParams iso; iso.ss_orders = {};
    CurvedKernelEntries R(q, k, e12, iso), A(q, k, ep), I(q, k, ep, iso);
    double wa[4] = {0, 0, 0, 0}, wi[4] = {0, 0, 0, 0};
    for (std::size_t i = 0; i < q.size(); i += 5)
        for (std::size_t j = 0; j < q.size(); ++j) {
            const Adjacency adj = R.adjacency(i, j);
            if (adj == Adjacency::None) continue;
            const int t = static_cast<int>(adj);
            const CurvedBlock Kr = R.lambda_exact(i, j), Ka = A.lambda_exact(i, j), Ki = I.lambda_exact(i, j);
            double nr = 0, da = 0, di = 0;
            for (int p = 0; p < 9; ++p)
                for (int c = 0; c < kCurvedComps; ++c) { nr += std::norm(Kr[p][c]); da += std::norm(Ka[p][c] - Kr[p][c]); di += std::norm(Ki[p][c] - Kr[p][c]); }
            wa[t] = std::max(wa[t], std::sqrt(da / nr)); wi[t] = std::max(wi[t], std::sqrt(di / nr));
        }
    std::printf("(7) Sauter-Schwab gegen isotrop 12: anisotrop Ecke %.1e, Kante %.1e, Selbstterm %.1e; isotrop 5: %.1e, %.1e, %.1e\n",
                wa[1], wa[2], wa[3], wi[1], wi[2], wi[3]);
    CHECK(wa[1] < 1e-6 && wa[2] < 2e-6 && wa[3] < 2e-6, "anisotrope Sauter-Schwab-Regel ungenau: %.2e %.2e %.2e", wa[1], wa[2], wa[3]);
    CHECK(wi[3] > 10 * wa[3], "Selbstterm: anisotrope Regel nicht genauer als isotrop 5 (%.2e gegen %.2e)", wa[3], wi[3]);
}

int main() {
    test_sauter_schwab_orders();
    test_far_blocks();
    test_subtraction();
    test_geometry();
    test_entries_flat();
    test_plemelj();
    test_scattering();
    REPORT();
}
