// Unstetig lineare Dichten auf ebenen Dreiecken (v0.45, Stufe 2a der gekruemmten Elemente):
// (1) analytische Innenintegrale fuer lambda_b gegen Brute-Force-Quadratur und gegen die Summenidentitaet;
// (2) Eintraege: sum_ab K_lambda = KernelEntries in allen Zweigen (fern, nah, Ecke, Kante, Selbstterm), psi-Basis;
// (3) Projektion: Mittelwert der linearen Projektion = konstante Projektion; Plemelj E b = b genauer als konstant;
// (4) Streuung an der Kugel gegen die unabhaengige Vorhersage aus Stufe 1 (prototype/curved, Galerkin-Projektion auf feinen
//     ebenen Unterteilungen, extrapoliert): eben/linear Glas -6,163 % (320) und -4,000 % (500), Gold -5,79 % / -3,87 %;
// (5) chirale Kugel: chi = 0 wie achiral, Spiegelsymmetrie sigma_s(chi) = sigma_-s(-chi), nahe der konstanten Rechnung;
// (6) zwei weit getrennte Kugeln: additiv.
#include <cstdio>

#include "check.hpp"
#include "cbem/geometry/quadrature.hpp"
#include "cbem/kernel/triangle_integrals.hpp"
#include "cbem/problems/linear_problem.hpp"
#include "cbem/sources/fields.hpp"

using namespace cbem;

static void test_integrals() {
    const std::array<Vec3, 3> p = {Vec3(0.1, -0.2, 0.05), Vec3(1.3, 0.1, -0.1), Vec3(0.4, 1.1, 0.2)};
    Vec3 n = cross(p[1] - p[0], p[2] - p[0]); n = n / norm(n);
    const Vec3 c = (p[0] + p[1] + p[2]) / 3.0;
    const Vec3 xs[] = {c + n * 0.3, c + n * 0.05, p[2] + n * 0.1, p[1] + (p[1] - c) * 0.4, c - n * 0.2 + Vec3(2, 1, 0)};
    const QuadRule R = QuadRule::subdivided(64);
    const double A = 0.5 * norm(cross(p[1] - p[0], p[2] - p[0]));
    double worst = 0, worst_sum = 0;
    for (const Vec3& x : xs) {
        std::array<Vec3, 3> Ig; std::array<real, 3> Ii; triangle_integrals_linear(x, p, n, Ig, Ii);
        Vec3 Ig0; real Ii0; triangle_integrals(x, p, n, Ig0, Ii0);
        worst_sum = std::max(worst_sum, norm(Ig[0] + Ig[1] + Ig[2] - Ig0) / norm(Ig0) + std::abs(Ii[0] + Ii[1] + Ii[2] - Ii0) / Ii0);
        for (int b = 0; b < 3; ++b) {
            Vec3 G{0, 0, 0}; double I = 0;
            for (std::size_t q = 0; q < R.w.size(); ++q) {
                const Vec3 y = p[0] * R.bary[q][0] + p[1] * R.bary[q][1] + p[2] * R.bary[q][2];
                const Vec3 z = x - y; const double r = norm(z), ww = R.w[q] * A * R.bary[q][b];
                G += z * (ww / (r * r * r)); I += ww / r;
            }
            worst = std::max(worst, norm(G - Ig[b]) / norm(Ig0) + std::abs(I - Ii[b]) / Ii0);
        }
    }
    std::printf("(1) Innenintegrale linear: gegen Brute Force %.1e, Summe gegen konstant %.1e\n", worst, worst_sum);
    CHECK(worst < 1e-8, "lineare Innenintegrale weichen von der Brute-Force-Quadratur ab: %.2e", worst);
    CHECK(worst_sum < 1e-13, "Summe der linearen Innenintegrale ungleich konstant: %.2e", worst_sum);
}

static void test_entries() {
    for (int pass = 0; pass < 2; ++pass) {
        const TriangleMesh m = pass == 0 ? make_icosphere(4) : make_cube_graded(2);
        EntryParams ep; ep.cache_near = false;
        KernelEntries E(m, cplx(1.3, 0.05), ep);
        LinearKernelEntries L(m, cplx(1.3, 0.05), ep);
        double worst[5] = {0, 0, 0, 0, 0}, worst_psi = 0;
        for (std::size_t i = 0; i < m.size(); i += 7)
            for (std::size_t j = 0; j < m.size(); ++j) {
                const int typ = !L.is_near(i, j) ? 0 : 1 + static_cast<int>(L.adjacency(i, j));
                if (typ == 0 && (j % 13)) continue;
                const KernelComp K = E.exact(i, j);
                const LinearBlock B = L.lambda_exact(i, j), P = L.block(i, j);
                double nk = 0, d = 0, dp = 0;
                for (int cc = 0; cc < 4; ++cc) {
                    cplx s = 0, sp = 0; for (int q = 0; q < 9; ++q) { s += B[q][cc]; sp += P[q][cc]; }
                    nk += std::norm(K[cc]); d += std::norm(s - K[cc]);
                    dp += std::norm(sp / 3.0 - K[cc] / std::sqrt(m.area[i] * m.area[j]));
                }
                worst[typ] = std::max(worst[typ], std::sqrt(d / nk));
                worst_psi = std::max(worst_psi, std::sqrt(dp / nk) * std::sqrt(m.area[i] * m.area[j]));
            }
        std::printf("(2) %s: sum K_lambda gegen KernelEntries: fern %.1e, nah %.1e, Ecke %.1e, Kante %.1e, Selbstterm %.1e; psi %.1e\n",
                    pass == 0 ? "Kugel" : "gradierter Wuerfel", worst[0], worst[1], worst[2], worst[3], worst[4], worst_psi);
        for (int t = 0; t < 4; ++t) CHECK(worst[t] < 1e-13, "Summenidentitaet verletzt (Typ %d): %.2e", t, worst[t]);
        // Selbstterm: Kugel Sauter-Schwab (exakt antisymmetrisiert), gestreckte Wuerfelelemente halbanalytisch (Quadraturfehler)
        CHECK(worst[4] < (pass == 0 ? 1e-13 : 3e-5), "Selbstterm: %.2e", worst[4]);
    }
}

static void test_projection_plemelj() {
    const TriangleMesh m = make_icosphere(4);
    const cplx k = 0.5 * 1.33;
    const auto bl = project_plane_wave_linear(m, k, 1.7689, Vec3(0, 0, 1), CVec3{1, 0, 0});
    const auto bc = project_plane_wave(m, k, 1.7689, Vec3(0, 0, 1), CVec3{1, 0, 0});
    const auto cl = linear_to_constant(m, bl);
    double nc = 0, dd = 0; for (std::size_t i = 0; i < bc.size(); ++i) { nc += std::norm(bc[i]); dd += std::norm(cl[i] - bc[i]); }
    LinearKernelEntries Le(m, k); KernelHMatrix Lh(Le, HMatrixParams{}); LinearCauchyOperator LE(m, Lh);
    KernelEntries Ce(m, k); KernelHMatrix Ch(Ce, HMatrixParams{}); CauchyOperator CE(m, Ch);
    std::vector<cplx> yl, yc; LE.apply(bl, yl); CE.apply(bc, yc);
    double nl = 0, dl = 0, dc = 0;
    for (std::size_t i = 0; i < bl.size(); ++i) { nl += std::norm(bl[i]); dl += std::norm(yl[i] - bl[i]); }
    for (std::size_t i = 0; i < bc.size(); ++i) dc += std::norm(yc[i] - bc[i]);
    const double el = std::sqrt(dl / nl), ec = std::sqrt(dc / nc);
    std::printf("(3) Projektion: Mittelwert linear = konstant bis %.1e; Plemelj |E b - b|/|b| linear %.2e, konstant %.2e\n",
                std::sqrt(dd / nc), el, ec);
    CHECK(std::sqrt(dd / nc) < 1e-14, "Mittelwert der linearen Projektion ungleich konstanter Projektion");
    CHECK(el < 1e-3 && el < 0.2 * ec, "Plemelj mit linearen Dichten nicht genauer: %.2e gegen %.2e", el, ec);
}

static void test_sphere_against_stage1() {
    struct Case { const char* name; cplx eps; real om, mie; int n; real pred, tol; };
    // Vorhersage aus Stufe 1 (results/curved_sphere_*.csv, eben/linear, extrapoliert); Toleranz: Extrapolationsunsicherheit
    const Case cs[] = {{"Glas", 2.25, 1.0, 0.2150978, 4, -6.163e-2, 1e-4}, {"Glas", 2.25, 1.0, 0.2150978, 5, -4.000e-2, 1e-4},
                       {"Gold", cplx(-11, 1.2), 0.5, 0.5900182, 4, -5.79e-2, 4e-4}, {"Gold", cplx(-11, 1.2), 0.5, 0.5900182, 5, -3.87e-2, 3e-4}};
    SolveOptions so; so.tol = 1e-8;
    for (const auto& c : cs) {
        const TriangleMesh m = make_icosphere(c.n);
        LinearScatteringProblem L({m}, {Medium{c.eps}}, c.om);
        const PlaneWaveResult r = L.solve_plane_wave(Vec3(0, 0, 1), CVec3{1, 0, 0}, so);
        const real e = r.sigma_ext / pi / c.mie - 1;
        std::printf("(4) %s, %zu Elemente: linear %+.4f %%, Vorhersage Stufe 1 %+.3f %% (%d It.)\n", c.name, m.size(), 100 * e, 100 * c.pred, r.iterations);
        CHECK(std::abs(e - c.pred) < c.tol, "%s %zu: %.5f gegen Vorhersage %.5f", c.name, m.size(), e, c.pred);
        CHECK(std::abs(r.sigma_ext - 4 * pi * std::real(r.forward) / (c.om * c.om)) < 1e-9 * r.sigma_ext, "optisches Theorem");
    }
}

static void test_chiral_and_multibody() {
    const TriangleMesh m = make_icosphere(4);
    SolveOptions so; so.tol = 1e-9;
    const CVec3 pp = circular_polarization(Vec3(0, 0, 1), +1), pm = circular_polarization(Vec3(0, 0, 1), -1);
    LinearScatteringProblem A({m}, {Medium{2.25}}, 1.0);
    LinearScatteringProblem C0({m}, {Medium{2.25, 1.0, 1e-30}}, 1.0);
    const real s_ach = A.solve_plane_wave(Vec3(0, 0, 1), pp, so).sigma_ext, s_c0 = C0.solve_plane_wave(Vec3(0, 0, 1), pp, so).sigma_ext;
    // Spiegelsymmetrie bis auf die ACA-Kompression (sie erhaelt die Symmetrie nicht): eps = 1e-4 -> 4e-6, 1e-6 -> 1e-8
    HMatrixParams hp6; hp6.eps = 1e-6;
    LinearScatteringProblem Cp({m}, {Medium{2.25, 1.0, 0.1}}, 1.0, Medium{}, hp6), Cm({m}, {Medium{2.25, 1.0, -0.1}}, 1.0, Medium{}, hp6);
    const real spp = Cp.solve_plane_wave(Vec3(0, 0, 1), pp, so).sigma_ext, spm = Cp.solve_plane_wave(Vec3(0, 0, 1), pm, so).sigma_ext;
    const real smm = Cm.solve_plane_wave(Vec3(0, 0, 1), pm, so).sigma_ext;
    ScatteringProblem Kc({m}, {Medium{2.25, 1.0, 0.1}}, 1.0);
    const real kpp = Kc.solve_plane_wave(Vec3(0, 0, 1), pp, so).sigma_ext;
    std::printf("(5) chiral: chi=0 %.8f / %.8f; chi=0,1: s+ %.6f, s- %.6f, s-(chi=-0,1) %.6f; konstant s+ %.6f\n", s_ach, s_c0, spp, spm, smm, kpp);
    CHECK(std::abs(s_ach - s_c0) < 1e-8 * s_ach, "chiraler Operator mit chi -> 0 ungleich achiral");
    CHECK(std::abs(spp - smm) < 1e-7 * spp, "Spiegelsymmetrie verletzt: %.2e", std::abs(spp - smm) / spp);
    CHECK(std::abs(spp - spm) > 1e-3 * spp, "kein Zirkulardichroismus");
    CHECK(std::abs(spp / kpp - 1) < 0.08, "linear und konstant weichen zu stark ab");
    // zwei weit getrennte Kugeln: additiv
    const TriangleMesh s = make_icosphere(3);
    LinearScatteringProblem one({s}, {Medium{cplx(-11, 1.2)}}, 0.5), two({s, translated(s, Vec3(30, 0, 0))}, {Medium{cplx(-11, 1.2)}, Medium{cplx(-11, 1.2)}}, 0.5);
    const real s1 = one.solve_plane_wave(Vec3(0, 0, 1), CVec3{1, 0, 0}, so).sigma_ext, s2 = two.solve_plane_wave(Vec3(0, 0, 1), CVec3{1, 0, 0}, so).sigma_ext;
    std::printf("(6) zwei Kugeln im Abstand 30: sigma / (2 sigma_1) = %.5f\n", s2 / (2 * s1));
    CHECK(std::abs(s2 / (2 * s1) - 1) < 0.05, "weit getrennte Koerper nicht additiv");
}

int main() {
    test_integrals();
    test_entries();
    test_projection_plemelj();
    test_sphere_against_stage1();
    test_chiral_and_multibody();
    REPORT();
}
