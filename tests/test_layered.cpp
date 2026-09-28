// Beschichtete Grenzflaechen: Parallelflaechen, Rueckfuehrung auf T_1 ohne Schicht, beschichtete Kugel gegen
// Aden-Kerker (tools/mie_coated.py), Phase der Vorwaertsamplitude, chirale Schicht (Spiegelsymmetrie).
#include <cmath>
#include <stdexcept>
#include "cbem/problems/layered_problem.hpp"
#include "cbem/sources/fields.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    // 1. Parallelflaechen: Wuerfel -> jede Seite genau um d verschoben; Kugel -> Radius 1 + d; zu grosser Versatz wirft
    {
        TriangleMesh c = offset_surface(make_cube_uniform(4), 0.1); real e = 0;
        for (auto& p : c.P) e = std::max(e, std::abs(std::max({std::abs(p.x), std::abs(p.y), std::abs(p.z)}) - 1.1));
        std::printf("  Wuerfel +0,1: max. Abweichung der Seiten %.1e\n", e);
        CHECK(e < 1e-12, "Wuerfel-Parallelflaeche: Seiten nicht um d verschoben (%.1e)", e);
        TriangleMesh ci = offset_surface(make_cube_graded(3), -0.05); e = 0;
        for (auto& p : ci.P) e = std::max(e, std::abs(std::max({std::abs(p.x), std::abs(p.y), std::abs(p.z)}) - 0.95));
        CHECK(e < 1e-12, "gradierter Wuerfel nach innen: Abweichung %.1e", e);
        TriangleMesh s = offset_surface(make_icosphere(6), 0.05); real er = 0;
        for (auto& p : s.P) er = std::max(er, std::abs(norm(p) - 1.05));
        std::printf("  Kugel +0,05: max. |r - 1,05| = %.1e\n", er);
        CHECK(er < 1e-3, "Kugel-Parallelflaeche: Radiusfehler %.1e", er);
        bool thrown = false; try { offset_surface(make_icosphere(3), -1.2); } catch (const std::runtime_error&) { thrown = true; }
        CHECK(thrown, "zu grosser Versatz nach innen nicht erkannt");
    }
    const real om = 0.5; const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10;
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, glass{2.25, 1.0, 0.0};
    // 2. ohne Schicht identisch mit ScatteringProblem (ein und zwei Koerper)
    {
        const TriangleMesh s = make_icosphere(3);
        LayeredGeometry g; add_body(g, s, gold);
        auto a = LayeredScatteringProblem(g, om, hp, EntryParams{}).solve_plane_wave(d, px, so);
        auto b = ScatteringProblem({s}, {gold}, om, {}, hp).solve_plane_wave(d, px, so);
        std::printf("  ohne Schicht: %.10f / %.10f, S(0) %.3e Abweichung\n", a.sigma_ext, b.sigma_ext, std::abs(a.forward - b.forward));
        CHECK(std::abs(a.sigma_ext - b.sigma_ext) < 1e-9 * b.sigma_ext && std::abs(a.forward - b.forward) < 1e-9 * std::abs(b.forward), "T_1 nicht reproduziert");
        const TriangleMesh s1 = translated(s, Vec3(-1.5, 0, 0)), s2 = translated(s, Vec3(1.5, 0, 0));
        LayeredGeometry g2; add_body(g2, s1, gold); add_body(g2, s2, glass);
        auto c = LayeredScatteringProblem(g2, om, hp, EntryParams{}).solve_plane_wave(d, px, so);
        auto e = ScatteringProblem({s1, s2}, {gold, glass}, om, {}, hp).solve_plane_wave(d, px, so);
        CHECK(std::abs(c.sigma_ext - e.sigma_ext) < 1e-9 * e.sigma_ext, "Mehrkoerper: %.10f statt %.10f", c.sigma_ext, e.sigma_ext);
    }
    // 3. Goldkern (Radius 1) mit Glasschale bis 1,2 gegen Aden-Kerker; Ordnung 2 in n
    {
        const real Qm = 0.95996079, Sre = 0.08639647, Sim = -0.30162585;          // tools/mie_coated.py
        real err[2]; int k = 0;
        for (int n : {3, 4}) {
            LayeredGeometry g; add_layered_body(g, {make_icosphere(n, 1.2), make_icosphere(n, 1.0)}, {glass, gold});
            auto r = LayeredScatteringProblem(g, om, hp).solve_plane_wave(d, px, so);
            const real Q = r.sigma_ext / (pi * 1.44);
            std::printf("  Gold@Glas n=%d: Q = %.6f (Mie %.6f), S(0) = %.5f%+.5fi (Mie %.5f%+.5fi), %d It.\n", n, Q, Qm, r.forward.real(), r.forward.imag(), Sre, Sim, r.iterations);
            err[k++] = std::abs(r.forward - cplx(Sre, Sim));
            CHECK(std::abs(Q - Qm) < 0.04 * Qm, "Q_ext zu weit von Aden-Kerker");
            CHECK(std::abs(4 * pi * r.forward.real() / (om * om) - r.sigma_ext) < 1e-8 * r.sigma_ext, "optisches Theorem: Re S(0) passt nicht zu sigma_ext");
        }
        std::printf("  |S - S_Mie|: n=3 %.2e, n=4 %.2e (Verhaeltnis %.2f, Ordnung 2: %.2f)\n", err[0], err[1], err[0] / err[1], 16.0 / 9.0);
        CHECK(err[0] / err[1] > 1.5, "keine Konvergenz der Vorwaertsamplitude");
    }
    // 4. Schale aus Aussenmedium ("neutral") aendert nichts bis auf die Diskretisierung: Abweichung faellt mit dem Netz
    {
        real rel[2]; int k = 0;
        for (int n : {3, 4}) {
            LayeredGeometry g0; add_body(g0, make_icosphere(n), gold);
            const real s0 = LayeredScatteringProblem(g0, om, hp).solve_plane_wave(d, px, so).sigma_ext;
            LayeredGeometry g1; add_layered_body(g1, {make_icosphere(n, 1.1), make_icosphere(n)}, {Medium{}, gold});
            const real s1 = LayeredScatteringProblem(g1, om, hp).solve_plane_wave(d, px, so).sigma_ext;
            rel[k++] = std::abs(s1 - s0) / s0;
            std::printf("  neutrale Schale n=%d: %.6f gegen %.6f (rel. %.1e)\n", n, s1, s0, rel[k - 1]);
        }
        CHECK(rel[1] < 0.02 && rel[1] < 0.8 * rel[0], "neutrale Schale: Abweichung faellt nicht mit dem Netz");
    }
    // 5. chirale Schale: sigma_s(chi) = sigma_{-s}(-chi)
    {
        auto sig = [&](real chi, int s) {
            LayeredGeometry g; add_coated_body(g, make_icosphere(3), gold, {Coating{0.15, Medium{2.25, 1.0, chi}}});
            return LayeredScatteringProblem(g, om, hp).solve_plane_wave(d, circular_polarization(d, s), so).sigma_ext; };
        const real a = sig(0.2, +1), b = sig(-0.2, -1), c = sig(0.2, -1);
        std::printf("  chirale Schale: s+(0,2) = %.8f, s-(-0,2) = %.8f, s-(0,2) = %.8f\n", a, b, c);
        CHECK(std::abs(a - b) < 1e-7 * a, "Spiegelsymmetrie verletzt");
        CHECK(std::abs(a - c) > 1e-4 * a, "chirale Schale ohne Zirkulardichroismus");
    }
    REPORT();
}
