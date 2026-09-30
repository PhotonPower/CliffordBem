// Nahfeld im Aussenraum (v0.24): Fernfeldgrenze, Goldkugel in Wasser gegen die Mie-Loesung (tools/mie_nearfield.py),
// Markierung innerer und zu naher Punkte, chirales Aussenmedium im Grenzfall chi -> 0, Nahfeld am Zweitor;
// H-Matrix-Auswertung gegen die direkte Summation (v0.25).
#include <cmath>
#include "cbem/problems/twoport_layer_problem.hpp"
#include "cbem/sources/near_field.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, water{1.7689, 1.0, 0.0};
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10;
    const Vec3 d(0, 0, 1); const CVec3 p = circular_polarization(d, +1);
    ScatteringProblem P({make_icosphere(8)}, {gold}, 0.5, water, hp);
    const auto r = P.solve_plane_wave(d, p, so);
    // 1. Fernfeldgrenze: R e^{-ikR} F_s(R x) -> F_inf(x)
    {
        const cplx k = water.k(0.5); const auto b = project_plane_wave(P.mesh(), k, water.eps, d, p);
        std::vector<cplx> hs(r.h.size()); for (std::size_t i = 0; i < hs.size(); ++i) hs[i] = r.h[i] - b[i];
        const Vec3 xh(0.6, 0.0, 0.8); const real R = 400;
        const Multivector Ff = far_field(P.mesh(), hs, k, xh), Fn = scattered_field(P.mesh(), hs, k, {xh * R})[0] * (R * std::exp(-cplx(0, 1) * k * R));
        const real e = std::abs(Fn.c[1] - Ff.c[1]) / std::abs(Ff.c[1]);
        std::printf("  Fernfeldgrenze: rel. Abweichung %.1e\n", e);
        CHECK(e < 5e-3, "Nahfeld geht nicht in das Fernfeld ueber");
    }
    // 2. gegen Mie (n = 8): |E|^2 und optische Chiralitaet an Punkten in 0,05 bis 1 Radien Abstand
    {
        const std::vector<Vec3> pts{Vec3(1.05, 0, 0), Vec3(0, 0, 1.2), Vec3(0.72, 0.576, 0.768), Vec3(2.0, 0, 0)};
        const real e2x[4] = {18.427022283, 0.197785875, 7.143908971, 1.687669249}, cx[4] = {-1.923944110, -0.310347638, -1.304088061, -1.015529397};
        const auto f = exterior_near_field(P.mesh(), r.h, water, 0.5, d, p, pts);
        for (int i = 0; i < 4; ++i) {
            std::printf("  (%.2f, %.2f, %.2f): |E|^2 %.4f (Mie %.4f, %+.1f %%), C/C0 %.4f (Mie %.4f)\n", pts[i].x, pts[i].y, pts[i].z,
                        f[i].enhancement, e2x[i], 100 * (f[i].enhancement / e2x[i] - 1), f[i].chirality, cx[i]);
            CHECK(!f[i].inside && !f[i].too_close, "Aussenpunkt falsch markiert");
            CHECK(std::abs(f[i].enhancement / e2x[i] - 1) < 0.08, "|E|^2 weicht von Mie ab");
            CHECK(std::abs(f[i].chirality - cx[i]) < 0.05 * std::abs(cx[i]) + 0.01, "Chiralitaet weicht von Mie ab");
        }
    }
    // 3. Markierungen: Punkt im Inneren, Punkt auf einem Knoten des Netzes
    {
        const auto f = exterior_near_field(P.mesh(), r.h, water, 0.5, d, p, {Vec3(0.3, 0.2, -0.1), P.mesh().P[5]});
        CHECK(f[0].inside, "innerer Punkt nicht markiert");
        CHECK(f[1].too_close, "Punkt auf dem Netz nicht markiert");
    }
    // 4. chirales Aussenmedium, chi -> 0: gleiches Nahfeld
    {
        const Medium w0{1.7689, 1.0, 1e-9};
        ScatteringProblem Q({make_icosphere(4)}, {gold}, 0.5, w0, hp), Q0({make_icosphere(4)}, {gold}, 0.5, water, hp);
        const auto rq = Q.solve_plane_wave(d, p, so), r0 = Q0.solve_plane_wave(d, p, so);
        const std::vector<Vec3> pts{Vec3(1.3, 0.2, 0.1)};
        const auto a = exterior_near_field(Q.mesh(), rq.h, w0, 0.5, d, p, pts), b = exterior_near_field(Q0.mesh(), r0.h, water, 0.5, d, p, pts);
        std::printf("  chi2 = 1e-9: |E|^2 %.8f gegen achiral %.8f\n", a[0].enhancement, b[0].enhancement);
        CHECK(std::abs(a[0].enhancement / b[0].enhancement - 1) < 1e-7, "chirales Aussenmedium: Grenzfall chi -> 0 verfehlt");
    }
    // 5. Zweitor: Glasschale auf Gold, Nahfeld ausserhalb der Schale gegen Mie (geschichtete Kugel)
    {
        const Medium glass{2.25, 1.0, 0.0};
        TwoPortLayerProblem T(make_icosphere(8), make_icosphere(8, 1.05), gold, glass, 0.05, 0.5, water, hp);
        const auto rt = T.solve_plane_wave(d, p, so);
        std::vector<cplx> h(rt.h.begin(), rt.h.begin() + 8 * T.outer_mesh().size());
        const auto f = exterior_near_field(T.outer_mesh(), h, water, 0.5, d, p, {Vec3(1.2, 0, 0)});
        const real ex = 10.700987268;                                     // tools/mie_nearfield.py, radii [1, 1.05]
        std::printf("  Zweitor, Glasschale: |E|^2 bei (1,2, 0, 0) %.4f (Mie %.4f, %+.1f %%)\n", f[0].enhancement, ex, 100 * (f[0].enhancement / ex - 1));
        CHECK(std::abs(f[0].enhancement / ex - 1) < 0.05, "Nahfeld am Zweitor weicht von Mie ab");
    }
    // 6. H-Matrix gegen direkte Summation (Gitter in der xz-Ebene, auch Punkte im Inneren)
    {
        const cplx k = water.k(0.5); const auto b = project_plane_wave(P.mesh(), k, water.eps, d, p);
        std::vector<cplx> hs(r.h.size()); for (std::size_t i = 0; i < hs.size(); ++i) hs[i] = r.h[i] - b[i];
        std::vector<Vec3> pts;
        for (int j = 0; j < 40; ++j) for (int i = 0; i < 80; ++i) pts.push_back(Vec3(-2.5 + 5.0 * i / 79, 0.0, -1.5 + 3.0 * j / 39));
        const auto Fd = scattered_field(P.mesh(), hs, k, pts);
        HMatrixParams q; q.eps = 1e-4; const NearFieldOperator H(P.mesh(), k, pts, q); const auto Fh = H.apply(hs);
        real emax = 0, fmax = 0;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            if (winding_number(P.mesh(), pts[i]) > 0.5) continue;
            real e = 0, f = 0; for (int c = 0; c < 8; ++c) { e += std::norm(Fh[i].c[c] - Fd[i].c[c]); f += std::norm(Fd[i].c[c]); }
            emax = std::max(emax, std::sqrt(e)); fmax = std::max(fmax, std::sqrt(f));
        }
        const auto& st = H.stats(); const real mem = (st.entries_dense + st.entries_lowrank) / (4.0 * pts.size() * P.mesh().size());
        std::printf("  H-Matrix (eps 1e-4): max Fehler / max|F| = %.1e, Speicher %.1f %% der dichten Matrix, Rang %.1f\n", emax / fmax, 100 * mem, st.mean_rank);
        CHECK(emax / fmax < 1e-3, "H-Matrix weicht von der direkten Summation ab");
        CHECK(mem < 0.6, "H-Matrix komprimiert nicht");
    }
    REPORT();
}
