// Nahfeld und Kraefte auf gekruemmten Elementen (v0.58): Projektion beliebiger Felder gegen die der ebenen Welle,
// Fernfeldgrenze, Goldkugel in Wasser gegen Mie (Nahfeld: Werte aus tests/test_near_field.cpp / tools/mie_nearfield.py;
// Strahlungsdruck: tests/test_optical_force.cpp / tools/mie_force.py), Kugel- und Parallelflaeche, Markierungen,
// Konvergenzordnung.
#include <cmath>
#include <cstdio>

#include "check.hpp"
#include "cbem/problems/curved_problem.hpp"
#include "cbem/sources/curved_near_field.hpp"

using namespace cbem;

int main() {
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, water{1.7689, 1.0, 0.0};
    const real om = 0.5;
    SolveOptions so; so.tol = 1e-10;
    const Vec3 d(0, 0, 1); const CVec3 pc = circular_polarization(d, +1), px{1.0, 0.0, 0.0};
    const QuadraticMesh q4 = quadratic_icosphere(4), q6 = quadratic_icosphere(6);
    // 1. Projektion eines allgemeinen Feldes = Projektion der ebenen Welle (gleiche Regel)
    {
        const auto b1 = project_incident_curved(q4, PlaneWaveField(water, om, d, px), water);
        const auto b2 = project_plane_wave_curved(q4, water.k(om), water.eps, d, px);
        double e = 0, n = 0; for (std::size_t i = 0; i < b1.size(); ++i) { e = std::max(e, std::abs(b1[i] - b2[i])); n = std::max(n, std::abs(b2[i])); }
        std::printf("  Projektion allgemein gegen ebene Welle: %.1e\n", e / n);
        CHECK(e < 1e-13 * n, "project_incident_curved weicht ab: %.2e", e / n);
    }
    CurvedScatteringProblem P4({q4}, {gold}, om, water), P6({q6}, {gold}, om, water);
    // 2. Fernfeldgrenze: R e^{-ikR} F_s(R x) -> F_inf(x)
    {
        const auto r = P4.solve_plane_wave(d, pc, so);
        const cplx k = water.k(om); const auto b = project_plane_wave_curved(q4, k, water.eps, d, pc);
        std::vector<cplx> hs(r.h.size()); for (std::size_t i = 0; i < hs.size(); ++i) hs[i] = r.h[i] - b[i];
        const Vec3 xh(0.6, 0.0, 0.8); const real R = 400;
        const Multivector Ff = far_field_curved(q4, hs, k, xh), Fn = scattered_field_curved(q4, hs, k, {xh * R})[0] * (R * std::exp(-cplx(0, 1) * k * R));
        const real e = std::abs(Fn.c[1] - Ff.c[1]) / std::abs(Ff.c[1]);
        std::printf("  Fernfeldgrenze: rel. Abweichung %.1e\n", e);
        CHECK(e < 5e-3, "Nahfeld geht nicht in das Fernfeld ueber: %.2e", e);
    }
    // 3. Nahfeld gegen Mie (zirkular), 320 und 720 Elemente: |E|^2 und optische Chiralitaet
    const std::vector<Vec3> pts{Vec3(1.05, 0, 0), Vec3(0, 0, 1.2), Vec3(0.72, 0.576, 0.768), Vec3(2.0, 0, 0)};
    const real e2x[4] = {18.427022283, 0.197785875, 7.143908971, 1.687669249}, cx[4] = {-1.923944110, -0.310347638, -1.304088061, -1.015529397};
    double we[2] = {0, 0};
    int idx = 0;
    for (CurvedScatteringProblem* P : {&P4, &P6}) {
        const auto r = P->solve_plane_wave(d, pc, so);
        const auto f = exterior_near_field_curved(P->mesh(), r.h, water, om, d, pc, pts);
        for (int i = 0; i < 4; ++i) {
            const double ee = std::abs(f[i].enhancement / e2x[i] - 1), ec = std::abs(f[i].chirality / cx[i] - 1);
            we[idx] = std::max(we[idx], std::max(ee, ec));
            CHECK(!f[i].inside && !f[i].too_close, "Aussenpunkt falsch markiert");
        }
        std::printf("  Nahfeld %zu Elemente: groesste Abweichung von Mie (|E|^2, C) %.3f %%\n", P->mesh().size(), 100 * we[idx]);
        ++idx;
    }
    CHECK(we[0] < 0.01, "Nahfeld (320 Elemente) weicht von Mie ab: %.3f %%", 100 * we[0]);
    CHECK(we[1] < 0.003, "Nahfeld (720 Elemente) weicht von Mie ab: %.3f %%", 100 * we[1]);
    CHECK(we[0] / we[1] > 3.0, "Nahfeld konvergiert nicht schneller als O(h^2): Faktor %.2f", we[0] / we[1]);   // O(h^4): 5,1
    // 4. Strahlungsdruck (linear) gegen Mie: Kugel und Parallelflaeche, keine Querkraft
    double wf[2] = {0, 0};
    idx = 0;
    const real ref = 10.629788036;
    for (CurvedScatteringProblem* P : {&P4, &P6}) {
        const auto r = P->solve_plane_wave(d, px, so);
        const NearFieldEval nf = make_plane_wave_eval_curved(P->mesh(), r.h, water, om, d, px);
        const Vec3 Fs = force_on_sphere(nf, water, Vec3(0, 0, 0), 1.5), Fo = force_on_offset(nf, water, P->mesh().flat, 0.2);
        wf[idx] = std::abs(Fs.z / ref - 1);
        std::printf("  Strahlungsdruck %zu Elemente: Kugel %+.4f %%, Parallelflaeche %+.4f %% gegen Mie\n", P->mesh().size(), 100 * (Fs.z / ref - 1),
                    100 * (Fo.z / ref - 1));
        CHECK(std::abs(Fo.z / Fs.z - 1) < 1e-3, "Kugel und Parallelflaeche stimmen nicht ueberein");
        CHECK(std::abs(Fs.x) + std::abs(Fs.y) < 1e-6 * Fs.z, "Querkraft auf die Kugel");
        ++idx;
    }
    CHECK(wf[0] < 3e-3 && wf[1] < 6e-4, "Strahlungsdruck weicht von Mie ab: %.2e / %.2e", wf[0], wf[1]);
    CHECK(wf[0] / wf[1] > 3.0, "Strahlungsdruck konvergiert nicht schneller als O(h^2): Faktor %.2f", wf[0] / wf[1]);
    // 5. Markierungen: Mittelpunkt innen; Punkt zwischen Sehne und gekruemmter Flaeche zu nah
    {
        const auto r = P4.solve_plane_wave(d, px, so);
        const Vec3 c = q4.flat.centroid[0], n = q4.flat.normal[0];
        const auto f = exterior_near_field_curved(q4, r.h, water, om, d, px, {Vec3(0, 0, 0), c + n * (0.5 * q4.bulge[0])});
        std::printf("  Markierungen: Mitte innen %d, Woelbung zu nah %d\n", int(f[0].inside), int(f[1].too_close));
        CHECK(f[0].inside, "Mittelpunkt nicht als innen markiert");
        CHECK(f[1].too_close, "Punkt in der Woelbung nicht als zu nah markiert");
    }
    // 6. chirales Aussenmedium (v0.59): Grenzfall chi -> 0 des Nahfelds (Streufeld je Helizitaet mit k_pm)
    {
        const Medium h0{1.7689, 1.0, 0.0}, h1{1.7689, 1.0, 1e-9};
        CurvedScatteringProblem A({q4}, {gold}, om, h0), B({q4}, {gold}, om, h1);
        const auto ra = A.solve_plane_wave(d, pc, so), rb = B.solve_plane_wave(d, pc, so);
        const std::vector<Vec3> x{Vec3(1.05, 0, 0), Vec3(0.72, 0.576, 0.768)};
        const auto fa = exterior_near_field_curved(q4, ra.h, h0, om, d, pc, x), fb = exterior_near_field_curved(q4, rb.h, h1, om, d, pc, x);
        double e = 0;
        for (int i = 0; i < 2; ++i) e = std::max(e, std::max(std::abs(fb[i].enhancement / fa[i].enhancement - 1), std::abs(fb[i].chirality / fa[i].chirality - 1)));
        std::printf("  chirales Aussenmedium, chi = 1e-9 gegen achiral: %.1e\n", e);
        CHECK(e < 1e-6, "Nahfeld im Grenzfall chi -> 0 verfehlt: %.2e", e);
    }
    // 7. H-Matrix fuer viele Punkte (v0.61) gegen die direkte Summation: Karte in y = 0 bis dicht an die Kugel, eps 1e-6 und
    // 1e-4 (Fehler bezogen auf max |F_s|), dazu das Gesamtfeld im chiralen Aussenmedium (zwei Helizitaeten)
    {
        std::vector<Vec3> x;
        for (int i = 0; i < 72; ++i)
            for (int j = 0; j < 72; ++j) { const Vec3 y(-2 + 4.0 * (i + 0.5) / 72, 0, -2 + 4.0 * (j + 0.5) / 72); if (norm(y) > 1.01) x.push_back(y); }
        const auto r = P4.solve_plane_wave(d, pc, so);
        const cplx k = water.k(om); const auto b = project_plane_wave_curved(q4, k, water.eps, d, pc);
        std::vector<cplx> hs(r.h.size()); for (std::size_t i = 0; i < hs.size(); ++i) hs[i] = r.h[i] - b[i];
        const auto F0 = scattered_field_curved(q4, hs, k, x);
        auto rel = [&](const std::vector<Multivector>& F) {
            real e = 0, n = 0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                real ei = 0, ni = 0; for (int c = 0; c < 8; ++c) { ei += std::norm(F[i].c[c] - F0[i].c[c]); ni += std::norm(F0[i].c[c]); }
                e = std::max(e, std::sqrt(ei)); n = std::max(n, std::sqrt(ni));
            }
            return e / n;
        };
        HStats st6, st4;
        const real e6 = rel(scattered_field_curved_hmatrix(q4, hs, k, x, 1e-6, &st6)), e4 = rel(scattered_field_curved_hmatrix(q4, hs, k, x, 1e-4, &st4));
        std::printf("  H-Matrix, %zu Punkte: eps 1e-6 Fehler %.1e (%zu Niedrigrangbloecke, Rang %.1f), eps 1e-4 Fehler %.1e (Rang %.1f)\n", x.size(), e6,
                    st6.n_lowrank, st6.mean_rank, e4, st4.mean_rank);
        CHECK(st6.n_lowrank > 0 && st6.n_dense > 0, "H-Matrix ohne Niedrigrang- oder dichte Bloecke");
        CHECK(e6 < 2e-6, "H-Matrix (eps 1e-6) weicht von der direkten Summation ab: %.2e", e6);
        CHECK(e4 < 2e-4 && e4 > e6, "H-Matrix (eps 1e-4): Fehler %.2e", e4);
        // chirales Aussenmedium: Gesamtfeld ueber exterior_near_field_curved, H-Matrix (Voreinstellung ab 4000 Punkten) gegen direkt
        const Medium hc{1.7689, 1.0, 0.05};
        CurvedScatteringProblem Pc({q4}, {gold}, om, hc);
        const auto rc = Pc.solve_plane_wave(d, pc, so);
        NearFieldOptions direct = curved_near_field_options(); direct.hmatrix_min_points = x.size() + 1;
        const auto fh = exterior_near_field_curved(q4, rc.h, hc, om, d, pc, x), fd = exterior_near_field_curved(q4, rc.h, hc, om, d, pc, x, direct);
        real ec = 0, nc = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
            for (int a = 0; a < 3; ++a) { ec = std::max(ec, std::abs(fh[i].E[a] - fd[i].E[a])); nc = std::max(nc, std::abs(fd[i].E[a])); }
        std::printf("  H-Matrix, chirales Aussenmedium: Gesamtfeld E %.1e\n", ec / nc);
        CHECK(x.size() >= curved_near_field_options().hmatrix_min_points, "zu wenige Punkte fuer die H-Matrix");
        CHECK(ec < 2e-6 * nc, "H-Matrix im chiralen Aussenmedium weicht ab: %.2e", ec / nc);
    }
    REPORT();
}
