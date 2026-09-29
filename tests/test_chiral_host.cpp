// Chirales Aussenmedium (v0.21): Helizitaetswellen mit k_pm, Aussenoperator P+ E_{k+} + P- E_{k-}, optisches Theorem je
// Kanal. Pruefungen: Kugel aus dem Aussenmedium unsichtbar, lineare Polarisation abgelehnt, Grenzfall chi -> 0,
// Spiegelsymmetrie, Goldkugel gegen tools/mie_chiral_layered.py (chi2), Zweitor gegen exakte Mehrschichtmethode.
#include <cmath>
#include <stdexcept>
#include "cbem/problems/layered_problem.hpp"
#include "cbem/problems/twoport_layer_problem.hpp"
#include "cbem/sources/fields.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    const real om = 0.5; const Vec3 d(0, 0, 1);
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10;
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, host{1.7689, 1.0, 0.05};
    const TriangleMesh s = make_icosphere(6);
    auto sig = [&](const Medium& core, const Medium& out, int h) {
        return ScatteringProblem({s}, {core}, om, out, hp).solve_plane_wave(d, circular_polarization(d, h), so).sigma_ext; };
    // 1. Kugel aus dem Aussenmedium ist unsichtbar
    {
        const real a = sig(host, host, +1);
        std::printf("  Kugel aus dem chiralen Aussenmedium: sigma+ = %.1e\n", a);
        CHECK(std::abs(a) < 1e-10, "Kugel aus dem Aussenmedium streut");
    }
    // 2. lineare Polarisation ist keine Eigenmode
    {
        bool thrown = false;
        try { ScatteringProblem({s}, {gold}, om, host, hp).solve_plane_wave(d, CVec3{1.0, 0.0, 0.0}, so); } catch (const std::invalid_argument&) { thrown = true; }
        CHECK(thrown, "lineare Polarisation im chiralen Aussenmedium angenommen");
    }
    // 3. Grenzfall chi -> 0 und Spiegelsymmetrie sigma_+(chi) = sigma_-(-chi)
    {
        const real a = sig(gold, Medium{1.7689, 1.0, 1e-9}, +1), b = sig(gold, Medium{1.7689, 1.0, 0.0}, +1);
        std::printf("  chi2 = 1e-9 gegen achiral: rel. %.1e\n", std::abs(a / b - 1));
        CHECK(std::abs(a / b - 1) < 1e-7, "Grenzfall chi -> 0 verfehlt");
        const real p = sig(gold, host, +1), m = sig(gold, Medium{1.7689, 1.0, -0.05}, -1);
        std::printf("  sigma+(chi2) - sigma-(-chi2) = %.1e\n", p - m);
        CHECK(std::abs(p - m) < 1e-6 * p, "Spiegelsymmetrie verletzt");
    }
    // 4. Goldkugel in chiralem Wasser gegen Mie (n = 6: sigma -2,4 %, CD -9,7 %, Ordnung 2)
    {
        const real mp = 11.765497482627687, mm = 11.91281031526871;
        const real p = sig(gold, host, +1), m = sig(gold, host, -1);
        std::printf("  Goldkugel: sigma+ %+.2f %%, sigma- %+.2f %%, CD %.4e (Mie %.4e, %+.1f %%)\n", 100 * (p / mp - 1), 100 * (m / mm - 1), p - m, mp - mm, 100 * ((p - m) / (mp - mm) - 1));
        CHECK(std::abs(p / mp - 1) < 0.035 && std::abs(m / mm - 1) < 0.035, "sigma im chiralen Aussenmedium zu ungenau");
        CHECK(std::abs((p - m) / (mp - mm) - 1) < 0.15, "CD im chiralen Aussenmedium zu ungenau");
    }
    // 5. Glasschale in chiralem Wasser: Zweitor gegen exakte Mehrschichtmethode (CD)
    {
        const Medium glass{2.25, 1.0, 0.0};
        TwoPortLayerProblem T(s, make_icosphere(6, 1.05), gold, glass, 0.05, om, host, hp);
        LayeredGeometry g(host); add_layered_body(g, {make_icosphere(6, 1.05), s}, {glass, gold});
        LayeredScatteringProblem L(g, om, hp);
        auto cd = [&](auto& P) { return P.solve_plane_wave(d, circular_polarization(d, +1), so).sigma_ext - P.solve_plane_wave(d, circular_polarization(d, -1), so).sigma_ext; };
        const real a = cd(T), b = cd(L);
        std::printf("  Glasschale in chiralem Wasser: CD Zweitor %.5e, exakt %.5e (rel. %.1e)\n", a, b, std::abs(a / b - 1));
        CHECK(std::abs(a / b - 1) < 0.02, "Zweitor und exakte Methode weichen ab");
    }
    REPORT();
}
