// Zweitor-Formulierung beschichteter Grenzflaechen (S-Matrix, v0.19): Grenzfall duenner neutraler Schicht, Genauigkeit
// gegen Aden-Kerker fuer duenne und fuer dicke Schichten (d > h, Stabilitaet), chirale Schicht (Spiegelsymmetrie, CD).
#include <cmath>
#include "cbem/problems/twoport_layer_problem.hpp"
#include "cbem/sources/fields.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    const real om = 0.5; const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10;
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, glass{2.25, 1.0, 0.0};
    const int n = 6;
    // 1. duenne neutrale Schicht (Vakuum, d = 1e-3): wie ohne Schicht bis auf die groessere Aussenflaeche
    {
        const auto r = TwoPortLayerProblem(make_icosphere(n), make_icosphere(n, 1.001), gold, Medium{}, 0.001, om, {}, hp).solve_plane_wave(d, px, so);
        const auto b = ScatteringProblem({make_icosphere(n)}, {gold}, om, {}, hp).solve_plane_wave(d, px, so);
        std::printf("  neutrale Schicht d = 0,001: sigma %.6f, ohne Schicht %.6f (rel. %.1e)\n", r.sigma_ext, b.sigma_ext, std::abs(r.sigma_ext / b.sigma_ext - 1));
        CHECK(std::abs(r.sigma_ext / b.sigma_ext - 1) < 1e-3, "neutrale Schicht weicht ab");
    }
    // 2. Glasschale d = 0,05 (d/h = 0,29) und d = 0,2 (d/h = 1,14) gegen Aden-Kerker (sigma exakt 2,336622 bzw. 4,342760)
    {
        const real dd[2] = {0.05, 0.2}, ex[2] = {2.336622483, 4.342760316}, tol[2] = {0.025, 0.01};
        for (int i = 0; i < 2; ++i) {
            const auto r = TwoPortLayerProblem(make_icosphere(n), make_icosphere(n, 1 + dd[i]), gold, glass, dd[i], om, {}, hp).solve_plane_wave(d, px, so);
            std::printf("  Glas d = %.2f: sigma %.6f (exakt %.6f, %+.2f %%), %d It.\n", dd[i], r.sigma_ext, ex[i], 100 * (r.sigma_ext / ex[i] - 1), r.iterations);
            CHECK(std::abs(r.sigma_ext / ex[i] - 1) < tol[i], "Zweitor zu ungenau bei d = %.2f", dd[i]);
            CHECK(r.iterations < 80, "Zweitor: zu viele Iterationen bei d = %.2f (%d)", dd[i], r.iterations);
        }
    }
    // 3. chirale Schale: sigma+(chi) = sigma-(-chi), CD gegen tools/mie_chiral_layered.py (d = 0,05: 7,120609e-3)
    {
        auto sig = [&](real chi, int s) {
            return TwoPortLayerProblem(make_icosphere(n), make_icosphere(n, 1.05), gold, Medium{2.25, 1.0, chi}, 0.05, om, {}, hp)
                .solve_plane_wave(d, circular_polarization(d, s), so).sigma_ext; };
        const real sp = sig(0.1, +1), sm = sig(0.1, -1), sp_m = sig(-0.1, -1), cdx = 7.120608982e-3;
        std::printf("  chirale Schale: CD %.4e (Mie %.4e, %+.1f %%), sigma+(chi) - sigma-(-chi) = %.1e\n", sp - sm, cdx, 100 * ((sp - sm) / cdx - 1), sp - sp_m);
        CHECK(std::abs(sp - sp_m) < 1e-6 * sp, "Spiegelsymmetrie verletzt");
        CHECK(std::abs((sp - sm) / cdx - 1) < 0.1, "CD zu ungenau");
    }
    REPORT();
}
