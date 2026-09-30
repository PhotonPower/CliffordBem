// Dipolanregung und Zerfallsraten (v0.26): Kugel aus dem Aussenmedium (Raten = 1), Goldkugel in Wasser gegen die
// Reihenloesung (tools/mie_dipole.py, radial und tangential), verdichtetes Netz bei kleinem Abstand.
#include <cmath>
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/dipole.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, water{1.7689, 1.0, 0.0};
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10;
    const CVec3 pr{1.0, 0.0, 0.0}, pt{0.0, 0.0, 1.0};
    auto rates = [&](const TriangleMesh& m, const Medium& core, real r0, const CVec3& p) {
        ScatteringProblem P({m}, {core}, 0.5, water, hp); const Vec3 x0(r0, 0, 0);
        const auto b = project_dipole(P.mesh(), water, 0.5, x0, p); const auto r = P.solve_rhs(b, so);
        return dipole_rates(P.mesh(), r.h, b, water, 0.5, x0, p); };
    // 1. Kugel aus dem Aussenmedium: keine Aenderung der Raten
    {
        const auto R = rates(make_icosphere(4), water, 1.3, pt);
        std::printf("  Kugel aus dem Aussenmedium: gamma_tot %.8f, gamma_rad %.8f\n", R.total, R.radiative);
        CHECK(std::abs(R.total - 1) < 1e-8 && std::abs(R.radiative - 1) < 1e-6, "Raten ohne Streuer nicht 1");
    }
    // 2. Goldkugel, Abstand 1 Radius (r0 = 2), gleichmaessiges Netz n = 8, gegen die Reihenloesung
    {
        const real ref[2][2] = {{3.236070, 2.967922}, {0.608173, 0.549501}};   // {gesamt, strahlend} radial, tangential
        for (int o = 0; o < 2; ++o) {
            const auto R = rates(make_icosphere(8), gold, 2.0, o == 0 ? pr : pt);
            std::printf("  r0 = 2, %s: gamma_tot %.4f (Ref %.4f, %+.1f %%), gamma_rad %.4f (Ref %.4f, %+.1f %%)\n", o == 0 ? "radial    " : "tangential",
                        R.total, ref[o][0], 100 * (R.total / ref[o][0] - 1), R.radiative, ref[o][1], 100 * (R.radiative / ref[o][1] - 1));
            CHECK(std::abs(R.total / ref[o][0] - 1) < 0.03 && std::abs(R.radiative / ref[o][1] - 1) < 0.03, "Raten bei r0 = 2 weichen ab");
        }
    }
    // 3. kleiner Abstand (d = 0,25): verdichtetes Netz trifft die Gesamtrate deutlich besser als das gleichmaessige
    {
        const real ref_t = 23.010256;                                         // radial, gesamt
        const auto U = rates(make_icosphere(8), gold, 1.25, pr), G = rates(make_icosphere_graded(8, Vec3(1, 0, 0), 0.3), gold, 1.25, pr);
        std::printf("  d = 0,25, radial: gleichmaessig %+.1f %%, verdichtet %+.1f %%\n", 100 * (U.total / ref_t - 1), 100 * (G.total / ref_t - 1));
        CHECK(std::abs(G.total / ref_t - 1) < 0.05, "verdichtetes Netz zu ungenau");
        CHECK(std::abs(G.total / ref_t - 1) < 0.3 * std::abs(U.total / ref_t - 1), "Verdichtung verbessert nichts");
    }
    REPORT();
}
