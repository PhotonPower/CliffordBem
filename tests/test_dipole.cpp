// Dipolanregung und Zerfallsraten (v0.26): Kugel aus dem Aussenmedium (Raten = 1), Goldkugel in Wasser gegen die
// Reihenloesung (tools/mie_dipole.py, radial und tangential), verdichtetes Netz bei kleinem Abstand; Fluoreszenzverstaerkung
// (v0.27) gegen Mie (Nahfeld bei der Anregung, Raten bei der Emission).
#include <cmath>
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/dipole.hpp"
#include "cbem/sources/near_field.hpp"
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
    // 4. Fluoreszenzverstaerkung, Emitter bei r0 = 1,5 (d = 0,5), Anregung x-polarisiert entlang z, q0 = 0,1 (Mie: 33,933)
    {
        const TriangleMesh m = make_icosphere_graded(8, Vec3(1, 0, 0), 0.5);
        ScatteringProblem P({m}, {gold}, 0.5, water, hp); const Vec3 x0(1.5, 0, 0), d(0, 0, 1); const CVec3 pe{1.0, 0.0, 0.0};
        const auto re = P.solve_plane_wave(d, pe, so);
        const auto nf = exterior_near_field(P.mesh(), re.h, water, 0.5, d, pe, {x0});
        real exc[3]; for (int a = 0; a < 3; ++a) exc[a] = std::norm(nf[0].E[a]);
        DipoleRates R[3];
        for (int a = 0; a < 3; ++a) {
            CVec3 p{0.0, 0.0, 0.0}; p[a] = 1.0;
            const auto b = project_dipole(P.mesh(), water, 0.5, x0, p); const auto r = P.solve_rhs(b, so);
            R[a] = dipole_rates(P.mesh(), r.h, b, water, 0.5, x0, p);
        }
        const real F = fluorescence_enhancement(exc, R, 0.1), Fx = 33.932694;
        std::printf("  Fluoreszenz (q0 = 0,1): Anregung %.3f (Mie 7.936), F/F0 %.3f (Mie %.3f, %+.1f %%)\n", exc[0], F, Fx, 100 * (F / Fx - 1));
        CHECK(std::abs(F / Fx - 1) < 0.06, "Fluoreszenzverstaerkung weicht von Mie ab");
        CHECK(std::abs(R[1].total / R[2].total - 1) < 0.02, "Kugelsymmetrie der tangentialen Raten verletzt");
    }
    REPORT();
}
