// Optische Kraefte ueber den Spannungstensor (v0.33): Strahlungsdruck auf eine Goldkugel in Wasser gegen Mie
// (tools/mie_force.py), Uebereinstimmung der drei Wege (Randspuren, Kugel, Parallelflaeche), verschwindende Querkraefte;
// Dimer: Bindungskraefte entgegengesetzt gleich, Summe = Kraft auf die Kugel um beide, Vorzeichen je Polarisation.
#include <cmath>
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/optical_force.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, water{1.7689, 1.0, 0.0};
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10; const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0}, py{0.0, 1.0, 0.0};
    {
        ScatteringProblem P({make_icosphere(8)}, {gold}, 0.5, water, hp); const auto r = P.solve_plane_wave(d, px, so);
        const real ref = 10.629788036;                                         // Mie: sigma_pr 1/2 eps |E0|^2
        const Vec3 Ft = force_from_traces(P.mesh(), r.h, water, {0, P.mesh().size()})[0];
        const Vec3 Fs = force_on_sphere(P.mesh(), r.h, water, 0.5, d, px, Vec3(0, 0, 0), 1.5);
        const Vec3 Fo = force_on_offset(P.mesh(), r.h, water, 0.5, d, px, P.mesh(), 0.2);
        std::printf("  Goldkugel: F_z Spuren %.4f, Kugel %.4f, Parallelflaeche %.4f (Mie %.4f)\n", Ft.z, Fs.z, Fo.z, ref);
        CHECK(std::abs(Fs.z / ref - 1) < 0.025, "Strahlungsdruck weicht von Mie ab");
        CHECK(std::abs(Ft.z / Fs.z - 1) < 0.005 && std::abs(Fo.z / Fs.z - 1) < 0.005, "Wege zur Kraft stimmen nicht ueberein");
        CHECK(std::abs(Fs.x) + std::abs(Fs.y) < 1e-4 * Fs.z, "Querkraft auf die Kugel");
    }
    {
        const real sx = 1.15; const TriangleMesh A = translated(make_icosphere(6), Vec3(-sx, 0, 0)), B = translated(make_icosphere(6), Vec3(sx, 0, 0));
        ScatteringProblem P({A, B}, {gold, gold}, 0.5, water, hp);
        for (int pol = 0; pol < 2; ++pol) {
            const CVec3 p = pol == 0 ? px : py; const auto r = P.solve_plane_wave(d, p, so);
            const Vec3 F1 = force_on_offset(P.mesh(), r.h, water, 0.5, d, p, A, 0.06), F2 = force_on_offset(P.mesh(), r.h, water, 0.5, d, p, B, 0.06);
            const Vec3 Fs = force_on_sphere(P.mesh(), r.h, water, 0.5, d, p, Vec3(0, 0, 0), 3.8);
            std::printf("  Dimer, Polarisation %s: F1 = (%+.3f, %.3f), F2 = (%+.3f, %.3f), Summe F_z %.4f, Kugel um beide %.4f\n",
                        pol == 0 ? "entlang " : "senkrecht", F1.x, F1.z, F2.x, F2.z, F1.z + F2.z, Fs.z);
            CHECK(std::abs(F1.x + F2.x) < 1e-3 * std::abs(F1.x), "Bindungskraefte nicht entgegengesetzt gleich");
            CHECK(std::abs((F1.z + F2.z) / Fs.z - 1) < 0.005, "Summe der Einzelkraefte weicht ab");
            CHECK(pol == 0 ? F1.x > 0 : F1.x < 0, "Vorzeichen der Bindungskraft falsch");
        }
    }
    REPORT();
}
