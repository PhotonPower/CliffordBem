// Mehrkoerperprobleme: grosser Abstand -> Summe der Einzelquerschnitte; blockdiagonaler gegen vereinigten
// Innenoperator; Vertauschungssymmetrie verschiedener Medien.
#include "cbem/problems/scattering_problem.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    const real om = 0.5; const TriangleMesh s = make_icosphere(3); const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10;
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, glass{2.25, 1.0, 0.0};
    real s1 = ScatteringProblem({s}, {gold}, om, {}, hp).solve_plane_wave(d, px, so).sigma_ext;
    auto pair = [&](real D, const Medium& a, const Medium& b, bool uni) {
        return ScatteringProblem({translated(s, Vec3(-D / 2, 0, 0)), translated(s, Vec3(D / 2, 0, 0))}, {a, b}, om, {}, hp, {}, uni).solve_plane_wave(d, px, so).sigma_ext; };
    real s40 = pair(40, gold, gold, false);
    std::printf("  D = 40: sigma = %.6f, 2 sigma_1 = %.6f (rel. %.1e)\n", s40, 2 * s1, std::abs(s40 - 2 * s1) / (2 * s1));
    CHECK(std::abs(s40 - 2 * s1) < 1e-3 * 2 * s1, "grosser Abstand: keine Additivitaet");
    real sb = pair(3, gold, gold, false), su = pair(3, gold, gold, true);
    std::printf("  D = 3: blockdiagonal %.6f, vereinigt %.6f (rel. %.1e)\n", sb, su, std::abs(sb - su) / su);
    CHECK(std::abs(sb - su) < 2e-3 * su, "blockdiagonal und vereinigt weichen ab");
    real sgg = pair(3, gold, glass, false), sgl = pair(3, glass, gold, false);
    std::printf("  Gold/Glas %.8f, Glas/Gold %.8f\n", sgg, sgl);
    CHECK(std::abs(sgg - sgl) < 1e-7 * sgg, "Vertauschungssymmetrie verletzt");
    REPORT();
}
