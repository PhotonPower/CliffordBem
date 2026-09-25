// Dimer aus zwei zu den Kanten gradierten Wuerfeln [-1,1]^3 (Abstand der Mittelpunkte D entlang x), Einfall
// entlang z: punktweise gegen Kanten-/Eck-Blockvorkonditionierung fuer mehrere Koerper.
// Beispiel: scatter_cube_dimer --L 3,5 --gap 0.5 --omega 0.5 --eps1 -11,1.2 --R 0.125
#include <chrono>
#include <cstdio>
#include <sstream>
#include <string>
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/fields.hpp"
using namespace cbem;
static std::vector<int> ilist(const std::string& s) { std::vector<int> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, ',')) v.push_back(std::stoi(t)); return v; }
int main(int argc, char** argv) {
    std::vector<int> Ls = {3, 5}; double gap = 0.5, om = 0.5, R = 0.125; cplx eps(-11, 1.2);
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--L") Ls = ilist(nxt()); else if (o == "--gap") gap = std::stod(nxt()); else if (o == "--omega") om = std::stod(nxt());
        else if (o == "--R") R = std::stod(nxt());
        else if (o == "--eps1") { std::string s = nxt(); auto c = s.find(','); eps = cplx(std::stod(s.substr(0, c)), std::stod(s.substr(c + 1))); }
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    const double D = 2 + gap; const Vec3 c1(-D / 2, 0, 0), c2(D / 2, 0, 0), d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    for (int L : Ls) {
        auto t0 = std::chrono::steady_clock::now();
        TriangleMesh cube = make_cube_graded(L);
        HMatrixParams hp; hp.eps = 1e-4; hp.sep_factor = 0; hp.exact_in_lowrank = true;
        ScatteringProblem P({translated(cube, c1), translated(cube, c2)}, {Medium{eps, 1.0, 0.0}, Medium{eps, 1.0, 0.0}}, om, {}, hp);
        double tb = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        auto r0 = P.solve_plane_wave(d, px);
        FeatureSet f = FeatureSet::cube(1.0, c1); f.append(FeatureSet::cube(1.0, c2));
        P.use_block_preconditioner(group_by_features(P.mesh(), f, R));
        auto r1 = P.solve_plane_wave(d, px);
        std::printf("L=%d, N=%zu (%zu Unbekannte), Spalt %.2f: punktweise %d It. (sigma %.6f), Kanten/Ecken R=%.3f %d It. (sigma %.6f); groesster Block %zu, Anteil %.2f; Aufbau H %.0f s, Bloecke %.0f s\n",
                    L, P.mesh().size(), 8 * P.mesh().size(), gap, r0.iterations, r0.sigma_ext, R, r1.iterations, r1.sigma_ext,
                    P.block_preconditioner()->max_block(), P.block_preconditioner()->fraction(), tb, P.block_preconditioner()->seconds());
        std::fflush(stdout);
    }
    return 0;
}
