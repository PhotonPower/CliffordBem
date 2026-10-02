// Kosten des Nahfelds auf gekruemmten Elementen (direkte Summation) fuer viele Auswertepunkte und der Kraft ueber die
// Kugelflaeche; Vergleich mit dem ebenen Pfad (NearFieldOperator ab 2000 Punkten). Goldkugel in Wasser.
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -Iinclude prototype/curved/near_field_cost.cpp build/libcbem.a -o near_field_cost
//   ./near_field_cost [n]   (Ikosaederkugel mit 20 n^2 Elementen, Voreinstellung 8)
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "cbem/problems/curved_problem.hpp"
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/curved_near_field.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

int main(int argc, char** argv) {
    const int n = argc > 1 ? std::atoi(argv[1]) : 8;
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, water{1.7689, 1.0, 0.0};
    const real om = 0.5; const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    SolveOptions so; so.tol = 1e-8;
    CurvedScatteringProblem C({quadratic_icosphere(n)}, {gold}, om, water);
    ScatteringProblem F({make_icosphere(n)}, {gold}, om, water);
    const auto rc = C.solve_plane_wave(d, px, so), rf = F.solve_plane_wave(d, px, so);
    std::printf("%zu Elemente\n", C.mesh().size());
    // Karte in der Ebene y = 0, |x| <= 2, Punkte ausserhalb der Kugel (r > 1,02)
    for (int M : {1000, 10000, 40000}) {
        std::vector<Vec3> pts;
        const int g = static_cast<int>(std::ceil(std::sqrt(16.0 * M / (16 - pi))));   // Gitter g x g im Quadrat, ohne Kugel
        for (int i = 0; i < g && static_cast<int>(pts.size()) < M; ++i)
            for (int j = 0; j < g && static_cast<int>(pts.size()) < M; ++j) {
                const Vec3 x(-2 + 4.0 * (i + 0.5) / g, 0, -2 + 4.0 * (j + 0.5) / g);
                if (norm(x) > 1.02) pts.push_back(x);
            }
        auto t0 = clk::now();
        const auto fc = exterior_near_field_curved(C.mesh(), rc.h, water, om, d, px, pts);
        const double tc = since(t0);
        t0 = clk::now();
        const auto ff = exterior_near_field(F.mesh(), rf.h, water, om, d, px, pts);
        const double tf = since(t0);
        std::printf("%6zu Punkte: gekruemmt direkt %.2f s (%.1f us je Punkt und Element), eben %.2f s (%s)\n", pts.size(), tc,
                    1e6 * tc / (pts.size() * C.mesh().size()), tf, pts.size() >= 2000 ? "H-Matrix" : "direkt");
    }
    auto t0 = clk::now();
    const Vec3 Fs = force_on_sphere(make_plane_wave_eval_curved(C.mesh(), rc.h, water, om, d, px), water, Vec3(0, 0, 0), 1.5);
    std::printf("Kraft ueber die Kugel (32 x 64 Punkte): %.2f s, F_z %.6f\n", since(t0), Fs.z);
}
