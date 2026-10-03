// Aufbau und Genauigkeit von CurvedLayeredScatteringProblem bei duennen Schichten mit verschiedenen Nahquadratur-Parametern:
// Goldkern (Radius 1) mit Glasschale der Dicke d, omega a = 0,5, gegen Aden-Kerker (tools/mie_coated.py); 12 Threads.
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -DCBEM_USE_OPENMP -Iinclude prototype/curved/layered_build.cpp build/libcbem.a -o layered_build
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "cbem/problems/curved_layered_problem.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;

int main() {
    const real om = 0.5;
    const Medium gold{cplx(-11, 1.2)}, glass{2.25};
    SolveOptions so; so.tol = 1e-10;
    struct Case { int n; real d, sig; };
    // sigma = Q pi (1 + d)^2 mit Q aus tools/mie_coated.py
    for (const Case& c : {Case{4, 0.01, 0.60647835 * pi * 1.01 * 1.01},
                          Case{8, 0.01, 0.60647835 * pi * 1.01 * 1.01}}) {
        CurvedLayeredGeometry g; add_layered_body(g, {quadratic_icosphere(c.n, 1 + c.d), quadratic_icosphere(c.n, 1.0)}, {glass, gold});
        CurvedNearParams v12 = curved_layered_near_params(); v12.subtract_outer_ratio = 1.5;
        CurvedNearParams pol = curved_layered_near_params(); pol.polar_radial = 0;
        std::printf("n = %d, d = %.2f:\n", c.n, c.d);
        // aeusseres Kriterium 1,5 (v0.62), 3 (v0.63 A), 3 mit polarer Korrektur (v0.63 A + B, Voreinstellung)
        for (const CurvedNearParams& np : {v12, pol, curved_layered_near_params()}) {
            auto t0 = clk::now();
            CurvedLayeredScatteringProblem P(g, om, curved_hmatrix_params(), EntryParams{}, np);
            const double tb = std::chrono::duration<double>(clk::now() - t0).count();
            const auto r = P.solve_plane_wave(Vec3(0, 0, 1), CVec3{1.0, 0.0, 0.0}, so);
            std::printf("  aeusseres Kriterium %.1f, polar %d x %d: Aufbau %.1f s (Nahfeld %.1f s), sigma %.9f, gegen Aden-Kerker %+.2e\n",
                        np.subtract_outer_ratio, np.polar_radial, np.polar_angular, tb, P.near_seconds(), r.sigma_ext, r.sigma_ext / c.sig - 1);
            std::fflush(stdout);
        }
    }
}
