// Aufwand der Nahquadratur gekruemmter Elemente ueber eine duenne Schicht (Voreinstellung curved_layered_near_params):
// Zahl der aeusseren Teilstuecke (Blaetter, je 25 Punkte), der subtrahierten Innenintegrale (je Aussenpunkt) und der Blaetter
// der Korrektur (je 25 Punkte) fuer das uebereinanderliegende Paar und fuer alle nahen Paare eines Elements; Zeit je Paar.
// Element 0 der Schale (Radius 1 + d) gegen die Elemente des Kerns (Radius 1), ein Thread.
// Bauen mit Zaehlern (aus dem Projektverzeichnis): instrumentierte Kopie von src/curved_entries.cpp, in der die Blattschleifen
// von outer und correction sowie inner_subtracted die globalen Zaehler g_outer_leaves, g_corr_leaves, g_inner_sub erhoehen:
//   g++ -std=c++17 -O3 -fopenmp -DCBEM_USE_OPENMP -fcx-fortran-rules -Iinclude prototype/curved/thin_layer_cost.cpp \
//       curved_entries_count.cpp build/libcbem.a -o thin_layer_cost
//   ./thin_layer_cost [aeusseres Kriterium] [n ...]   (Voreinstellung 1,5; n = 4 8)
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "cbem/problems/curved_layered_problem.hpp"

namespace cbem { extern long g_outer_leaves, g_inner_sub, g_corr_leaves; }
using namespace cbem;
using clk = std::chrono::steady_clock;

int main(int argc, char** argv) {
    const real outer = argc > 1 ? std::atof(argv[1]) : 1.5;
    std::vector<int> ns; for (int i = 2; i < argc; ++i) ns.push_back(std::atoi(argv[i]));
    if (ns.empty()) ns = {4, 8};
    const cplx k = 0.5 * std::sqrt(cplx(2.25));
    EntryParams ep; ep.cache_near = false;
    for (int n : ns)
        for (real d : {0.05, 0.02, 0.01, 0.005}) {
            const QuadraticMesh m = merge_quadratic({quadratic_icosphere(n, 1 + d), quadratic_icosphere(n, 1.0)});
            const std::size_t N = m.size() / 2;
            CurvedNearParams np = curved_layered_near_params(); np.subtract_outer_ratio = outer;
            const CurvedKernelEntries E(m, k, ep, np);
            std::printf("n = %d, d = %.3f (d/h = %.3f):\n", n, d, d / m.flat.hmax[0]);
            auto run = [&](const char* name, const std::vector<std::size_t>& js) {
                g_outer_leaves = g_inner_sub = g_corr_leaves = 0;
                auto t0 = clk::now();
                for (std::size_t j : js) E.lambda_near(0, j);
                const double t = std::chrono::duration<double>(clk::now() - t0).count();
                const double np = double(js.size());
                std::printf("  %-22s %3zu Paare: je Paar %8.1f aeussere Blaetter, %8.0f Aussenpunkte, %6.1f Korrekturblaetter je Aussenpunkt; %.2f ms je Paar\n",
                            name, js.size(), g_outer_leaves / np, g_inner_sub / np, double(g_corr_leaves) / std::max(1L, g_inner_sub), 1e3 * t / np);
            };
            run("uebereinander", {N});
            std::vector<std::size_t> js;
            for (std::size_t j = N; j < 2 * N; ++j) if (j != N && E.is_near(0, j)) js.push_back(j);
            run("uebrige nahe Paare", js);
            std::fflush(stdout);
        }
}
