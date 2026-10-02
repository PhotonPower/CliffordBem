// Kosten und Genauigkeit der Nahquadratur gekruemmter Elemente fuer Paare ueber eine duenne Schicht: Element 0 der
// Schalenflaeche (Radius 1 + d) gegen alle nahen Elemente der Kernflaeche (Radius 1), ein Thread. Varianten der Kriterien
// (CurvedNearParams): Voreinstellung, Randabstand (adapt_to_boundary), lockere Korrektur, lockere Aussenregel. Abweichung
// bezogen auf den groessten Eintrag, gegen die Voreinstellung.
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -DCBEM_USE_OPENMP -Iinclude prototype/curved/layered_pairs.cpp build/libcbem.a -o layered_pairs
//   ./layered_pairs [n]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "cbem/assembly/curved_entries.hpp"
#include "cbem/geometry/quadratic_mesh.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;

int main(int argc, char** argv) {
    const int n = argc > 1 ? std::atoi(argv[1]) : 4;
    const cplx k = 0.5 * std::sqrt(cplx(2.25));
    EntryParams ep; ep.cache_near = false;
    struct Var { const char* name; bool adapt; real outer, corr; };
    for (real d : {0.2, 0.05, 0.02, 0.01}) {
        const QuadraticMesh m = merge_quadratic({quadratic_icosphere(n, 1 + d), quadratic_icosphere(n, 1.0)});
        const std::size_t N = m.size() / 2;
        std::vector<std::size_t> js;
        {
            const CurvedKernelEntries E0(m, k, ep);
            for (std::size_t j = N; j < 2 * N; ++j) if (E0.is_near(0, j)) js.push_back(j);
        }
        std::printf("d = %.2f (d/h = %.2f): %zu nahe Paare ueber die Schicht\n", d, d / m.flat.hmax[0], js.size());
        std::vector<CurvedBlock> ref;
        for (const Var& v : {Var{"Voreinstellung", false, 1.5, 1.5}, Var{"Randabstand", true, 1.5, 1.5}, Var{"Randabstand, Korrektur 3", true, 1.5, 3.0},
                             Var{"Randabstand, Aussen 3", true, 3.0, 1.5}}) {
            if (!v.adapt && d < 0.015) { std::printf("  %-26s (uebersprungen)\n", v.name); continue; }
            CurvedNearParams np; np.adapt_to_boundary = v.adapt; np.subtract_outer_ratio = v.outer; np.correction_ratio = v.corr;
            const CurvedKernelEntries E(m, k, ep, np);
            std::vector<CurvedBlock> K(js.size());
            auto t0 = clk::now();
            for (std::size_t a = 0; a < js.size(); ++a) K[a] = E.lambda_near(0, js[a]);
            const double t = std::chrono::duration<double>(clk::now() - t0).count();
            if (ref.empty()) ref = K;
            real e = 0, s = 0;
            for (std::size_t a = 0; a < js.size(); ++a)
                for (int q = 0; q < 9; ++q) for (int c = 0; c < kCurvedComps; ++c) { e = std::max(e, std::abs(K[a][q][c] - ref[a][q][c])); s = std::max(s, std::abs(ref[a][q][c])); }
            const double tstack = [&] { auto t1 = clk::now(); E.lambda_near(0, N); return std::chrono::duration<double>(clk::now() - t1).count(); }();
            std::printf("  %-26s %.3f s je Element (%.1f ms je Paar, uebereinander %.1f ms), Abweichung %.1e\n", v.name, t, 1e3 * t / js.size(), 1e3 * tstack, e / s);
            std::fflush(stdout);
        }
    }
}
