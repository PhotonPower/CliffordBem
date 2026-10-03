// Kriterien und Blattregeln der Nahquadratur ueber eine duenne Schicht (adapt_to_boundary): Zeit je Element (Element 0 der
// Schale gegen alle nahen Elemente des Kerns, ein Thread) und Fehler gegen eine verschaerfte Referenz (Kriterien 0,6/0,6,
// Gauss 6 x 6), je Paar bezogen auf dessen groessten Eintrag (max ueber die Paare) und auf den groessten Eintrag aller Paare.
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -DCBEM_USE_OPENMP -Iinclude prototype/curved/thin_layer_rules.cpp build/libcbem.a -o thin_layer_rules
//   ./thin_layer_rules [n d]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "cbem/problems/curved_layered_problem.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;

int main(int argc, char** argv) {
    struct Geo { int n; real d; };
    std::vector<Geo> geos;
    if (argc > 2) geos.push_back({std::atoi(argv[1]), std::atof(argv[2])});
    else geos = {{4, 0.02}, {4, 0.01}, {8, 0.01}};
    const cplx k = 0.5 * std::sqrt(cplx(2.25));
    EntryParams ep; ep.cache_near = false;
    struct Var { real outer, corr; int ro, rc; };
    for (const Geo& G : geos) {
        const QuadraticMesh m = merge_quadratic({quadratic_icosphere(G.n, 1 + G.d), quadratic_icosphere(G.n, 1.0)});
        const std::size_t N = m.size() / 2;
        std::vector<std::size_t> js;
        { const CurvedKernelEntries E(m, k, ep); for (std::size_t j = N; j < 2 * N; ++j) if (E.is_near(0, j)) js.push_back(j); }
        auto compute = [&](const Var& v, double& t) {
            CurvedNearParams np = curved_layered_near_params();
            np.subtract_outer_ratio = v.outer; np.correction_ratio = v.corr; np.outer_rule = v.ro; np.correction_rule = v.rc;
            const CurvedKernelEntries E(m, k, ep, np);
            std::vector<CurvedBlock> K(js.size());
            auto t0 = clk::now();
            for (std::size_t a = 0; a < js.size(); ++a) K[a] = E.lambda_near(0, js[a]);
            t = std::chrono::duration<double>(clk::now() - t0).count();
            return K;
        };
        double tref;
        const auto R = compute({0.6, 0.6, 6, 6}, tref);
        real gmax = 0;
        for (auto& B : R) for (auto& c : B) for (auto& v : c) gmax = std::max(gmax, std::abs(v));
        std::printf("n = %d, d = %.3f (d/h = %.3f), %zu Paare; Referenz %.1f s\n", G.n, G.d, G.d / m.flat.hmax[0], js.size(), tref);
        for (const Var& v : {Var{1.5, 1.5, 5, 5}, Var{3, 1.5, 5, 5}, Var{5, 1.5, 5, 5}, Var{8, 1.5, 5, 5}, Var{1.5, 3, 5, 5}, Var{3, 3, 5, 5},
                             Var{5, 3, 5, 5}, Var{3, 2, 5, 5}, Var{5, 2, 5, 5}, Var{3, 3, 4, 4}, Var{5, 3, 4, 5}, Var{5, 3, 5, 4}}) {
            double t;
            const auto K = compute(v, t);
            real ep_ = 0, eg = 0;
            for (std::size_t a = 0; a < K.size(); ++a) {
                real e = 0, s = 0;
                for (int q = 0; q < 9; ++q) for (int c = 0; c < kCurvedComps; ++c) { e = std::max(e, std::abs(K[a][q][c] - R[a][q][c])); s = std::max(s, std::abs(R[a][q][c])); }
                ep_ = std::max(ep_, e / s); eg = std::max(eg, e / gmax);
            }
            std::printf("  aussen %.1f, Korrektur %.1f, Regeln %d/%d: %.3f s je Element, Fehler je Paar %.1e, global %.1e\n", v.outer, v.corr, v.ro, v.rc, t, ep_, eg);
            std::fflush(stdout);
        }
    }
}
