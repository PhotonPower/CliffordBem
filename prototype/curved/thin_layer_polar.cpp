// Korrektur der Singularitaetssubtraktion in Polarkoordinaten um den Fusspunkt (CurvedNearParams::polar_radial/_angular)
// gegen die adaptive Korrektur: Zeit und Fehler gegen eine verschaerfte Referenz (Kriterien 0,6/0,6, Gauss 6 x 6, adaptiv).
//  (a) duenne Schicht: Element 0 der Schale (Radius 1 + d) gegen alle nahen Elemente des Kerns (aeusseres Kriterium 3);
//  (b) eine Flaeche: Element 0 gegen alle nahen, nicht benachbarten Elemente derselben Kugel;
//  (c) Nahfeld: point_integrals an Punkten im Abstand 0,003 ... 0,3 ueber Element 0 (Schwerpunkt und nahe einer Kante).
// Fehler je Paar bzw. Punkt bezogen auf dessen groessten Eintrag, Maximum ueber Paare/Punkte; ein Thread.
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -DCBEM_USE_OPENMP -Iinclude prototype/curved/thin_layer_polar.cpp build/libcbem.a -o thin_layer_polar
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "cbem/problems/curved_layered_problem.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;

struct Var { const char* name; real outer, corr; int ro, rc, nr, na; real below; };

static CurvedNearParams params(const Var& v) {
    CurvedNearParams np = curved_layered_near_params();
    np.subtract_outer_ratio = v.outer; np.correction_ratio = v.corr; np.outer_rule = v.ro; np.correction_rule = v.rc;
    np.polar_radial = v.nr; np.polar_angular = v.na; np.polar_below = v.below;
    return np;
}

template <class F> static void compare(const char* title, const std::vector<Var>& vars, F&& compute) {
    std::printf("%s\n", title);
    double tref;
    const auto R = compute(params(vars[0]), tref);
    std::printf("  %-28s %.3f s\n", vars[0].name, tref);
    for (std::size_t i = 1; i < vars.size(); ++i) {
        double t;
        const auto K = compute(params(vars[i]), t);
        real e = 0;
        for (std::size_t a = 0; a < K.size(); ++a) {
            real ea = 0, sa = 0;
            for (std::size_t q = 0; q < K[a].size(); ++q) { ea = std::max(ea, std::abs(K[a][q] - R[a][q])); sa = std::max(sa, std::abs(R[a][q])); }
            e = std::max(e, ea / sa);
        }
        std::printf("  %-28s %.3f s, Fehler %.1e\n", vars[i].name, t, e);
        std::fflush(stdout);
    }
}

int main() {
    const cplx k = 0.5 * std::sqrt(cplx(2.25));
    EntryParams ep; ep.cache_near = false;
    const real all = 1e300;
    const std::vector<Var> vars = {{"Referenz", 0.6, 0.6, 6, 6, 0, 0, all},
                                   {"adaptiv (Voreinstellung)", 3, 1.5, 5, 5, 0, 0, all},
                                   {"polar 6 x 4", 3, 1.5, 5, 5, 6, 4, all},
                                   {"polar 8 x 6", 3, 1.5, 5, 5, 8, 6, all},
                                   {"polar 10 x 8", 3, 1.5, 5, 5, 10, 8, all},
                                   {"polar 12 x 10", 3, 1.5, 5, 5, 12, 10, all},
                                   {"polar 8 x 6 unter 0,3 h", 3, 1.5, 5, 5, 8, 6, 0.3},
                                   {"polar 10 x 8 unter 0,3 h", 3, 1.5, 5, 5, 10, 8, 0.3},
                                   {"polar 10 x 8 unter 0,1 h", 3, 1.5, 5, 5, 10, 8, 0.1}};
    auto flat = [](const CurvedBlock& B) { std::vector<cplx> v; for (auto& c : B) for (auto& x : c) v.push_back(x); return v; };
    // (a) duenne Schicht
    for (auto [n, d] : {std::pair<int, real>{4, 0.02}, {4, 0.01}, {8, 0.01}, {8, 0.002}}) {
        const QuadraticMesh m = merge_quadratic({quadratic_icosphere(n, 1 + d), quadratic_icosphere(n, 1.0)});
        const std::size_t N = m.size() / 2;
        std::vector<std::size_t> js;
        { const CurvedKernelEntries E(m, k, ep); for (std::size_t j = N; j < 2 * N; ++j) if (E.is_near(0, j)) js.push_back(j); }
        char title[128]; std::snprintf(title, sizeof title, "(a) Schicht n = %d, d = %.3f (d/h = %.3f), %zu Paare", n, d, d / m.flat.hmax[0], js.size());
        compare(title, vars, [&](const CurvedNearParams& np, double& t) {
            const CurvedKernelEntries E(m, k, ep, np);
            std::vector<std::vector<cplx>> K;
            auto t0 = clk::now();
            for (std::size_t j : js) K.push_back(flat(E.lambda_near(0, j)));
            t = std::chrono::duration<double>(clk::now() - t0).count();
            return K;
        });
    }
    // (b) eine Flaeche
    {
        const QuadraticMesh m = quadratic_icosphere(4);
        std::vector<std::size_t> js;
        { const CurvedKernelEntries E(m, k, ep); for (std::size_t j = 1; j < m.size(); ++j) if (E.is_near(0, j) && E.adjacency(0, j) == Adjacency::None) js.push_back(j); }
        char title[128]; std::snprintf(title, sizeof title, "(b) eine Flaeche n = 4, %zu nahe, nicht benachbarte Paare", js.size());
        compare(title, vars, [&](const CurvedNearParams& np, double& t) {
            const CurvedKernelEntries E(m, k, ep, np);
            std::vector<std::vector<cplx>> K;
            auto t0 = clk::now();
            for (std::size_t j : js) K.push_back(flat(E.lambda_near(0, j)));
            t = std::chrono::duration<double>(clk::now() - t0).count();
            return K;
        });
    }
    // (c) Nahfeld: Punkte ueber Element 0
    {
        const QuadraticMesh m = quadratic_icosphere(4);
        std::vector<Vec3> pts;
        const Vec3 c = m.X(0, {1.0 / 3, 1.0 / 3, 1.0 / 3}), e = m.X(0, {0.48, 0.48, 0.04}), f = m.X(0, {0.495, 0.495, 0.01}), g = m.X(0, {0.9, 0.05, 0.05});
        for (real h : {0.001, 0.003, 0.01, 0.03, 0.1, 0.3}) for (const Vec3& p : {c, e, f, g}) pts.push_back(p * (1 + h));
        compare("(c) Nahfeld: 24 Punkte im Abstand 0,001 ... 0,3 ueber Element 0 (Mitte, nahe Kante, an Kante, nahe Ecke; n = 4)", vars, [&](const CurvedNearParams& np, double& t) {
            const CurvedKernelEntries E(m, k, ep, np);
            std::vector<std::vector<cplx>> K;
            auto t0 = clk::now();
            for (const Vec3& x : pts) { const auto G = E.point_integrals(x, 0); std::vector<cplx> v; for (auto& g : G) for (auto& z : g) v.push_back(z); K.push_back(v); }
            t = std::chrono::duration<double>(clk::now() - t0).count();
            return K;
        });
    }
}
