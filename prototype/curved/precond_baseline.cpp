// Vorkonditionierung gekruemmter Elemente, Ausgangslage: Iterationen mit dem elementweisen Vorkonditionierer 2 (1 + J_G)^{-1}
// und ohne, fuer Glas und Gold, mehrere Netze und Toleranzen; Kosten je Iteration (Operator, Vorkonditionierer, GMRES).
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -Iinclude prototype/curved/precond_baseline.cpp build/libcbem.a -o precond_baseline
#include <chrono>
#include <cstdio>
#include <vector>

#include "cbem/problems/curved_problem.hpp"
#include "cbem/solvers/gmres.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

int main() {
    const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    struct Case { const char* name; Medium in, out; real om; };
    for (const Case& c : {Case{"Glas", Medium{2.25}, Medium{}, 1.0}, Case{"Gold", Medium{cplx(-11, 1.2)}, Medium{}, 0.5},
                          Case{"Gold/Wasser", Medium{cplx(-11, 1.2)}, Medium{1.7689}, 0.5}}) {
        for (int n : {4, 8}) {
            CurvedScatteringProblem P({quadratic_icosphere(n)}, {c.in}, c.om, c.out);
            const auto& T = P.T();
            const auto b = project_plane_wave_curved(P.mesh(), c.out.k(c.om), c.out.eps, d, px);
            LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.apply(x, y); };
            LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.precondition(x, y); };
            std::printf("%-12s %4zu Elemente:", c.name, P.mesh().size());
            for (real tol : {1e-6, 1e-10}) {
                std::vector<cplx> x1, x0;                                       // je Lauf ein eigener Startvektor (null)
                const GmresResult g1 = gmres(A, b, x1, &M, tol, 300, 3000);
                const GmresResult g0 = gmres(A, b, x0, nullptr, tol, 300, 3000);
                std::printf("  tol %.0e: %3d It. (ohne Vorkond. %3d)", tol, g1.iterations, g0.iterations);
            }
            // Kosten je Iteration
            std::vector<cplx> y;
            auto t0 = clk::now(); for (int r = 0; r < 3; ++r) T.apply(b, y); const double ta = since(t0) / 3;
            t0 = clk::now(); for (int r = 0; r < 3; ++r) T.precondition(b, y); const double tp = since(t0) / 3;
            std::vector<cplx> x; t0 = clk::now();
            const GmresResult g = gmres(A, b, x, &M, 1e-10, 300, 3000);
            const double tg = since(t0);
            std::printf("\n    je Iteration: Operator %.3f s, Vorkond. %.4f s, GMRES gesamt %.3f s (Rest %.3f s)\n", ta, tp, tg / g.iterations,
                        tg / g.iterations - ta - tp);
        }
    }
}
