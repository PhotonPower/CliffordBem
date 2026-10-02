// Krylov-Recycling fuer gekruemmte Elemente: Iterationen und Zeit fuer eine Folge ebener Wellen (6 Richtungen x 2
// Polarisationen, wie bei der Orientierungsmittelung) mit GMRES, RecyclingGmres und GCRO-DR; Goldkugel bei eps = -11 + 1,2i
// und nahe der Dipolresonanz.
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -Iinclude prototype/curved/recycling_rhs.cpp build/libcbem.a -o recycling_rhs
//   ./recycling_rhs [n] [tol]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>

#include "cbem/problems/curved_problem.hpp"
#include "cbem/solvers/recycling_gmres.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;

int main(int argc, char** argv) {
    const int n = argc > 1 ? std::atoi(argv[1]) : 4;
    const real tol = argc > 2 ? std::atof(argv[2]) : 1e-6;
    // Richtungen: Oktaeder (6), je zwei Polarisationen senkrecht dazu
    std::vector<std::pair<Vec3, CVec3>> rhs;
    for (const Vec3& d : {Vec3(0, 0, 1), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0.6, 0.0, 0.8), Vec3(0.0, -0.6, 0.8), Vec3(0.48, 0.6, -0.64)}) {
        const Vec3 dn = d / norm(d);
        Vec3 u = std::abs(dn.z) < 0.9 ? cross(dn, Vec3(0, 0, 1)) : cross(dn, Vec3(1, 0, 0)); u = u / norm(u);
        const Vec3 v = cross(dn, u);
        rhs.push_back({dn, CVec3{u.x, u.y, u.z}});
        rhs.push_back({dn, CVec3{v.x, v.y, v.z}});
    }
    struct Case { const char* name; Medium in; real om; };
    for (const Case& c : {Case{"Gold eps -11+1,2i", Medium{cplx(-11, 1.2)}, 0.5}, Case{"nahe Resonanz eps -2,3+0,2i", Medium{cplx(-2.3, 0.2)}, 0.5}}) {
        CurvedScatteringProblem P({quadratic_icosphere(n)}, {c.in}, c.om);
        const auto& T = P.T();
        LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.apply(x, y); };
        LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.precondition(x, y); };
        std::vector<std::vector<cplx>> B;
        for (auto& [d, p] : rhs) B.push_back(project_plane_wave_curved(P.mesh(), c.om, 1.0, d, p));
        std::printf("%s, %zu Elemente, tol %.0e\n", c.name, P.mesh().size(), tol);
        std::vector<std::vector<cplx>> X0(B.size());
        for (int v = 0; v < 3; ++v) {
            std::unique_ptr<MultiRhsSolver> S;
            if (v == 1) S = std::make_unique<RecyclingGmres>(A, &M, 120);
            if (v == 2) S = std::make_unique<GcroDr>(A, &M, 20, 80);
            auto t0 = clk::now(); int tot = 0; double dmax = 0;
            std::printf("  %-16s", v == 0 ? "GMRES" : (v == 1 ? "Recycling 120" : "GCRO-DR 20/80"));
            for (std::size_t r = 0; r < B.size(); ++r) {
                std::vector<cplx> x;
                const GmresResult g = v == 0 ? gmres(A, B[r], x, &M, tol, 300, 3000) : S->solve(B[r], x, tol, 300, 3000);
                tot += g.iterations; std::printf(" %d", g.iterations);
                if (v == 0) X0[r] = x;
                else { double nn = 0, dd = 0; for (std::size_t i = 0; i < x.size(); ++i) { nn += std::norm(X0[r][i]); dd += std::norm(x[i] - X0[r][i]); } dmax = std::max(dmax, std::sqrt(dd / nn)); }
            }
            std::printf("  | gesamt %d It., %.1f s%s\n", tot, std::chrono::duration<double>(clk::now() - t0).count(),
                        v ? (" , Abweichung von GMRES " + std::to_string(dmax).substr(0, 8)).c_str() : "");
        }
    }
}
