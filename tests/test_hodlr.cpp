// Hierarchische Faktorisierung (HODLR): Einzeleintraege gegen den Operator, feine Toleranz = direkter Loeser,
// grobe Toleranz als Vorkonditionierer (gleiche Loesung, wenige Iterationen); Mehrkoerper mit chiralem Medium.
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/fields.hpp"
#include <random>
#include "check.hpp"
using namespace cbem;
int main() {
    TriangleMesh s2 = make_icosphere(3); HMatrixParams hp; hp.eps = 1e-10; SolveOptions so; so.tol = 1e-10;
    ScatteringProblem P({translated(s2, Vec3(-1.3, 0, 0)), translated(s2, Vec3(1.3, 0, 0))},
                        {Medium{cplx(-4, 0.3), 1.0, 0.0}, Medium{2.25, 1.0, 0.15}}, 0.7, {}, hp);
    const std::size_t N = P.mesh().size();
    // (1) system_entry gegen T.apply auf Einheitsvektoren
    std::mt19937 g(2); std::normal_distribution<real> nd; std::vector<cplx> x(8 * N), y, yref(8 * N, 0.0);
    for (auto& v : x) v = cplx(nd(g), nd(g));
    P.T().apply(x, y);
    real e1 = 0, n1 = 0;
    for (std::size_t i = 0; i < N; i += 37) {
        for (std::size_t j = 0; j < N; ++j) { Mat8 T = P.system_entry(i, j); for (int r = 0; r < 8; ++r) for (int q = 0; q < 8; ++q) yref[8 * i + r] += T[r * 8 + q] * x[8 * j + q]; }
        for (int r = 0; r < 8; ++r) { e1 += std::norm(yref[8 * i + r] - y[8 * i + r]); n1 += std::norm(y[8 * i + r]); }
    }
    std::printf("  Einzeleintraege gegen Operator: rel. Abw. %.1e\n", std::sqrt(e1 / n1));
    CHECK(std::sqrt(e1 / n1) < 1e-7, "system_entry passt nicht zum Operator");
    const Vec3 d = Vec3(0.3, 0.4, 0.866) / norm(Vec3(0.3, 0.4, 0.866)); const CVec3 pol = circular_polarization(d, +1);
    auto r0 = P.solve_plane_wave(d, pol, so);
    for (real eps : {1e-9, 1e-2}) {
        HodlrParams hpar; hpar.leaf = 32; hpar.eps = eps; P.use_hodlr_preconditioner(hpar);
        auto r = P.solve_plane_wave(d, pol, so);
        std::printf("  HODLR eps %.0e: GMRES %d It. (punktweise %d), max. Rang %zu, %.1f MB, sigma %.10f / %.10f\n", eps, r.iterations, r0.iterations,
                    P.hodlr()->max_rank(), P.hodlr()->bytes() / 1048576.0, r.sigma_ext, r0.sigma_ext);
        CHECK(std::abs(r.sigma_ext - r0.sigma_ext) < 1e-8 * std::abs(r0.sigma_ext), "Loesung weicht ab");
        if (eps < 1e-6) CHECK(r.iterations <= 3, "feine Toleranz sollte direkt loesen");
        else CHECK(r.iterations < r0.iterations / 2, "grobe Toleranz sollte Iterationen stark senken");
    }
    REPORT();
}
