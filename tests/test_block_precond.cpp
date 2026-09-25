// Block-Vorkonditionierer: gleiche Loesung wie punktweise, weniger Iterationen; volle Abdeckung = exakte Inverse.
#include "cbem/solvers/block_preconditioner.hpp"
#include "cbem/solvers/gmres.hpp"
#include "cbem/sources/fields.hpp"
#include "cbem/problems/scattering_problem.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    TriangleMesh m = make_cube_uniform(3);
    Medium in{cplx(-4, 0.3), 1.0}, out{1.0, 1.0}; const real om = 0.5; const cplx k1 = om * std::sqrt(in.eps), k2 = om;
    KernelEntries E1e(m, k1), E2e(m, k2); HMatrixParams p; p.eps = 1e-10; p.sep_factor = 0; p.exact_in_lowrank = true;
    KernelHMatrix H1(E1e, p), H2(E2e, p); CauchyOperator E1(m, H1), E2(m, H2); TransmissionOperator T(m, E1, E2, in, out);
    auto b = project_plane_wave(m, k2, 1.0, Vec3(0, 0, 1), Vec3(1, 0, 0));
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.apply(x, y); };
    LinOp Mp = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.precondition(x, y); };
    std::vector<cplx> h0; GmresResult r0 = gmres(A, b, h0, &Mp, 1e-10, 300, 2000);
    auto G = group_by_features(m, FeatureSet::cube(), 0.6); BlockPreconditioner P(m, E1e, E2e, T, G);
    LinOp Mb = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { P.apply(x, y); };
    std::vector<cplx> h1; GmresResult r1 = gmres(A, b, h1, &Mb, 1e-10, 300, 2000);
    real e = 0, s = 0; for (std::size_t i = 0; i < h0.size(); ++i) { e += std::norm(h1[i] - h0[i]); s += std::norm(h0[i]); }
    std::printf("  punktweise %d It., Bloecke %d It. (Anteil %.2f), Loesungsdifferenz %.1e\n", r0.iterations, r1.iterations, P.fraction(), std::sqrt(e / s));
    CHECK(r1.iterations < r0.iterations, "Bloecke sollten weniger Iterationen brauchen");
    CHECK(std::sqrt(e / s) < 1e-8, "Loesungen verschieden");
    // ein einziger Block ueber alles = exakte Inverse von T (bis auf H-Fehler 1e-10): 1-2 Iterationen
    std::vector<std::vector<std::size_t>> all(1); for (std::size_t t = 0; t < m.size(); ++t) all[0].push_back(t);
    BlockPreconditioner Pf(m, E1e, E2e, T, all);
    LinOp Mf = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { Pf.apply(x, y); };
    std::vector<cplx> h2; GmresResult r2 = gmres(A, b, h2, &Mf, 1e-10, 300, 2000);
    std::printf("  voller Block: %d It.\n", r2.iterations);
    CHECK(r2.iterations <= 3, "voller Block sollte exakt invertieren");
    // Mehrere Koerper, einer chiral: innerer Block aus ScatteringProblem::inner_block
    {
        TriangleMesh s2 = make_icosphere(2); HMatrixParams hp; hp.eps = 1e-10; SolveOptions so; so.tol = 1e-10;
        ScatteringProblem P({translated(s2, Vec3(-1.3, 0, 0)), translated(s2, Vec3(1.3, 0, 0))},
                            {Medium{cplx(-4, 0.3), 1.0, 0.0}, Medium{2.25, 1.0, 0.15}}, 0.7, {}, hp);
        const Vec3 d(0.3, 0.4, 0.866); const CVec3 pol = circular_polarization(d / norm(d), +1);
        auto r0 = P.solve_plane_wave(d / norm(d), pol, so);
        P.use_block_preconditioner(group_by_clusters(P.mesh(), 32)); auto r1 = P.solve_plane_wave(d / norm(d), pol, so);
        std::vector<std::vector<std::size_t>> all(1); for (std::size_t t = 0; t < P.mesh().size(); ++t) all[0].push_back(t);
        P.use_block_preconditioner(all); auto r2 = P.solve_plane_wave(d / norm(d), pol, so);
        std::printf("  Dimer (chiral + plasmonisch): punktweise %d It., Cluster (32) %d It., voller Block %d It.; sigma %.10f / %.10f / %.10f\n",
                    r0.iterations, r1.iterations, r2.iterations, r0.sigma_ext, r1.sigma_ext, r2.sigma_ext);
        CHECK(r2.iterations <= 3, "voller Block (Mehrkoerper, chiral) sollte exakt invertieren");
        CHECK(r1.iterations < r0.iterations, "Clusterbloecke sollten Iterationen sparen");
        CHECK(std::abs(r1.sigma_ext - r0.sigma_ext) < 1e-8 * std::abs(r0.sigma_ext), "Loesungen verschieden");
    }
    REPORT();
}
