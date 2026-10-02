// Geschichtete Koerper auf gekruemmten Elementen, Messung vor dem Einbau: Goldkern (Radius 1, eps = -11 + 1,2i) mit
// Glasschale (eps = 2,25) bzw. neutraler Schale (Aussenmedium) der Dicke d, omega a = 0,5, gegen Aden-Kerker
// (tools/mie_coated.py). System wie LayeredTransmissionOperator (layered_problem.hpp), aus den Bausteinen der gekruemmten
// Elemente zusammengesetzt: je Gebiet R ein CurvedCauchyOperator auf seinem Rand (Flaechen mit ihren eigenen
// Aussennormalen, Vorzeichen sigma in der Dichte), je Flaeche J_G = Galerkin-Projektion der Transmissionsabbildung;
// Zeile s: 1/2 (h_s - [E_out(sigma u)]_s) + 1/2 (J_s h_s - [E_in(sigma u)]_s) = b_s. Vergleich mit LayeredScatteringProblem
// (konstante Dichten, ebene Dreiecke) auf denselben Ikosaedernetzen.
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -DCBEM_USE_OPENMP -Iinclude prototype/curved/layered.cpp build/libcbem.a -o layered
//   ./layered [Randabstand 0/1] [n ...]   (Voreinstellung 1, n = 4 6 8; Randabstand: CurvedNearParams::adapt_to_boundary)
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>

#include "cbem/problems/curved_problem.hpp"
#include "cbem/problems/layered_problem.hpp"
#include "cbem/solvers/gmres.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

struct NullOp : BoundaryOperator {
    std::size_t n;
    explicit NullOp(std::size_t n_) : n(n_) {}
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const override { y.assign(x.size(), cplx(0)); }
    std::size_t size() const override { return n; }
};

struct CurvedRun { real sigma; int it; double build, solve; };

// Kern (Radius 1) mit Schale bis 1 + d; zwei Flaechen, drei Gebiete (0 aussen, 1 Schale, 2 Kern)
static CurvedRun curved_coated(int n, real d, const Medium& shell, const Medium& core, real om, bool adapt) {
    const Medium vac{};
    auto t0 = clk::now();
    const QuadraticMesh s0 = quadratic_icosphere(n, 1 + d), s1 = quadratic_icosphere(n, 1.0);
    const std::size_t N0 = s0.size(), N1 = s1.size();
    const QuadraticMesh shellmesh = merge_quadratic({s0, s1});
    // Gebietsoperatoren
    struct Reg { std::unique_ptr<CurvedKernelEntries> E; std::unique_ptr<CurvedHMatrix> H; std::unique_ptr<CurvedCauchyOperator> C; };
    auto make = [&](const QuadraticMesh& m, const Medium& med) {
        CurvedNearParams np; np.adapt_to_boundary = adapt;
        Reg r; r.E = std::make_unique<CurvedKernelEntries>(m, med.k(om), EntryParams{}, np); r.H = std::make_unique<CurvedHMatrix>(*r.E, curved_hmatrix_params());
        r.C = std::make_unique<CurvedCauchyOperator>(m, *r.H); return r;
    };
    Reg R0 = make(s0, vac), R1 = make(shellmesh, shell), R2 = make(s1, core);
    // Transmissionsabbildungen je Flaeche (innen/aussen) und Vorkonditionierer 2 (1 + J_G)^{-1}
    const NullOp z0(24 * N0), z1(24 * N1);
    const CurvedTransmissionOperator J0(s0, curved_psi_matrices(s0), z0, z0, std::vector<Medium>(N0, shell), vac);
    const CurvedTransmissionOperator J1(s1, curved_psi_matrices(s1), z1, z1, std::vector<Medium>(N1, core), shell);
    const double tb = since(t0);
    const std::size_t A = 24 * N0, B = 24 * N1;
    LinOp T = [&](const std::vector<cplx>& x, std::vector<cplx>& y) {
        std::vector<cplx> h0(x.begin(), x.begin() + A), h1(x.begin() + A, x.end()), Jh0, Jh1, w0, w1, w2;
        J0.apply_J(h0, Jh0); J1.apply_J(h1, Jh1);
        std::vector<cplx> m0(A); for (std::size_t i = 0; i < A; ++i) m0[i] = -h0[i];          // Aussenraum: sigma = -1 an Flaeche 0
        std::vector<cplx> m1(A + B);                                                         // Schale: +J0 h0 an 0, -h1 an 1
        for (std::size_t i = 0; i < A; ++i) m1[i] = Jh0[i];
        for (std::size_t i = 0; i < B; ++i) m1[A + i] = -h1[i];
        R0.C->apply(m0, w0); R1.C->apply(m1, w1); R2.C->apply(Jh1, w2);                      // Kern: +J1 h1 an 1
        y.assign(A + B, cplx(0));
        for (std::size_t i = 0; i < A; ++i) y[i] = 0.5 * (h0[i] - w0[i]) + 0.5 * (Jh0[i] - w1[i]);
        for (std::size_t i = 0; i < B; ++i) y[A + i] = 0.5 * (h1[i] - w1[A + i]) + 0.5 * (Jh1[i] - w2[i]);
    };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) {
        std::vector<cplx> a(x.begin(), x.begin() + A), b(x.begin() + A, x.end()), pa, pb;
        J0.precondition(a, pa); J1.precondition(b, pb);
        y = pa; y.insert(y.end(), pb.begin(), pb.end());
    };
    const Vec3 dd(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    const cplx k0 = vac.k(om);
    std::vector<cplx> b = project_plane_wave_curved(s0, k0, vac.eps, dd, px); b.resize(A + B, cplx(0));
    t0 = clk::now();
    std::vector<cplx> h;
    const GmresResult g = gmres(T, b, h, &M, 1e-8, 300, 3000);
    const double ts = since(t0);
    std::vector<cplx> hs(h.begin(), h.begin() + A); for (std::size_t i = 0; i < A; ++i) hs[i] -= b[i];
    return {extinction_cross_section_curved(s0, hs, k0, vac.eps, dd, px), g.iterations, tb, ts};
}

static CurvedRun const_coated(int n, real d, const Medium& shell, const Medium& core, real om) {
    auto t0 = clk::now();
    LayeredGeometry g;
    add_layered_body(g, {make_icosphere(n, 1 + d), make_icosphere(n, 1.0)}, {shell, core});
    HMatrixParams hp; hp.eps = 1e-6;
    LayeredScatteringProblem P(g, om, hp);
    const double tb = since(t0);
    SolveOptions so; so.tol = 1e-8;
    t0 = clk::now();
    const auto r = P.solve_plane_wave(Vec3(0, 0, 1), CVec3{1.0, 0.0, 0.0}, so);
    return {r.sigma_ext, r.iterations, tb, since(t0)};
}

int main(int argc, char** argv) {
    const bool adapt = argc > 1 ? std::atoi(argv[1]) != 0 : true;
    std::vector<int> ns;
    for (int i = 2; i < argc; ++i) ns.push_back(std::atoi(argv[i]));
    if (ns.empty()) ns = {4, 6, 8};
    const real om = 0.5;
    const Medium gold{cplx(-11, 1.2)}, glass{2.25}, vac{};
    // Aden-Kerker (tools/mie_coated.py), Q_ext bezogen auf pi (1 + d)^2; ohne Schicht: Q = 0,53516393 bezogen auf pi 1,05^2
    const real sig_bare = 0.53516393 * pi * 1.05 * 1.05;
    struct Ref { real d, q; };
    for (const Ref& rf : {Ref{0.2, 0.95996079}, Ref{0.05, 0.67462135}, Ref{0.02, 0.62317153}, Ref{0.01, 0.60647835}}) {
        const real sig = rf.q * pi * (1 + rf.d) * (1 + rf.d), dsig = sig - sig_bare;
        std::printf("d = %.2f: sigma (Aden-Kerker) %.6f, ohne Schicht %.6f, Schichtwirkung %.6f\n", rf.d, sig, sig_bare, dsig);
        for (int n : ns) {
            for (int curved : {1, 0}) {
                const CurvedRun c = curved ? curved_coated(n, rf.d, glass, gold, om, adapt) : const_coated(n, rf.d, glass, gold, om);
                const CurvedRun z = curved ? curved_coated(n, rf.d, vac, gold, om, adapt) : const_coated(n, rf.d, vac, gold, om);
                std::printf("  n %2d %s: beschichtet %+.4f %%, neutral %+.4f %%, Schichtwirkung gegen ohne %+.3f %%, gegen neutral %+.3f %%"
                            " (It. %d/%d, Aufbau %.1f s, Loesen %.1f s)\n", n, curved ? "gekruemmt" : "konstant ", 100 * (c.sigma / sig - 1),
                            100 * (z.sigma / sig_bare - 1), 100 * ((c.sigma - sig_bare) / dsig - 1), 100 * ((c.sigma - z.sigma) / dsig - 1), c.it, z.it,
                            c.build, c.solve);
                std::fflush(stdout);
            }
        }
    }
}
