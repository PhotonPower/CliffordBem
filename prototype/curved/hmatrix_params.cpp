// H-Matrix gekruemmter Elemente: Aufbauzeit, Speicher, Raenge und Fehler des Produkts in Abhaengigkeit von den
// Partitionsparametern (leaf, eta, sep_factor) und der ACA-Toleranz. Referenz: dieselbe Matrix mit eps = 1e-8.
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -Iinclude prototype/curved/hmatrix_params.cpp build/libcbem.a -o hmatrix_params
//   ./hmatrix_params [n]   (Ikosaederkugel mit 20 n^2 Elementen, Voreinstellung 8)
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "cbem/assembly/curved_entries.hpp"
#include "cbem/geometry/quadratic_mesh.hpp"
#include "cbem/hmatrix/curved_hmatrix.hpp"

using namespace cbem;

int main(int argc, char** argv) {
    const QuadraticMesh q = quadratic_icosphere(argc > 1 ? std::atoi(argv[1]) : 8);
    const cplx k = 0.5 * std::sqrt(cplx(-11, 1.2));
    CurvedKernelEntries E(q, k);
    const std::size_t N = E.size();
    std::mt19937 g(3); std::normal_distribution<double> nd;
    std::vector<cplx> Z(N * kCurvedComps * 8);
    for (auto& z : Z) z = cplx(nd(g), nd(g));
    auto prod = [&](const CurvedHMatrix& H) { std::vector<cplx> Y(N * 8, cplx(0)); H.apply(Z, Y); return Y; };
    HMatrixParams pr; pr.eps = 1e-8;
    const std::vector<cplx> Yr = prod(CurvedHMatrix(E, pr));
    double nr = 0; for (auto y : Yr) nr += std::norm(y);
    std::printf("%zu Elemente, Nahfeld %.2f s\n", q.size(), E.near_seconds());
    struct V { std::size_t leaf; real eta, sep, eps; };
    for (const V& v : std::vector<V>{{32, 1.0, 3.0, 1e-5}, {32, 2.0, 3.0, 1e-5}, {64, 2.0, 3.0, 1e-5}, {32, 1.5, 3.0, 1e-5},
                                     {32, 2.0, 2.0, 1e-5}, {64, 2.0, 2.0, 1e-5}, {48, 2.0, 3.0, 1e-5}}) {
        HMatrixParams hp; hp.leaf = v.leaf; hp.eta = v.eta; hp.sep_factor = v.sep; hp.eps = v.eps;
        auto t0 = std::chrono::steady_clock::now();
        CurvedHMatrix H(E, hp);
        const double tb = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        t0 = std::chrono::steady_clock::now();
        const std::vector<cplx> Y = prod(H);
        const double ta = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        double d = 0; for (std::size_t i = 0; i < Y.size(); ++i) d += std::norm(Y[i] - Yr[i]);
        const HStats& s = H.stats();
        std::printf("leaf %3zu eta %.1f sep %.1f eps %.0e: Aufbau %.2f s, Produkt %.3f s, %.0f MB (dicht %.0f MB), Bloecke %zu/%zu, Rang %.1f/%zu, Fehler %.1e\n",
                    v.leaf, v.eta, v.sep, v.eps, tb, ta, s.bytes() / 1e6, 16.0 * s.entries_dense / 1e6, s.n_dense, s.n_lowrank, s.mean_rank,
                    s.max_rank, std::sqrt(d / nr));
    }
}
