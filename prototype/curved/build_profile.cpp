// Aufbau gekruemmter Elemente: Aufteilung der Zeit (Nahfeld nach Sauter-Schwab / getrennte Nahpaare, H-Matrizen, Rest)
// fuer eine Kugel aus Gold (zwei Wellenzahlen: aussen und innen).
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -Iinclude prototype/curved/build_profile.cpp build/libcbem.a -o build_profile
//   ./build_profile [n]   (Ikosaederkugel mit 20 n^2 Elementen, Voreinstellung 8)
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "cbem/hmatrix/curved_hmatrix.hpp"
#include "cbem/problems/curved_problem.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

int main(int argc, char** argv) {
    const int n = argc > 1 ? std::atoi(argv[1]) : 8;
    const real om = 0.5; const Medium gold{cplx(-11, 1.2)};
    const QuadraticMesh q = quadratic_icosphere(n);
    auto t0 = clk::now();
    CurvedScatteringProblem P({q}, {gold}, om);
    const double total = since(t0);
    std::printf("%zu Elemente: Aufbau CurvedScatteringProblem %.2f s (Nahfeld-Cache %.2f s)\n", q.size(), total, P.near_seconds());
    double hsum = 0, ssum = 0, nsum = 0;
    for (cplx k : {Medium{}.k(om), gold.k(om)}) {
        EntryParams ep; ep.cache_near = false;
        CurvedKernelEntries E(q, k, ep);
        // Nahpaare nach Art
        std::size_t nss[4] = {0, 0, 0, 0}, nsep = 0;
        double tss[4] = {0, 0, 0, 0}, tsep = 0;
        for (std::size_t i = 0; i < q.size(); ++i)
            for (std::size_t j = 0; j < q.size(); ++j) {
                if (!E.is_near(i, j)) continue;
                const Adjacency a = E.adjacency(i, j);
                auto t1 = clk::now();
                if (a == Adjacency::None) { volatile auto K = E.lambda_near(i, j)[0][0]; (void)K; tsep += since(t1); ++nsep; }
                else { volatile auto K = E.lambda_sauter_schwab(i, j, a)[0][0]; (void)K; tss[static_cast<int>(a)] += since(t1); ++nss[static_cast<int>(a)]; }
            }
        // H-Matrix mit Nahfeld-Cache (wie im Problem)
        CurvedKernelEntries Ec(q, k);
        t0 = clk::now();
        CurvedHMatrix H(Ec);
        const double th = since(t0);
        std::printf("k = %.3f%+.3fi: Sauter-Schwab %.2f s (Ecke %zu: %.2f s, Kante %zu: %.2f s, Selbst %zu: %.2f s), getrennt %zu: %.2f s, "
                    "H-Matrix %.2f s (%.0f MB, Rang %.1f)\n", k.real(), k.imag(), tss[1] + tss[2] + tss[3], nss[1], tss[1], nss[2], tss[2], nss[3], tss[3],
                    nsep, tsep, th, H.stats().bytes() / 1e6, H.stats().mean_rank);
        hsum += th; ssum += tss[1] + tss[2] + tss[3]; nsum += tsep;
    }
    std::printf("Summe: Sauter-Schwab %.2f s, getrennt %.2f s, H-Matrizen %.2f s, Rest (Transmission, psi, Quadratur, ...) %.2f s\n", ssum, nsum, hsum,
                total - P.near_seconds() - hsum);
}
