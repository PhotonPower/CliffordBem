// Nahe getrennte Paare: Fehler gegen eine strenge Referenz und Zeit je Paar fuer Regeln und Kriterien der aeusseren
// Integration und der Korrektur
// (v0.52, docs/results_curved.md). Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -Iinclude prototype/curved/near_rules.cpp build/libcbem.a -o near_rules
//   ./near_rules [n | netz.msh] [Schritt]   (Ikosaeder 20 n^2 oder Gmsh zweiter Ordnung; jedes Schritt-te Element i)
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "cbem/assembly/curved_entries.hpp"
#include "cbem/geometry/quadratic_mesh.hpp"
#include "cbem/geometry/gmsh_io.hpp"
#include <string>

using namespace cbem;

int main(int argc, char** argv) {
    const std::string arg = argc > 1 ? argv[1] : "4";
    const QuadraticMesh q = arg.size() > 3 ? read_gmsh_quadratic(arg)[0].mesh : quadratic_icosphere(std::atoi(arg.c_str()));
    const int step = argc > 2 ? std::atoi(argv[2]) : 7;
    const cplx k(1.3, 0.05);
    EntryParams ep; ep.cache_near = false;
    std::vector<std::pair<std::size_t, std::size_t>> P;
    {
        CurvedKernelEntries E(q, k, ep);
        for (std::size_t i = 0; i < q.size(); i += step)
            for (std::size_t j = 0; j < q.size(); ++j)
                if (E.is_near(i, j) && E.adjacency(i, j) == Adjacency::None) P.push_back({i, j});
    }
    auto run = [&](const CurvedNearParams& np, double& sec) {
        CurvedKernelEntries D(q, k, ep, np);
        std::vector<CurvedBlock> B; B.reserve(P.size());
        auto t0 = std::chrono::steady_clock::now();
        for (auto [i, j] : P) B.push_back(D.lambda_near(i, j));
        sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / P.size();
        return B;
    };
    auto err = [&](const std::vector<CurvedBlock>& A, const std::vector<CurvedBlock>& B, double& rms) {
        double w = 0, s2 = 0;
        for (std::size_t n = 0; n < P.size(); ++n) {
            double na = 0, d = 0;
            for (int p = 0; p < 9; ++p) for (int c = 0; c < kCurvedComps; ++c) { na += std::norm(A[n][p][c]); d += std::norm(B[n][p][c] - A[n][p][c]); }
            w = std::max(w, std::sqrt(d / na)); s2 += d / na;
        }
        rms = std::sqrt(s2 / P.size());
        return w;
    };
    double t, rms;
    // Referenz A: doppelt adaptiv streng (wie test_curved); Referenz B: Subtraktion mit Gauss 8 x 8 und 0,3/0,3
    CurvedNearParams ra; ra.subtract = false; ra.outer_ratio = 0.15; ra.inner_ratio = 0.08; ra.inner_depth = 20;
    CurvedNearParams rb; rb.outer_rule = 8; rb.correction_rule = 8;
    const auto A = run(ra, t);
    const auto Bref = run(rb, t);
    std::printf("%zu Paare; Referenzen gegeneinander: max %.1e\n", P.size(), err(A, Bref, rms));
    struct V { int orule; real oratio; int crule; real cratio; };
    for (const V& v : std::vector<V>{{0, 0.3, 0, 0.3}, {0, 0.5, 0, 0.5}, {4, 1.0, 4, 1.0}, {4, 0.7, 4, 0.7}, {5, 1.5, 5, 1.5}, {5, 2.5, 5, 2.5}}) {
        CurvedNearParams np; np.outer_rule = v.orule; np.subtract_outer_ratio = v.oratio; np.correction_rule = v.crule; np.correction_ratio = v.cratio;
        const auto B = run(np, t);
        const double e = err(Bref, B, rms);
        std::printf("aussen %s%d %.1f, Korrektur %s%d %.1f: max %.1e, rms %.1e, %.0f us je Paar\n", v.orule ? "G" : "D", v.orule ? v.orule : 7,
                    v.oratio, v.crule ? "G" : "D", v.crule ? v.crule : 7, v.cratio, e, rms, 1e6 * t);
    }
}
