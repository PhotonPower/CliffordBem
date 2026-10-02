// ACA auf Quadraturpunkten statt auf Elementpaaren: Rang, Kernauswertungen und Zeit fuer zufaellige zulaessige Clusterpaare
// (je nc Elemente um zwei Zentren; Abstand wie bei eta = 2, sep_factor = 3). Elementebene: ACA+ ueber die 7 gestapelten
// Komponenten mit block_far je Elementpaar (wie CurvedHMatrix); Punktebene: ACA+ ueber die 4 Kernkomponenten (s, v z) an
// den 7 Dunavant-Punkten je Element, ohne Projektion auf die Basis (untere Schranke des Aufwands).
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -Iinclude prototype/curved/point_aca.cpp build/libcbem.a -o point_aca
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <numeric>
#include <random>
#include <vector>

#include "cbem/assembly/curved_entries.hpp"
#include "cbem/geometry/quadratic_mesh.hpp"
#include "cbem/geometry/quadrature.hpp"
#include "cbem/hmatrix/aca.hpp"
#include "cbem/kernel/dirac_kernel.hpp"

using namespace cbem;
using clk = std::chrono::steady_clock;

int main(int argc, char** argv) {
    const QuadraticMesh q = quadratic_icosphere(8);
    const cplx k = 0.5 * std::sqrt(cplx(-11, 1.2));
    const std::size_t nc = argc > 1 ? std::atoi(argv[1]) : 22;
    const real eps = 1e-6;
    CurvedKernelEntries E(q, k);
    CurvedQuadrature Q(q, QuadRule::dunavant7());
    const int nq = Q.q;
    std::mt19937 g(5); std::uniform_int_distribution<std::size_t> u(0, q.size() - 1);
    auto cluster = [&](std::size_t c0) {
        std::vector<std::size_t> id(q.size()); std::iota(id.begin(), id.end(), 0);
        std::partial_sort(id.begin(), id.begin() + nc, id.end(), [&](std::size_t a, std::size_t b) {
            return norm(q.flat.centroid[a] - q.flat.centroid[c0]) < norm(q.flat.centroid[b] - q.flat.centroid[c0]); });
        id.resize(nc); return id;
    };
    auto diam = [&](const std::vector<std::size_t>& c) { real d = 0; for (auto a : c) for (auto b : c) d = std::max(d, norm(q.flat.centroid[a] - q.flat.centroid[b])); return d + 2 * q.flat.hmax[c[0]]; };
    auto dist = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) { real d = 1e300; for (auto x : a) for (auto y : b) d = std::min(d, norm(q.flat.centroid[x] - q.flat.centroid[y])); return d - 2 * q.flat.hmax[a[0]]; };
    double te = 0, tp = 0; double re = 0, rp = 0; long ee = 0, ep = 0; int cnt = 0;
    while (cnt < 30) {
        const auto A = cluster(u(g)), B = cluster(u(g));
        const real d = dist(A, B), da = diam(A), db = diam(B);
        if (!(std::min(da, db) <= 2.0 * d && d > 3 * q.flat.hmax[A[0]] && std::max(da, db) <= 4.0 * d)) continue;   // wie eta = 2, nicht zu weit
        ++cnt;
        // Elementebene
        {
            const std::size_t m = 3 * nc, n = 3 * nc; long evals = 0;
            std::vector<std::vector<cplx>> rowc(nc), colc(nc);
            auto t0 = clk::now();
            RowFn row = [&](std::size_t i, cplx* out) {
                auto& rc = rowc[i / 3];
                if (rc.empty()) {
                    rc.assign(3 * 7 * n, 0);
                    for (std::size_t j = 0; j < nc; ++j) { const CurvedBlock K = E.block_far(A[i / 3], B[j]); ++evals;
                        for (int p = 0; p < 3; ++p) for (int qq = 0; qq < 3; ++qq) for (int c = 0; c < 7; ++c) rc[p * 7 * n + c * n + 3 * j + qq] = K[p * 3 + qq][c]; }
                }
                std::copy(rc.begin() + (i % 3) * 7 * n, rc.begin() + (i % 3 + 1) * 7 * n, out);
            };
            ColFn col = [&](std::size_t J, cplx* out) {
                const std::size_t c = J / n, j = J % n; auto& cc = colc[j / 3];
                if (cc.empty()) {
                    cc.assign(3 * 7 * m, 0);
                    for (std::size_t i = 0; i < nc; ++i) { const CurvedBlock K = E.block_far(A[i], B[j / 3]); ++evals;
                        for (int p = 0; p < 3; ++p) for (int qq = 0; qq < 3; ++qq) for (int c7 = 0; c7 < 7; ++c7) cc[(qq * 7 + c7) * m + 3 * i + p] = K[p * 3 + qq][c7]; }
                }
                std::copy(cc.begin() + ((j % 3) * 7 + c) * m, cc.begin() + ((j % 3) * 7 + c + 1) * m, out);
            };
            LowRank f = aca_plus(row, col, m, 7 * n, eps); recompress(f, eps);
            te += std::chrono::duration<double>(clk::now() - t0).count(); re += f.rank(); ee += evals * nq * nq;
        }
        // Punktebene (4 Komponenten: s und v z), Zeilen: Punkte x, Spalten: (Komponente, Punkt y)
        {
            const std::size_t m = nc * nq, n = nc * nq; long evals = 0;
            auto pt = [&](const std::vector<std::size_t>& C, std::size_t p) { return Q.points(C[p / nq])[p % nq]; };
            auto t0 = clk::now();
            RowFn row = [&](std::size_t i, cplx* out) {
                const Vec3 x = pt(A, i);
                for (std::size_t j = 0; j < n; ++j) { const Vec3 z = x - pt(B, j); const KernelValue kv = dirac_kernel_fast(z, k); ++evals;
                    out[j] = kv.s; out[n + j] = kv.vcoef * z.x; out[2 * n + j] = kv.vcoef * z.y; out[3 * n + j] = kv.vcoef * z.z; }
            };
            std::vector<std::vector<cplx>> colc(n);
            ColFn col = [&](std::size_t J, cplx* out) {
                const std::size_t c = J / n, j = J % n; auto& cc = colc[j];
                if (cc.empty()) {
                    cc.assign(4 * m, 0); const Vec3 y = pt(B, j);
                    for (std::size_t i = 0; i < m; ++i) { const Vec3 z = pt(A, i) - y; const KernelValue kv = dirac_kernel_fast(z, k); ++evals;
                        cc[i] = kv.s; cc[m + i] = kv.vcoef * z.x; cc[2 * m + i] = kv.vcoef * z.y; cc[3 * m + i] = kv.vcoef * z.z; }
                }
                std::copy(cc.begin() + c * m, cc.begin() + (c + 1) * m, out);
            };
            LowRank f = aca_plus(row, col, m, 4 * n, eps); recompress(f, eps);
            tp += std::chrono::duration<double>(clk::now() - t0).count(); rp += f.rank(); ep += evals;
        }
    }
    std::printf("%d Clusterpaare zu %zu Elementen, eps %.0e\n", cnt, nc, eps);
    std::printf("Elementebene: Rang %.1f, Kernauswertungen %.0f, %.2f ms je Block\n", re / cnt, double(ee) / cnt, 1e3 * te / cnt);
    std::printf("Punktebene:   Rang %.1f, Kernauswertungen %.0f, %.2f ms je Block (ohne Projektion auf die Basis)\n", rp / cnt, double(ep) / cnt, 1e3 * tp / cnt);
}
