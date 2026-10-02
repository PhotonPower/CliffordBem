// Sauter-Schwab auf gekruemmten Elementen: Fehler je Nachbarschaftstyp (Selbstterm, Kante, Ecke) gegen eine hohe Ordnung
// und Zeit je Paar, abhaengig von EntryParams::ss_order (Gauss-Punkte je Richtung im 4D-Wuerfel).
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -Iinclude prototype/curved/ss_orders.cpp build/libcbem.a -o ss_orders
//   ./ss_orders [n | netz.msh] [Referenzordnung] [xi,eta1,eta2,eta3 ...]
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "cbem/assembly/curved_entries.hpp"
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/geometry/quadratic_mesh.hpp"

using namespace cbem;

int main(int argc, char** argv) {
    const std::string arg = argc > 1 ? argv[1] : "4";
    const QuadraticMesh q = arg.size() > 3 ? read_gmsh_quadratic(arg)[0].mesh : quadratic_icosphere(std::atoi(arg.c_str()));
    const int ref_order = argc > 2 ? std::atoi(argv[2]) : 12;
    const cplx k = std::getenv("CBEM_K_RE") ? cplx(std::atof(std::getenv("CBEM_K_RE")), std::atof(std::getenv("CBEM_K_IM"))) : cplx(0.09, 1.66);   // Voreinstellung: Gold innen bei omega = 0,5
    EntryParams ep; ep.cache_near = false;
    std::vector<std::pair<std::size_t, std::size_t>> P[4];
    {
        CurvedKernelEntries E(q, k, ep);
        for (std::size_t i = 0; i < q.size(); i += 3)
            for (std::size_t j = 0; j < q.size(); ++j) {
                const Adjacency a = E.adjacency(i, j);
                if (a != Adjacency::None) P[static_cast<int>(a)].push_back({i, j});
            }
    }
    const char* name[4] = {"", "Ecke", "Kante", "Selbstterm"};
    std::vector<CurvedBlock> ref[4];
    {
        EntryParams er = ep; er.ss_order = ref_order;
        CurvedKernelEntries R(q, k, er);
        for (int t = 1; t < 4; ++t) for (auto [i, j] : P[t]) ref[t].push_back(R.lambda_exact(i, j));
    }
    std::printf("%zu Elemente; Paare: Ecke %zu, Kante %zu, Selbstterm %zu; Referenzordnung %d\n", q.size(), P[1].size(), P[2].size(),
                P[3].size(), ref_order);
    // Ordnungen je Richtung (xi, eta1, eta2, eta3), fuer alle drei Typen gleich; Punkte je Teilgebiet = Produkt
    std::vector<std::array<int, 4>> sets = {{5, 5, 5, 5}, {6, 6, 6, 6}, {8, 5, 5, 5}, {5, 8, 5, 5}, {5, 5, 8, 5}, {5, 5, 5, 8}};
    if (argc > 3) {   // weitere Saetze als Argumente "a,b,c,d"
        sets.clear();
        for (int a = 3; a < argc; ++a) { std::array<int, 4> o; std::sscanf(argv[a], "%d,%d,%d,%d", &o[0], &o[1], &o[2], &o[3]); sets.push_back(o); }
    }
    for (const auto& o : sets) {
        CurvedNearParams np; np.ss_orders = {o, o, o};
        CurvedKernelEntries D(q, k, ep, np);
        std::printf("Ordnungen %d,%d,%d,%d (%3d):", o[0], o[1], o[2], o[3], o[0] * o[1] * o[2] * o[3]);
        for (int t = 1; t < 4; ++t) {
            double w = 0;
            auto t0 = std::chrono::steady_clock::now();
            std::vector<CurvedBlock> B; for (auto [i, j] : P[t]) B.push_back(D.lambda_exact(i, j));
            const double us = 1e6 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / P[t].size();
            std::vector<double> es;
            for (std::size_t n = 0; n < P[t].size(); ++n) {
                double na = 0, d = 0;
                for (int p = 0; p < 9; ++p) for (int c = 0; c < kCurvedComps; ++c) { na += std::norm(ref[t][n][p][c]); d += std::norm(B[n][p][c] - ref[t][n][p][c]); }
                es.push_back(std::sqrt(d / na)); w = std::max(w, es.back());
            }
            if (std::getenv("CBEM_WORST")) {   // schlechteste Paare mit Formmassen
                std::vector<std::size_t> id(es.size()); for (std::size_t n = 0; n < id.size(); ++n) id[n] = n;
                std::sort(id.begin(), id.end(), [&](std::size_t a, std::size_t b) { return es[a] > es[b]; });
                auto shape = [&](std::size_t e, real& minang, real& ratio) {
                    const auto v = q.flat.vertices(e); minang = 10;
                    real lmax = 0;
                    for (int c = 0; c < 3; ++c) {
                        const Vec3 a1 = v[(c + 1) % 3] - v[c], a2 = v[(c + 2) % 3] - v[c];
                        minang = std::min(minang, std::acos(dot(a1, a2) / (norm(a1) * norm(a2))));
                        lmax = std::max(lmax, norm(a1));
                    }
                    ratio = lmax * lmax / (2 * q.flat.area[e]);   // lmax / Hoehe auf lmax; gleichseitig 1,15
                };
                for (std::size_t r = 0; r < 6 && r < id.size(); ++r) {
                    const auto [i, j] = P[t][id[r]];
                    real ai, ri, aj, rj; shape(i, ai, ri); shape(j, aj, rj);
                    std::printf("\n    %s %.1e: Winkel min %.0f/%.0f Grad, lmax/Hoehe %.2f/%.2f, h %.3f/%.3f, Woelbung/h %.3f/%.3f", name[t], es[id[r]],
                                ai * 180 / pi, aj * 180 / pi, ri, rj, q.flat.hmax[i], q.flat.hmax[j], q.bulge[i] / q.flat.hmax[i], q.bulge[j] / q.flat.hmax[j]);
                }
                std::printf("\n");
            }
            std::sort(es.begin(), es.end());
            if (std::getenv("CBEM_QUANTILES"))
                std::printf("  %s Median %.1e, 90 %% %.1e, 99 %% %.1e, max %.1e (%.0f us)", name[t], es[es.size() / 2], es[es.size() * 9 / 10],
                            es[es.size() * 99 / 100], w, us);
            else std::printf("  %s %.1e (%.0f us)", name[t], w, us);
        }
        std::printf("\n");
    }
}
