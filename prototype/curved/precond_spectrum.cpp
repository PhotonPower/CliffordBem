// Spektrum des vorkonditionierten Operators gekruemmter Elemente (Arnoldi): Ritz-Werte von T M (M = 2 (1 + J_G)^{-1} je
// Element, oder ohne) und der GMRES-Residuenverlauf aus der Hessenberg-Matrix. Zeigt, ob wenige Ausreisser (Resonanzen) die
// Iterationen bestimmen (dann helfen Grobgitter-Korrektur oder Deflation) oder ein breites Spektrum (dann lokale Vorkond.).
// Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -Iinclude prototype/curved/precond_spectrum.cpp build/libcbem.a -o precond_spectrum
//   ./precond_spectrum [n] [Arnoldi-Schritte]
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "cbem/linalg/dense.hpp"
#include "cbem/problems/curved_problem.hpp"

using namespace cbem;

static cplx dotc(const std::vector<cplx>& a, const std::vector<cplx>& b) { cplx s = 0; for (std::size_t i = 0; i < a.size(); ++i) s += std::conj(a[i]) * b[i]; return s; }

int main(int argc, char** argv) {
    const int n = argc > 1 ? std::atoi(argv[1]) : 4, m = argc > 2 ? std::atoi(argv[2]) : 80;
    const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    struct Case { const char* name; Medium in, out; real om; };
    for (const Case& c : {Case{"Glas", Medium{2.25}, Medium{}, 1.0}, Case{"Gold", Medium{cplx(-11, 1.2)}, Medium{}, 0.5}}) {
        CurvedScatteringProblem P({quadratic_icosphere(n)}, {c.in}, c.om, c.out);
        const auto& T = P.T();
        const auto b = project_plane_wave_curved(P.mesh(), c.out.k(c.om), c.out.eps, d, px);
        for (int pre : {1, 0}) {
            // Arnoldi auf B = T M (Rechtsvorkonditionierung wie gmres), Start b / |b|
            std::vector<std::vector<cplx>> V(1, b);
            const real nb = std::sqrt(std::real(dotc(b, b)));
            for (auto& v : V[0]) v /= nb;
            Matrix H(m + 1, m);
            int k = 0;
            std::vector<cplx> w, t;
            for (; k < m; ++k) {
                if (pre) { T.precondition(V[k], t); T.apply(t, w); } else T.apply(V[k], w);
                for (int j = 0; j <= k; ++j) { const cplx h = dotc(V[j], w); H(j, k) = h; for (std::size_t i = 0; i < w.size(); ++i) w[i] -= h * V[j][i]; }
                const real hn = std::sqrt(std::real(dotc(w, w)));
                H(k + 1, k) = hn;
                if (hn < 1e-14) { ++k; break; }
                for (auto& v : w) v /= hn;
                V.push_back(w);
            }
            // GMRES-Residuen: min |beta e1 - H_j y| fuer j = 1..k (Givens)
            std::vector<real> res;
            {
                Matrix R(k + 1, k); for (int j = 0; j < k; ++j) for (int i = 0; i <= k; ++i) R(i, j) = H(i, j);
                std::vector<cplx> g(k + 1, 0); g[0] = 1.0;
                for (int j = 0; j < k; ++j) {
                    const cplx a = R(j, j), bb = R(j + 1, j);
                    const real r = std::sqrt(std::norm(a) + std::norm(bb));
                    const cplx cs = a / r, sn = bb / r;
                    for (int jj = j; jj < k; ++jj) {
                        const cplx x0 = R(j, jj), x1 = R(j + 1, jj);
                        R(j, jj) = std::conj(cs) * x0 + std::conj(sn) * x1; R(j + 1, jj) = -sn * x0 + cs * x1;
                    }
                    const cplx g0 = g[j], g1 = g[j + 1];
                    g[j] = std::conj(cs) * g0 + std::conj(sn) * g1; g[j + 1] = -sn * g0 + cs * g1;
                    res.push_back(std::abs(g[j + 1]));
                }
            }
            // Ritz-Werte
            Matrix Hk(k, k); for (int j = 0; j < k; ++j) for (int i = 0; i < k; ++i) Hk(i, j) = H(i, j);
            std::vector<cplx> lam; Matrix Vv; eig_complex(Hk, lam, Vv);
            std::sort(lam.begin(), lam.end(), [](cplx a, cplx b) { return std::abs(a - 1.0) > std::abs(b - 1.0); });
            int it6 = -1, it10 = -1;
            for (int j = 0; j < static_cast<int>(res.size()); ++j) { if (it6 < 0 && res[j] < 1e-6) it6 = j + 1; if (it10 < 0 && res[j] < 1e-10) it10 = j + 1; }
            std::printf("%-5s %zu Elemente %s: GMRES 1e-6 nach %d, 1e-10 nach %d It.; Residuum nach 5/10/20/40 It.: %.1e %.1e %.1e %.1e\n", c.name,
                        P.mesh().size(), pre ? "mit Vorkond. " : "ohne Vorkond.", it6, it10, res[std::min<int>(4, res.size() - 1)],
                        res[std::min<int>(9, res.size() - 1)], res[std::min<int>(19, res.size() - 1)], res[std::min<int>(39, res.size() - 1)]);
            std::printf("    Ritz-Werte, am weitesten von 1:");
            for (int j = 0; j < std::min<int>(8, lam.size()); ++j) std::printf(" %.3f%+.3fi", lam[j].real(), lam[j].imag());
            real rmin = 1e300, rmax = 0; for (auto l : lam) { rmin = std::min(rmin, std::abs(l)); rmax = std::max(rmax, std::abs(l)); }
            std::printf("\n    |lambda| von %.3f bis %.3f\n", rmin, rmax);
        }
    }
}
