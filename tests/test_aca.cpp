#include <random>
#include "cbem/assembly/kernel_entries.hpp"
#include "cbem/hmatrix/aca.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    TriangleMesh m = make_icosphere(8);
    KernelEntries E(m, cplx(1.5, 0.2));
    std::vector<std::size_t> R, C;
    for (std::size_t t = 0; t < m.size(); ++t) { if (m.centroid[t].z > 0.8) R.push_back(t); if (m.centroid[t].z < -0.8) C.push_back(t); }
    for (int c = 0; c < 4; ++c)
        for (real eps : {1e-3, 1e-6}) {
            RowFn row = [&](std::size_t i, cplx* o) { for (std::size_t j = 0; j < C.size(); ++j) o[j] = E.far(R[i], C[j])[c]; };
            ColFn col = [&](std::size_t j, cplx* o) { for (std::size_t i = 0; i < R.size(); ++i) o[i] = E.far(R[i], C[j])[c]; };
            LowRank lr = aca_partial(row, col, R.size(), C.size(), eps); recompress(lr, eps);
            real num = 0, den = 0;
            for (std::size_t i = 0; i < R.size(); ++i) for (std::size_t j = 0; j < C.size(); ++j) {
                cplx a = E.far(R[i], C[j])[c], s = 0; for (std::size_t k = 0; k < lr.rank(); ++k) s += lr.U(i, k) * lr.V(j, k);
                num += std::norm(a - s); den += std::norm(a);
            }
            real err = std::sqrt(num / den);
            CHECK(err < 20 * eps, "Komponente %d, eps %.0e: Fehler %.1e (Rang %zu)", c, eps, err, lr.rank());
            std::printf("  Komponente %d eps %.0e: Rang %zu, Fehler %.1e\n", c, eps, lr.rank(), err);
        }
    // ACA+ (v0.30): A = u1 v1^T + 0,3 u2 v2^T mit disjunkten Traegern -- die teilpivotisierte ACA sieht den zweiten Summanden nie
    {
        const std::size_t m = 200, n = 150; std::mt19937 rng(1); std::normal_distribution<double> nd(0, 1);
        std::vector<cplx> u1(m, 0.0), v1(n, 0.0), u2(m, 0.0), v2(n, 0.0);
        for (std::size_t i = 0; i < 120; ++i) u1[i] = nd(rng); for (std::size_t i = 120; i < m; ++i) u2[i] = nd(rng);
        for (std::size_t j = 0; j < 90; ++j) v1[j] = nd(rng); for (std::size_t j = 90; j < n; ++j) v2[j] = nd(rng);
        auto A = [&](std::size_t i, std::size_t j) { return u1[i] * v1[j] + 0.3 * u2[i] * v2[j]; };
        RowFn row = [&](std::size_t i, cplx* o) { for (std::size_t j = 0; j < n; ++j) o[j] = A(i, j); };
        ColFn col = [&](std::size_t j, cplx* o) { for (std::size_t i = 0; i < m; ++i) o[i] = A(i, j); };
        real err[2];
        for (int plus = 0; plus < 2; ++plus) {
            LowRank f = aca_select(plus, row, col, m, n, 1e-8);
            double e = 0, s2 = 0;
            for (std::size_t i = 0; i < m; ++i) for (std::size_t j = 0; j < n; ++j) { cplx x = 0; for (std::size_t k = 0; k < f.rank(); ++k) x += f.U(i, k) * f.V(j, k); e += std::norm(x - A(i, j)); s2 += std::norm(A(i, j)); }
            err[plus] = std::sqrt(e / s2);
        }
        std::printf("  disjunkte Traeger: ACA partiell Fehler %.1e, ACA+ %.1e\n", err[0], err[1]);
        CHECK(err[1] < 1e-10, "ACA+ erfasst den Block nicht");
    }
    REPORT();
}
