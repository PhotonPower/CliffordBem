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
    REPORT();
}
