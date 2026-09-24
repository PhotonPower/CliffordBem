#include "cbem/operators/cauchy_operator.hpp"
#include <random>
#include "check.hpp"
using namespace cbem;
// Referenz: E x mit exakten Eintraegen fuer alle Paare
static std::vector<cplx> dense_apply(const KernelEntries& E, const CauchyOperator& op, const std::vector<cplx>& x) {
    const TriangleMesh& m = E.mesh(); const std::size_t N = m.size();
    std::vector<cplx> Z; op.make_Z(x, Z); std::vector<cplx> y(8 * N, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) { KernelComp K = E.exact(i, j); for (int c = 0; c < 4; ++c) for (int q = 0; q < 8; ++q) y[8 * i + q] += K[c] * Z[32 * j + 8 * c + q]; }
        for (int q = 0; q < 8; ++q) y[8 * i + q] *= -2.0 / std::sqrt(m.area[i]);
    }
    return y;
}
int main() {
    std::mt19937 g(3); std::normal_distribution<real> nd;
    TriangleMesh m = make_icosphere(6);
    KernelEntries E(m, 1.5);
    std::vector<cplx> x(8 * m.size()); for (auto& v : x) v = cplx(nd(g), nd(g));
    for (AcaMode mode : {AcaMode::Joint, AcaMode::Componentwise})
        for (real eps : {1e-4, 1e-7}) {
            HMatrixParams p; p.eps = eps; p.mode = mode;
            KernelHMatrix H(E, p); CauchyOperator op(m, H);
            std::vector<cplx> y; op.apply(x, y); auto y0 = dense_apply(E, op, x);
            real num = 0, den = 0; for (std::size_t i = 0; i < y.size(); ++i) { num += std::norm(y[i] - y0[i]); den += std::norm(y0[i]); }
            real err = std::sqrt(num / den);
            std::printf("  %s eps %.0e: %zu dicht / %zu NR, Fehler %.1e\n", mode == AcaMode::Joint ? "joint" : "comp ", eps, H.stats().n_dense, H.stats().n_lowrank, err);
            CHECK(err < 10 * eps, "Fehler %.1e bei eps %.0e", err, eps);
        }
    REPORT();
}
