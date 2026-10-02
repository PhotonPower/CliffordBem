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
    for (AcaMode mode : {AcaMode::Joint, AcaMode::Componentwise, AcaMode::Multivector})
        for (real eps : {1e-4, 1e-7}) {
            HMatrixParams p; p.eps = eps; p.mode = mode;
            KernelHMatrix H(E, p); CauchyOperator op(m, H);
            std::vector<cplx> y; op.apply(x, y); auto y0 = dense_apply(E, op, x);
            real num = 0, den = 0; for (std::size_t i = 0; i < y.size(); ++i) { num += std::norm(y[i] - y0[i]); den += std::norm(y0[i]); }
            real err = std::sqrt(num / den);
            std::printf("  %s eps %.0e: %zu dicht / %zu NR, Fehler %.1e\n", mode == AcaMode::Joint ? "joint" : (mode == AcaMode::Multivector ? "mv   " : "comp "), eps, H.stats().n_dense, H.stats().n_lowrank, err);
            CHECK(err < 10 * eps, "Fehler %.1e bei eps %.0e", err, eps);
        }
    // einfache Genauigkeit (v0.56, Voreinstellung) gegen double: Rundungsfehler der Eintraege, halber Speicher
    for (AcaMode mode : {AcaMode::Joint, AcaMode::Componentwise}) {
        HMatrixParams pd; pd.eps = 1e-7; pd.mode = mode; pd.single_precision = false;
        HMatrixParams ps = pd; ps.single_precision = true;
        KernelHMatrix Hd(E, pd), Hs(E, ps); CauchyOperator od(m, Hd), os(m, Hs);
        std::vector<cplx> yd, ys; od.apply(x, yd); os.apply(x, ys);
        real num = 0, den = 0; for (std::size_t i = 0; i < yd.size(); ++i) { num += std::norm(ys[i] - yd[i]); den += std::norm(yd[i]); }
        const real err = std::sqrt(num / den), ratio = Hs.stats().bytes() / Hd.stats().bytes();
        std::printf("  %s float gegen double: %.1e, Speicher %.2f\n", mode == AcaMode::Joint ? "joint" : "comp ", err, ratio);
        CHECK(err < 3e-7 && err > 0, "einfache Genauigkeit: Abweichung %.1e", err);
        CHECK(std::abs(ratio - 0.5) < 1e-12, "einfache Genauigkeit: Speicher %.3f statt 0,5", ratio);
    }
    REPORT();
}
