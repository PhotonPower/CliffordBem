#include "cbem/operators/cauchy_operator.hpp"
#include <cmath>

namespace cbem {

CauchyOperator::CauchyOperator(const TriangleMesh& m, const KernelHMatrix& H) : N_(m.size()), H_(H) {
    Ln_.resize(4 * N_); isqA_.resize(N_);
    for (std::size_t j = 0; j < N_; ++j) {
        real sa = std::sqrt(m.area[j]); isqA_[j] = 1.0 / sa;
        Multivector n = Multivector::vector(m.normal[j]);
        Multivector e[4] = {Multivector::blade(0), Multivector::blade(1), Multivector::blade(2), Multivector::blade(4)};
        for (int c = 0; c < 4; ++c) {
            Mat8 L = (e[c] * n).left_matrix();
            for (auto& v : L) v /= sa;
            Ln_[4 * j + c] = L;
        }
    }
}

void CauchyOperator::make_Z(const std::vector<cplx>& x, std::vector<cplx>& Z) const {
    Z.assign(32 * N_, cplx(0));
    for (std::size_t j = 0; j < N_; ++j)
        for (int c = 0; c < 4; ++c) cbem::apply(Ln_[4 * j + c], &x[8 * j], &Z[32 * j + 8 * c]);
}

void CauchyOperator::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    std::vector<cplx> Z; make_Z(x, Z);
    y.assign(8 * N_, cplx(0));
    H_.apply(Z, y);
    for (std::size_t i = 0; i < N_; ++i) for (int q = 0; q < 8; ++q) y[8 * i + q] *= -2.0 * isqA_[i];
}

}  // namespace cbem
