#include "cbem/operators/linear_cauchy_operator.hpp"

#include <stdexcept>

namespace cbem {

LinearCauchyOperator::LinearCauchyOperator(const TriangleMesh& m, const KernelHMatrix& H) : NI_(3 * m.size()), H_(H) {
    if (H.size() != NI_) throw std::invalid_argument("LinearCauchyOperator: H-Matrix mit 3 N Indizes erwartet (LinearKernelEntries)");
    Ln_.resize(4 * NI_);
    const Multivector e[4] = {Multivector::blade(0), Multivector::blade(1), Multivector::blade(2), Multivector::blade(4)};
    for (std::size_t t = 0; t < m.size(); ++t) {
        const Multivector n = Multivector::vector(m.normal[t]);
        for (int c = 0; c < 4; ++c) {
            const Mat8 L = (e[c] * n).left_matrix();
            for (int a = 0; a < 3; ++a) Ln_[4 * (3 * t + a) + c] = L;
        }
    }
}

void LinearCauchyOperator::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    std::vector<cplx> Z(32 * NI_, cplx(0));
    for (std::size_t J = 0; J < NI_; ++J)
        for (int c = 0; c < 4; ++c) cbem::apply(Ln_[4 * J + c], &x[8 * J], &Z[32 * J + 8 * c]);
    y.assign(8 * NI_, cplx(0));
    H_.apply(Z, y);
    for (auto& v : y) v *= -2.0;
}

}  // namespace cbem
