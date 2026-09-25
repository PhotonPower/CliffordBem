#include "cbem/operators/chiral_cauchy_operator.hpp"

namespace cbem {

Mat8 helicity_projector(int sign) {
    Multivector P = Multivector::blade(0, 0.5) + Multivector::blade(7, cplx(0, 0.5 * sign));
    return P.left_matrix();
}

void ChiralCauchyOperator::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    const std::size_t n = x.size();
    std::vector<cplx> xp(n), xm(n), yp, ym;
    for (std::size_t t = 0; t < n / 8; ++t) { cbem::apply(Pp_, &x[8 * t], &xp[8 * t]); cbem::apply(Pm_, &x[8 * t], &xm[8 * t]); }
    Ep_.apply(xp, yp); Em_.apply(xm, ym);
    y.resize(n);
    for (std::size_t i = 0; i < n; ++i) y[i] = yp[i] + ym[i];     // E_k P = P E_k (P zentral)
}

}  // namespace cbem
