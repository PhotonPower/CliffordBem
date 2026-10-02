#include "cbem/operators/curved_operators.hpp"

#include <stdexcept>

#include "cbem/linalg/dense.hpp"
#include "cbem/operators/transmission_operator.hpp"

namespace cbem {

CurvedCauchyOperator::CurvedCauchyOperator(const QuadraticMesh& m, const CurvedHMatrix& H) : NI_(3 * m.size()), H_(H) {
    if (H.size() != NI_) throw std::invalid_argument("CurvedCauchyOperator: H-Matrix mit 3 N Indizes erwartet");
    for (int c = 0; c < kCurvedComps; ++c) Lb_[c] = Multivector::blade(kCurvedBlade[c]).left_matrix();
}

void CurvedCauchyOperator::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    std::vector<cplx> Z(8 * kCurvedComps * NI_, cplx(0));
    for (std::size_t J = 0; J < NI_; ++J)
        for (int c = 0; c < kCurvedComps; ++c) cbem::apply(Lb_[c], &x[8 * J], &Z[8 * kCurvedComps * J + 8 * c]);
    y.assign(8 * NI_, cplx(0));
    H_.apply(Z, y);
    for (auto& v : y) v *= -2.0;
}

CurvedTransmissionOperator::CurvedTransmissionOperator(const QuadraticMesh& m, const std::vector<std::array<real, 9>>& S,
                                                       const BoundaryOperator& E1, const BoundaryOperator& E2,
                                                       const std::vector<Medium>& in, const Medium& out)
    : N_(m.size()), E1_(E1), E2_(E2) {
    if (in.size() != N_ || S.size() != N_) throw std::invalid_argument("CurvedTransmissionOperator: je Element ein Medium und eine Basis");
    CurvedQuadrature Q(m, QuadRule::subdivided(2));
    JG_.assign(N_, std::vector<cplx>(576, cplx(0)));
    P_.assign(N_, std::vector<cplx>(576, cplx(0)));
    for (std::size_t t = 0; t < N_; ++t) {
        std::vector<cplx>& G = JG_[t];
        for (int p = 0; p < Q.q; ++p) {
            const Mat8 J = transmission_map(Q.normals(t)[p], in[t], out);
            real psi[3];
            for (int a = 0; a < 3; ++a) { psi[a] = 0; for (int k = 0; k < 3; ++k) psi[a] += S[t][a * 3 + k] * Q.lam[p][k]; }
            const real w = Q.weights(t)[p];
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
                const real f = w * psi[a] * psi[b];
                for (int r = 0; r < 8; ++r) for (int c = 0; c < 8; ++c) G[(8 * a + r) * 24 + 8 * b + c] += f * J[r * 8 + c];
            }
        }
        // P = 2 (1 + J_G)^{-1}
        Matrix A(24, 24);
        for (int r = 0; r < 24; ++r) for (int c = 0; c < 24; ++c) A(r, c) = G[r * 24 + c] + (r == c ? 1.0 : 0.0);
        std::vector<std::size_t> piv; lu_factor(A, piv);
        for (int c = 0; c < 24; ++c) {
            std::vector<cplx> e(24, cplx(0)); e[c] = 2.0;
            lu_solve(A, piv, e.data());
            for (int r = 0; r < 24; ++r) P_[t][r * 24 + c] = e[r];
        }
    }
}

namespace {
void block_apply(const std::vector<std::vector<cplx>>& B, const std::vector<cplx>& x, std::vector<cplx>& y) {
    y.assign(x.size(), cplx(0));
    for (std::size_t t = 0; t < B.size(); ++t) {
        const cplx* xi = &x[24 * t]; cplx* yi = &y[24 * t]; const std::vector<cplx>& M = B[t];
        for (int r = 0; r < 24; ++r) { cplx s = 0; for (int c = 0; c < 24; ++c) s += M[r * 24 + c] * xi[c]; yi[r] = s; }
    }
}
}  // namespace

void CurvedTransmissionOperator::apply_J(const std::vector<cplx>& x, std::vector<cplx>& y) const { block_apply(JG_, x, y); }

void CurvedTransmissionOperator::precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const { block_apply(P_, x, y); }

void CurvedTransmissionOperator::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    std::vector<cplx> e2, jx, e1;
    E2_.apply(x, e2);
    apply_J(x, jx);
    E1_.apply(jx, e1);
    y.resize(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) y[i] = 0.5 * (x[i] + e2[i]) + 0.5 * (jx[i] - e1[i]);
}

}  // namespace cbem
