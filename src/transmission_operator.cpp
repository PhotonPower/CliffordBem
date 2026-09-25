#include "cbem/operators/transmission_operator.hpp"
#include <cmath>

namespace cbem {


Mat8 transmission_map(const Vec3& n, const Medium& in, const Medium& out) {
    const cplx se = std::sqrt(in.eps) / std::sqrt(out.eps), sm = std::sqrt(in.mu) / std::sqrt(out.mu);
    const cplx a = sm, b = se;                  // Standardwahl
    // Normalanteile: Nn = D1 C1^{-1} C2 D2^{-1}
    const cplx I1(0, 1);
    const cplx c1[2][2] = {{in.eps, I1 * in.chi}, {-I1 * in.chi, in.mu}}, c2[2][2] = {{out.eps, I1 * out.chi}, {-I1 * out.chi, out.mu}};
    const cplx det = c1[0][0] * c1[1][1] - c1[0][1] * c1[1][0];
    const cplx ci[2][2] = {{c1[1][1] / det, -c1[0][1] / det}, {-c1[1][0] / det, c1[0][0] / det}};
    const cplx d1[2] = {std::sqrt(in.eps), std::sqrt(in.mu)}, d2[2] = {std::sqrt(out.eps), std::sqrt(out.mu)};
    cplx Nn[2][2];
    for (int r = 0; r < 2; ++r) for (int c = 0; c < 2; ++c) {
        cplx s = 0; for (int q = 0; q < 2; ++q) s += ci[r][q] * c2[q][c];
        Nn[r][c] = d1[r] * s / d2[c];
    }
    Mat8 J{};
    const int VEC[3] = {1, 2, 4};               // e1, e2, e3
    const int BIV[3] = {6, 5, 3};               // I e1 = e23, I e2 = -e13, I e3 = e12
    const real BS[3] = {1, -1, 1};
    for (int col = 0; col < 8; ++col) {
        cplx h[8] = {}; h[col] = 1.0;
        CVec3 v, Ib;
        for (int d = 0; d < 3; ++d) { v[d] = h[VEC[d]]; Ib[d] = h[BIV[d]] * BS[d]; }
        cplx vn = v[0] * n.x + v[1] * n.y + v[2] * n.z, bn = Ib[0] * n.x + Ib[1] * n.y + Ib[2] * n.z;
        cplx out8[8] = {};
        out8[0] = a * h[0]; out8[7] = b * h[7];
        const cplx vn2 = Nn[0][0] * vn + Nn[0][1] * bn, bn2 = Nn[1][0] * vn + Nn[1][1] * bn;
        for (int d = 0; d < 3; ++d) {
            cplx vt = v[d] - vn * n[d], bt = Ib[d] - bn * n[d];
            out8[VEC[d]] = se * vt + vn2 * n[d];
            out8[BIV[d]] = (sm * bt + bn2 * n[d]) * BS[d];
        }
        for (int r = 0; r < 8; ++r) J[r * 8 + col] = out8[r];
    }
    return J;
}

TransmissionOperator::TransmissionOperator(const TriangleMesh& m, const BoundaryOperator& E1, const BoundaryOperator& E2,
                                           const Medium& in, const Medium& out)
    : N_(m.size()), E1_(E1), E2_(E2) { setup(m, std::vector<Medium>(m.size(), in), out); }

TransmissionOperator::TransmissionOperator(const TriangleMesh& m, const BoundaryOperator& E1, const BoundaryOperator& E2,
                                           const std::vector<Medium>& in, const Medium& out)
    : N_(m.size()), E1_(E1), E2_(E2) { setup(m, in, out); }

void TransmissionOperator::setup(const TriangleMesh& m, const std::vector<Medium>& in, const Medium& out) {
    J_.resize(N_); P_.resize(N_);
    for (std::size_t t = 0; t < N_; ++t) {
        J_[t] = transmission_map(m.normal[t], in[t], out);
        Mat8 A = J_[t]; for (int i = 0; i < 8; ++i) A[i * 8 + i] += 1.0;
        P_[t] = inverse8(A); for (auto& v : P_[t]) v *= 2.0;
    }
}

void TransmissionOperator::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    std::vector<cplx> Jx(8 * N_), E1Jx, E2x;
    for (std::size_t t = 0; t < N_; ++t) cbem::apply(J_[t], &x[8 * t], &Jx[8 * t]);
    E2_.apply(x, E2x); E1_.apply(Jx, E1Jx);
    y.resize(8 * N_);
    for (std::size_t i = 0; i < 8 * N_; ++i) y[i] = 0.5 * (x[i] + E2x[i]) + 0.5 * (Jx[i] - E1Jx[i]);
}

void TransmissionOperator::precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    y.resize(8 * N_);
    for (std::size_t t = 0; t < N_; ++t) cbem::apply(P_[t], &x[8 * t], &y[8 * t]);
}

}  // namespace cbem
