#include "cbem/operators/transmission_operator.hpp"
#include <cmath>

namespace cbem {

namespace {
// Gauss-Jordan-Inversion einer 8x8-Matrix
Mat8 inverse8(Mat8 A) {
    Mat8 I{}; for (int i = 0; i < 8; ++i) I[i * 8 + i] = 1.0;
    for (int c = 0; c < 8; ++c) {
        int p = c; for (int r = c + 1; r < 8; ++r) if (std::abs(A[r * 8 + c]) > std::abs(A[p * 8 + c])) p = r;
        for (int k = 0; k < 8; ++k) { std::swap(A[c * 8 + k], A[p * 8 + k]); std::swap(I[c * 8 + k], I[p * 8 + k]); }
        cplx d = A[c * 8 + c];
        for (int k = 0; k < 8; ++k) { A[c * 8 + k] /= d; I[c * 8 + k] /= d; }
        for (int r = 0; r < 8; ++r) if (r != c) {
            cplx f = A[r * 8 + c];
            for (int k = 0; k < 8; ++k) { A[r * 8 + k] -= f * A[c * 8 + k]; I[r * 8 + k] -= f * I[c * 8 + k]; }
        }
    }
    return I;
}
}  // namespace

Mat8 transmission_map(const Vec3& n, const Medium& in, const Medium& out) {
    const cplx se = std::sqrt(in.eps) / std::sqrt(out.eps), sm = std::sqrt(in.mu) / std::sqrt(out.mu);
    const cplx a = sm, b = se;                  // Standardwahl
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
        for (int d = 0; d < 3; ++d) {
            cplx vt = v[d] - vn * n[d], bt = Ib[d] - bn * n[d];
            out8[VEC[d]] = se * vt + vn * n[d] / se;
            out8[BIV[d]] = (sm * bt + bn * n[d] / sm) * BS[d];
        }
        for (int r = 0; r < 8; ++r) J[r * 8 + col] = out8[r];
    }
    return J;
}

TransmissionOperator::TransmissionOperator(const TriangleMesh& m, const CauchyOperator& E1, const CauchyOperator& E2,
                                           const Medium& in, const Medium& out)
    : N_(m.size()), E1_(E1), E2_(E2), J_(m.size()), P_(m.size()) {
    for (std::size_t t = 0; t < N_; ++t) {
        J_[t] = transmission_map(m.normal[t], in, out);
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
