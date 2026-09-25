#include "cbem/clifford/multivector.hpp"
#include <cmath>

namespace cbem {

Mat8 inverse8(Mat8 A) {
    Mat8 I{}; for (int i = 0; i < 8; ++i) I[i * 8 + i] = 1.0;
    for (int c = 0; c < 8; ++c) {
        int p = c; for (int r = c + 1; r < 8; ++r) if (std::abs(A[r * 8 + c]) > std::abs(A[p * 8 + c])) p = r;
        for (int k = 0; k < 8; ++k) { std::swap(A[c * 8 + k], A[p * 8 + k]); std::swap(I[c * 8 + k], I[p * 8 + k]); }
        const cplx d = A[c * 8 + c];
        for (int k = 0; k < 8; ++k) { A[c * 8 + k] /= d; I[c * 8 + k] /= d; }
        for (int r = 0; r < 8; ++r) if (r != c) {
            const cplx f = A[r * 8 + c];
            for (int k = 0; k < 8; ++k) { A[r * 8 + k] -= f * A[c * 8 + k]; I[r * 8 + k] -= f * I[c * 8 + k]; }
        }
    }
    return I;
}

cplx det8(Mat8 A) {
    cplx det = 1.0;
    for (int c = 0; c < 8; ++c) {
        int p = c; for (int r = c + 1; r < 8; ++r) if (std::abs(A[r * 8 + c]) > std::abs(A[p * 8 + c])) p = r;
        if (A[p * 8 + c] == cplx(0)) return 0.0;
        if (p != c) { for (int k = 0; k < 8; ++k) std::swap(A[c * 8 + k], A[p * 8 + k]); det = -det; }
        det *= A[c * 8 + c];
        for (int r = c + 1; r < 8; ++r) { const cplx f = A[r * 8 + c] / A[c * 8 + c]; for (int k = c; k < 8; ++k) A[r * 8 + k] -= f * A[c * 8 + k]; }
    }
    return det;
}

Multivector mv_inverse(const Multivector& M, bool* ok) {
    const Mat8 L = M.left_matrix();
    const real sc = std::sqrt(mv_norm2(M));
    const cplx d = det8(L);
    const bool good = sc > 0 && std::abs(d) > 1e-24 * std::pow(sc, 8);
    if (ok) *ok = good;
    Multivector R; if (!good) return R;
    const Mat8 Li = inverse8(L);
    for (int b = 0; b < 8; ++b) R.c[b] = Li[b * 8 + 0];          // M^{-1} = L(M)^{-1} * 1
    return R;
}

}  // namespace cbem
