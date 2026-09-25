#include "cbem/linalg/dense.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace cbem {

namespace {
cplx dotc(const cplx* x, const cplx* y, std::size_t n) { cplx s = 0; for (std::size_t i = 0; i < n; ++i) s += std::conj(x[i]) * y[i]; return s; }
real nrm(const cplx* x, std::size_t n) { real s = 0; for (std::size_t i = 0; i < n; ++i) s += std::norm(x[i]); return std::sqrt(s); }
}  // namespace

real frobenius(const Matrix& A) { return nrm(A.a.data(), A.a.size()); }

void qr_cgs2(const Matrix& A, Matrix& Q, Matrix& R) {
    const std::size_t m = A.rows, n = A.cols;
    Q = A; R = Matrix(n, n);
    for (std::size_t j = 0; j < n; ++j) {
        cplx* qj = Q.col(j);
        for (int pass = 0; pass < 2; ++pass)
            for (std::size_t i = 0; i < j; ++i) {
                cplx h = dotc(Q.col(i), qj, m); R(i, j) += h;
                const cplx* qi = Q.col(i); for (std::size_t r = 0; r < m; ++r) qj[r] -= h * qi[r];
            }
        real nj = nrm(qj, m); R(j, j) = nj;
        if (nj > 0) for (std::size_t r = 0; r < m; ++r) qj[r] /= nj;
    }
}

void svd_jacobi(const Matrix& A, Matrix& U, std::vector<real>& s, Matrix& V, real tol, int max_sweeps) {
    const std::size_t m = A.rows, n = A.cols;
    U = A; V = Matrix(n, n); for (std::size_t i = 0; i < n; ++i) V(i, i) = 1.0;
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        real off = 0;
        for (std::size_t p = 0; p + 1 < n; ++p)
            for (std::size_t q = p + 1; q < n; ++q) {
                cplx* up = U.col(p); cplx* uq = U.col(q);
                real alpha = std::real(dotc(up, up, m)), beta = std::real(dotc(uq, uq, m));
                cplx gamma = dotc(up, uq, m); real ag = std::abs(gamma);
                if (ag <= tol * std::sqrt(alpha * beta) || ag == 0) continue;
                off = std::max(off, ag / std::sqrt(alpha * beta));
                // Rotation im Raum der Spalten (p,q): Phase von gamma abspalten, dann reelle Jacobi-Rotation
                cplx ph = gamma / ag;
                real zeta = (beta - alpha) / (2 * ag);
                real t = (zeta >= 0 ? 1.0 : -1.0) / (std::abs(zeta) + std::sqrt(1 + zeta * zeta));
                real c = 1 / std::sqrt(1 + t * t), sn = c * t;
                for (std::size_t r = 0; r < m; ++r) {
                    cplx a = up[r], b = uq[r] * std::conj(ph);
                    up[r] = c * a - sn * b; uq[r] = (sn * a + c * b) * ph;
                }
                cplx* vp = V.col(p); cplx* vq = V.col(q);
                for (std::size_t r = 0; r < n; ++r) {
                    cplx a = vp[r], b = vq[r] * std::conj(ph);
                    vp[r] = c * a - sn * b; vq[r] = (sn * a + c * b) * ph;
                }
            }
        if (off < tol) break;
    }
    s.assign(n, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
        s[j] = nrm(U.col(j), m);
        if (s[j] > 0) for (std::size_t r = 0; r < m; ++r) U(r, j) /= s[j];
    }
    std::vector<std::size_t> ord(n); std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return s[a] > s[b]; });
    Matrix U2(m, n), V2(n, n); std::vector<real> s2(n);
    for (std::size_t k = 0; k < n; ++k) {
        s2[k] = s[ord[k]];
        std::copy(U.col(ord[k]), U.col(ord[k]) + m, U2.col(k));
        std::copy(V.col(ord[k]), V.col(ord[k]) + n, V2.col(k));
    }
    U = std::move(U2); V = std::move(V2); s = std::move(s2);
}

void lu_factor(Matrix& A, std::vector<std::size_t>& piv) {
    const std::size_t n = A.rows; piv.resize(n);
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k; real best = std::abs(A(k, k));
        for (std::size_t i = k + 1; i < n; ++i) if (std::abs(A(i, k)) > best) { best = std::abs(A(i, k)); p = i; }
        piv[k] = p;
        if (p != k) for (std::size_t j = 0; j < n; ++j) std::swap(A(k, j), A(p, j));
        const cplx d = A(k, k);
        for (std::size_t i = k + 1; i < n; ++i) A(i, k) /= d;
        for (std::size_t j = k + 1; j < n; ++j) {
            const cplx a = A(k, j); if (a == cplx(0)) continue;
            cplx* cj = A.col(j); const cplx* ck = A.col(k);
            for (std::size_t i = k + 1; i < n; ++i) cj[i] -= ck[i] * a;
        }
    }
}

void lu_solve(const Matrix& LU, const std::vector<std::size_t>& piv, cplx* b) {
    const std::size_t n = LU.rows;
    for (std::size_t k = 0; k < n; ++k) if (piv[k] != k) std::swap(b[k], b[piv[k]]);
    for (std::size_t j = 0; j < n; ++j) { const cplx bj = b[j]; const cplx* c = LU.col(j); for (std::size_t i = j + 1; i < n; ++i) b[i] -= c[i] * bj; }
    for (std::size_t j = n; j-- > 0;) { b[j] /= LU(j, j); const cplx bj = b[j]; const cplx* c = LU.col(j); for (std::size_t i = 0; i < j; ++i) b[i] -= c[i] * bj; }
}

}  // namespace cbem
