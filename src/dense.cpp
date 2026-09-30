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

bool eig_complex(const Matrix& A0, std::vector<cplx>& lambda, Matrix& V, int max_iter_per_eig) {
    const std::size_t n = A0.rows;
    Matrix H = A0, Z(n, n);
    for (std::size_t i = 0; i < n; ++i) Z(i, i) = 1.0;
    // Hessenberg-Form mit Householder-Spiegelungen: H = Q^H A Q, Z = Q
    for (std::size_t j = 0; j + 2 < n; ++j) {
        std::vector<cplx> v(n - j - 1); real xn = 0;
        for (std::size_t i = j + 1; i < n; ++i) { v[i - j - 1] = H(i, j); xn += std::norm(H(i, j)); }
        xn = std::sqrt(xn); if (xn == 0) continue;
        const cplx ph = std::abs(v[0]) > 0 ? v[0] / std::abs(v[0]) : cplx(1.0);
        v[0] += ph * xn;
        real vn = 0; for (auto& x : v) vn += std::norm(x); vn = std::sqrt(vn); for (auto& x : v) x /= vn;
        for (std::size_t c = 0; c < n; ++c) {                              // H <- (I - 2vv^H) H
            cplx s = 0; for (std::size_t i = j + 1; i < n; ++i) s += std::conj(v[i - j - 1]) * H(i, c);
            for (std::size_t i = j + 1; i < n; ++i) H(i, c) -= 2.0 * v[i - j - 1] * s;
        }
        for (std::size_t r = 0; r < n; ++r) {                              // H <- H (I - 2vv^H), Z <- Z (I - 2vv^H)
            cplx s = 0, t = 0;
            for (std::size_t i = j + 1; i < n; ++i) { s += H(r, i) * v[i - j - 1]; t += Z(r, i) * v[i - j - 1]; }
            for (std::size_t i = j + 1; i < n; ++i) { H(r, i) -= 2.0 * s * std::conj(v[i - j - 1]); Z(r, i) -= 2.0 * t * std::conj(v[i - j - 1]); }
        }
        for (std::size_t i = j + 2; i < n; ++i) H(i, j) = 0;
    }
    // einfach verschobener QR-Algorithmus auf dem aktiven Block [lo, hi]
    const real eps = 1e-15;
    long hi = static_cast<long>(n) - 1; int it = 0, total = 0;
    while (hi > 0) {
        long l = hi;
        for (; l > 0; --l) if (std::abs(H(l, l - 1)) <= eps * (std::abs(H(l, l)) + std::abs(H(l - 1, l - 1)))) { H(l, l - 1) = 0; break; }
        if (l == hi) { --hi; it = 0; continue; }                            // Eigenwert abgespalten
        if (++it > max_iter_per_eig || ++total > max_iter_per_eig * static_cast<int>(n)) return false;
        // Wilkinson-Shift aus dem 2x2-Block unten rechts; gelegentlich ausnahmsweise Shift
        const cplx a = H(hi - 1, hi - 1), b = H(hi - 1, hi), c = H(hi, hi - 1), d = H(hi, hi);
        const cplx tr = a + d, det = a * d - b * c, disc = std::sqrt(tr * tr - 4.0 * det);
        const cplx e1 = 0.5 * (tr + disc), e2 = 0.5 * (tr - disc);
        cplx mu = std::abs(e1 - d) < std::abs(e2 - d) ? e1 : e2;
        if (it % 11 == 10) mu = d + std::abs(H(hi, hi - 1));
        // QR-Schritt auf [l, hi] mit Givens-Rotationen (implizit ueber den Buckel)
        cplx x = H(l, l) - mu, y = H(l + 1, l);
        for (long k = l; k < hi; ++k) {
            const real r = std::sqrt(std::norm(x) + std::norm(y));
            const cplx cs = r > 0 ? x / r : cplx(1.0), sn = r > 0 ? y / r : cplx(0.0);
            // G = [[conj(cs), conj(sn)], [-sn, cs]] auf die Zeilen k, k+1 (alle Spalten ab max(l, k-1))
            for (std::size_t col = static_cast<std::size_t>(std::max(l, k - 1)); col < n; ++col) {
                const cplx h1 = H(k, col), h2 = H(k + 1, col);
                H(k, col) = std::conj(cs) * h1 + std::conj(sn) * h2; H(k + 1, col) = -sn * h1 + cs * h2;
            }
            for (std::size_t row = 0; row <= static_cast<std::size_t>(std::min(hi, k + 2)); ++row) {   // G^H auf die Spalten k, k+1
                const cplx h1 = H(row, k), h2 = H(row, k + 1);
                H(row, k) = cs * h1 + sn * h2; H(row, k + 1) = -std::conj(sn) * h1 + std::conj(cs) * h2;
            }
            for (std::size_t row = 0; row < n; ++row) {
                const cplx z1 = Z(row, k), z2 = Z(row, k + 1);
                Z(row, k) = cs * z1 + sn * z2; Z(row, k + 1) = -std::conj(sn) * z1 + std::conj(cs) * z2;
            }
            if (k + 1 < hi) { x = H(k + 1, k); y = H(k + 2, k); }
        }
    }
    // Eigenwerte und Eigenvektoren aus der Schur-Form
    lambda.resize(n); for (std::size_t i = 0; i < n; ++i) lambda[i] = H(i, i);
    real tn = 0; for (auto& v : H.a) tn = std::max(tn, std::abs(v));
    V = Matrix(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<cplx> x(n, cplx(0)); x[i] = 1.0;
        for (long j = static_cast<long>(i) - 1; j >= 0; --j) {
            cplx s = 0; for (std::size_t l = j + 1; l <= i; ++l) s += H(j, l) * x[l];
            cplx den = H(j, j) - lambda[i]; if (std::abs(den) < 1e-14 * tn) den = 1e-14 * tn;
            x[j] = -s / den;
        }
        real vn = 0;
        for (std::size_t r = 0; r < n; ++r) { cplx s = 0; for (std::size_t l = 0; l <= i; ++l) s += Z(r, l) * x[l]; V(r, i) = s; vn += std::norm(s); }
        vn = std::sqrt(vn); for (std::size_t r = 0; r < n; ++r) V(r, i) /= vn;
    }
    return true;
}

void lu_solve(const Matrix& LU, const std::vector<std::size_t>& piv, cplx* b) {
    const std::size_t n = LU.rows;
    for (std::size_t k = 0; k < n; ++k) if (piv[k] != k) std::swap(b[k], b[piv[k]]);
    for (std::size_t j = 0; j < n; ++j) { const cplx bj = b[j]; const cplx* c = LU.col(j); for (std::size_t i = j + 1; i < n; ++i) b[i] -= c[i] * bj; }
    for (std::size_t j = n; j-- > 0;) { b[j] /= LU(j, j); const cplx bj = b[j]; const cplx* c = LU.col(j); for (std::size_t i = 0; i < j; ++i) b[i] -= c[i] * bj; }
}

}  // namespace cbem
