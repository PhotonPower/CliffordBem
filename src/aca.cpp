#include "cbem/hmatrix/aca.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace cbem {

LowRank aca_partial(const RowFn& row, const ColFn& col, std::size_t m, std::size_t n, real eps, std::size_t rmax) {
    std::vector<std::vector<cplx>> Us, Vs;
    std::vector<char> used(m, 0);
    std::vector<cplx> r(n), u(m);
    std::size_t i = 0; real normS2 = 0;
    const std::size_t kmax = std::min({rmax, m, n});
    for (std::size_t it = 0; it < kmax; ++it) {
        row(i, r.data());
        for (std::size_t k = 0; k < Us.size(); ++k) { cplx ui = Us[k][i]; for (std::size_t c = 0; c < n; ++c) r[c] -= ui * Vs[k][c]; }
        used[i] = 1;
        std::size_t j = 0; real best = -1;
        for (std::size_t c = 0; c < n; ++c) {
            if (std::abs(r[c]) > best) { best = std::abs(r[c]); j = c; }
        }
        if (best <= 1e-300) {
            std::size_t cand = m;
            for (std::size_t q = 0; q < m; ++q) {
                if (!used[q]) { cand = q; break; }
            }
            if (cand == m) break;
            i = cand;
            continue;
        }
        std::vector<cplx> v(n); cplx piv = r[j]; for (std::size_t c = 0; c < n; ++c) v[c] = r[c] / piv;
        col(j, u.data());
        for (std::size_t k = 0; k < Us.size(); ++k) { cplx vj = Vs[k][j]; for (std::size_t q = 0; q < m; ++q) u[q] -= Us[k][q] * vj; }
        real nu = 0, nv = 0; for (auto x : u) nu += std::norm(x); for (auto x : v) nv += std::norm(x); nu = std::sqrt(nu); nv = std::sqrt(nv);
        real cross_terms = 0;
        for (std::size_t k = 0; k < Us.size(); ++k) {
            cplx a = 0, b = 0;
            for (std::size_t q = 0; q < m; ++q) a += std::conj(Us[k][q]) * u[q];
            for (std::size_t c = 0; c < n; ++c) b += std::conj(Vs[k][c]) * v[c];
            cross_terms += 2 * std::real(a * b);
        }
        normS2 += nu * nu * nv * nv + cross_terms;
        Us.push_back(u); Vs.push_back(v);
        if (nu * nv <= eps * std::sqrt(std::abs(normS2))) break;
        std::size_t inext = m; real bu = -1;
        for (std::size_t q = 0; q < m; ++q) if (!used[q] && std::abs(u[q]) > bu) { bu = std::abs(u[q]); inext = q; }
        if (inext == m) break;
        i = inext;
    }
    LowRank lr; lr.U = Matrix(m, Us.size()); lr.V = Matrix(n, Vs.size());
    for (std::size_t k = 0; k < Us.size(); ++k) {
        std::copy(Us[k].begin(), Us[k].end(), lr.U.col(k)); std::copy(Vs[k].begin(), Vs[k].end(), lr.V.col(k));
    }
    return lr;
}

void recompress(LowRank& lr, real eps) {
    const std::size_t r = lr.rank();
    if (r == 0) return;
    Matrix Qu, Ru, Qv, Rv;
    qr_cgs2(lr.U, Qu, Ru); qr_cgs2(lr.V, Qv, Rv);
    Matrix M(r, r);                                   // M = Ru Rv^T
    for (std::size_t a = 0; a < r; ++a)
        for (std::size_t b = 0; b < r; ++b) { cplx s = 0; for (std::size_t c = 0; c < r; ++c) s += Ru(a, c) * Rv(b, c); M(a, b) = s; }
    Matrix W, Z; std::vector<real> s;
    svd_jacobi(M, W, s, Z);                           // M = W diag(s) Z^H
    std::size_t rk = 0; while (rk < r && s[rk] > eps * s[0]) ++rk; rk = std::max<std::size_t>(rk, 1);
    Matrix U(lr.U.rows, rk), V(lr.V.rows, rk);        // U = Qu W S,  V = Qv conj(Z)
    for (std::size_t k = 0; k < rk; ++k) {
        for (std::size_t q = 0; q < lr.U.rows; ++q) { cplx t = 0; for (std::size_t a = 0; a < r; ++a) t += Qu(q, a) * W(a, k); U(q, k) = t * s[k]; }
        for (std::size_t q = 0; q < lr.V.rows; ++q) { cplx t = 0; for (std::size_t a = 0; a < r; ++a) t += Qv(q, a) * std::conj(Z(a, k)); V(q, k) = t; }
    }
    lr.U = std::move(U); lr.V = std::move(V);
}

}  // namespace cbem
