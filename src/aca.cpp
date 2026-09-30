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

LowRank aca_plus(const RowFn& row, const ColFn& col, std::size_t m, std::size_t n, real eps, std::size_t rmax) {
    std::vector<std::vector<cplx>> Us, Vs;
    std::vector<char> urow(m, 0), ucol(n, 0);
    const std::size_t kmax = std::min({rmax, m, n});
    // Rest einer Zeile bzw. Spalte gegen die bisherigen Kreuze
    auto res_row = [&](std::size_t i, std::vector<cplx>& r) { r.resize(n); row(i, r.data()); for (std::size_t k = 0; k < Us.size(); ++k) { const cplx a = Us[k][i]; for (std::size_t c = 0; c < n; ++c) r[c] -= a * Vs[k][c]; } };
    auto res_col = [&](std::size_t j, std::vector<cplx>& u) { u.resize(m); col(j, u.data()); for (std::size_t k = 0; k < Us.size(); ++k) { const cplx a = Vs[k][j]; for (std::size_t q = 0; q < m; ++q) u[q] -= Us[k][q] * a; } };
    auto nrm2 = [](const std::vector<cplx>& x) { real s = 0; for (auto v : x) s += std::norm(v); return s; };
    // Referenzen: gleichmaessig verteilte, noch nicht verwendete Indizes (deterministisch), deren Rest nicht verschwindet
    std::size_t nref_r = 0, nref_c = 0;
    auto next_ref = [](std::size_t len, std::size_t& counter, const std::vector<char>& used, std::size_t avoid) {
        for (std::size_t t = 0; t < len; ++t) {
            const std::size_t idx = (counter * 7919 + len / 2) % len; ++counter;
            if (!used[idx] && idx != avoid) return idx;
        }
        for (std::size_t q = 0; q < len; ++q) if (!used[q] && q != avoid) return q;
        return len;
    };
    std::size_t ir = next_ref(m, nref_r, urow, m), jr = next_ref(n, nref_c, ucol, n);
    std::vector<cplx> rref, cref, u, v;
    if (ir < m) res_row(ir, rref);
    if (jr < n) res_col(jr, cref);
    // Masstab fuer "verschwunden": relativ zum groessten bisher gesehenen Eintrag (Rundungsrauschen nach exakter Erfassung ist
    // nicht null; ein absoluter Grenzwert liesse die Referenz auf dem Rauschen pivotieren)
    real scale = 0;
    auto maxabs = [](const std::vector<cplx>& x) { real a = 0; for (auto v0 : x) a = std::max(a, std::abs(v0)); return a; };
    scale = std::max(maxabs(rref), maxabs(cref));
    auto vanished = [&](const std::vector<cplx>& x) { return maxabs(x) <= 1e-12 * scale; };
    real S2 = 0;
    for (std::size_t it = 0; it < kmax; ++it) {
        // Referenzen ersetzen, deren Rest verschwunden ist oder die schon Pivot waren
        for (int tries = 0; ir < m && (urow[ir] || vanished(rref)); ++tries) {
            if (tries == 8) { ir = m; break; }                                  // keine Referenzzeile mit Rest mehr
            ir = next_ref(m, nref_r, urow, ir); if (ir < m) { res_row(ir, rref); scale = std::max(scale, maxabs(rref)); }
        }
        for (int tries = 0; jr < n && (ucol[jr] || vanished(cref)); ++tries) {
            if (tries == 8) { jr = n; break; }
            jr = next_ref(n, nref_c, ucol, jr); if (jr < n) { res_col(jr, cref); scale = std::max(scale, maxabs(cref)); }
        }
        if (ir == m && jr == n) break;                                          // beide Referenzen ohne Rest: Block erfasst
        std::size_t js = n, is = m; real a = -1, b = -1;
        if (ir < m) for (std::size_t c = 0; c < n; ++c) if (!ucol[c] && std::abs(rref[c]) > a) { a = std::abs(rref[c]); js = c; }
        if (jr < n) for (std::size_t q = 0; q < m; ++q) if (!urow[q] && std::abs(cref[q]) > b) { b = std::abs(cref[q]); is = q; }
        if (std::max(a, b) <= 1e-12 * scale) break;
        std::size_t i, j; cplx piv;
        if (a >= b) {                                                           // Spalte aus der Referenzzeile
            j = js; res_col(j, u);
            i = m; real bu = -1; for (std::size_t q = 0; q < m; ++q) if (!urow[q] && std::abs(u[q]) > bu) { bu = std::abs(u[q]); i = q; }
            if (i == m || bu <= 1e-12 * scale) { ucol[j] = 1; continue; }
            res_row(i, v);
        } else {                                                                // Zeile aus der Referenzspalte
            i = is; res_row(i, v);
            j = n; real bv = -1; for (std::size_t c = 0; c < n; ++c) if (!ucol[c] && std::abs(v[c]) > bv) { bv = std::abs(v[c]); j = c; }
            if (j == n || bv <= 1e-12 * scale) { urow[i] = 1; continue; }
            res_col(j, u);
        }
        piv = v[j];
        for (auto& x : v) x /= piv;                                             // u v^T mit v(j) = 1
        const real nu = std::sqrt(nrm2(u)), nv = std::sqrt(nrm2(v));
        real cross = 0;
        for (std::size_t k = 0; k < Us.size(); ++k) {
            cplx pa = 0, pb = 0;
            for (std::size_t q = 0; q < m; ++q) pa += std::conj(Us[k][q]) * u[q];
            for (std::size_t c = 0; c < n; ++c) pb += std::conj(Vs[k][c]) * v[c];
            cross += 2 * std::real(pa * pb);
        }
        S2 += nu * nu * nv * nv + cross;
        Us.push_back(u); Vs.push_back(v); urow[i] = 1; ucol[j] = 1;
        if (ir < m) { const cplx c0 = u[ir]; for (std::size_t c = 0; c < n; ++c) rref[c] -= c0 * v[c]; }   // Referenzreste nachfuehren
        if (jr < n) { const cplx c0 = v[jr]; for (std::size_t q = 0; q < m; ++q) cref[q] -= u[q] * c0; }
        const real S = std::sqrt(std::abs(S2));
        const real er = ir < m ? std::sqrt(nrm2(rref) * m) : 0, ec = jr < n ? std::sqrt(nrm2(cref) * n) : 0;   // auf den Block hochgerechnet
        if (nu * nv <= eps * S && er <= eps * S && ec <= eps * S) break;
    }
    LowRank lr; lr.U = Matrix(m, Us.size()); lr.V = Matrix(n, Vs.size());
    for (std::size_t k = 0; k < Us.size(); ++k) { std::copy(Us[k].begin(), Us[k].end(), lr.U.col(k)); std::copy(Vs[k].begin(), Vs[k].end(), lr.V.col(k)); }
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
