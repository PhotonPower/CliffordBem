#include "cbem/solvers/recycling_gmres.hpp"
#include <cmath>

namespace cbem {

namespace {
cplx dotc(const std::vector<cplx>& a, const std::vector<cplx>& b) { cplx s = 0; for (std::size_t i = 0; i < a.size(); ++i) s += std::conj(a[i]) * b[i]; return s; }
real nrm(const std::vector<cplx>& a) { return std::sqrt(std::real(dotc(a, a))); }
void axpy(cplx a, const std::vector<cplx>& x, std::vector<cplx>& y) { for (std::size_t i = 0; i < y.size(); ++i) y[i] += a * x[i]; }
}  // namespace

void RecyclingGmres::applyB(const std::vector<cplx>& v, std::vector<cplx>& out) {
    if (hasM_) { M_(v, tmp_); A_(tmp_, out); } else A_(v, out);
}

GmresResult RecyclingGmres::solve(const std::vector<cplx>& b, std::vector<cplx>& x, real tol, int restart, int max_iter) {
    const std::size_t n = b.size();
    GmresResult res; const real bn = nrm(b);
    std::vector<cplx> z(n, cplx(0)), r = b, w;
    if (bn == 0) { x.assign(n, cplx(0)); res.converged = true; return res; }
    // 1. Projektion auf den gespeicherten Unterraum (zweifach)
    for (int pass = 0; pass < 2; ++pass)
        for (std::size_t i = 0; i < C_.size(); ++i) { const cplx c = dotc(C_[i], r); axpy(-c, C_[i], r); axpy(c, U_[i], z); }
    int total = 0;
    while (true) {
        const real beta = nrm(r); res.rel_residual = beta / bn;
        if (res.rel_residual <= tol || total >= max_iter) break;
        const int m = restart;
        const std::size_t k = C_.size();
        std::vector<std::vector<cplx>> V; V.reserve(m + 1);
        V.push_back(r); for (auto& v : V[0]) v /= beta;
        std::vector<std::vector<cplx>> Hraw(m + 1, std::vector<cplx>(m, cplx(0))), H = Hraw, E(k, std::vector<cplx>(m, cplx(0)));
        std::vector<cplx> cs(m), sn(m), g(m + 1, cplx(0)); g[0] = beta;
        int j = 0; bool stop = false;
        for (; j < m && total < max_iter; ++j, ++total) {
            applyB(V[j], w);
            for (int pass = 0; pass < 2; ++pass)                               // (I - C C^H) B v, zweifach
                for (std::size_t i = 0; i < k; ++i) { const cplx e = dotc(C_[i], w); E[i][j] += e; axpy(-e, C_[i], w); }
            for (int l = 0; l <= j; ++l) { const cplx h = dotc(V[l], w); Hraw[l][j] = h; axpy(-h, V[l], w); }
            for (int l = 0; l <= j; ++l) { const cplx h = dotc(V[l], w); Hraw[l][j] += h; axpy(-h, V[l], w); }   // Nachorthogonalisierung
            const real hn = nrm(w); Hraw[j + 1][j] = hn;
            for (int l = 0; l <= j + 1; ++l) H[l][j] = Hraw[l][j];
            for (int l = 0; l < j; ++l) {                                      // bisherige Givens-Rotationen
                const cplx t = std::conj(cs[l]) * H[l][j] + std::conj(sn[l]) * H[l + 1][j];
                H[l + 1][j] = -sn[l] * H[l][j] + cs[l] * H[l + 1][j]; H[l][j] = t;
            }
            const cplx a = H[j][j], bb = H[j + 1][j]; const real den = std::sqrt(std::norm(a) + std::norm(bb));
            cs[j] = a / den; sn[j] = bb / den;
            H[j][j] = std::conj(cs[j]) * a + std::conj(sn[j]) * bb; H[j + 1][j] = 0;
            g[j + 1] = -sn[j] * g[j]; g[j] = std::conj(cs[j]) * g[j];
            res.rel_residual = std::abs(g[j + 1]) / bn;
            if (hn > 0) { V.emplace_back(w); for (auto& v : V.back()) v /= hn; } else { stop = true; ++j; ++total; break; }
            if (res.rel_residual <= tol) { ++j; ++total; break; }
        }
        const int mc = j;                                                    // Spalten dieses Zyklus
        // 3. y aus dem rotierten Dreieckssystem, z += (V - U E) y, r -= V_{m+1} Hraw y
        std::vector<cplx> y(mc);
        for (int i = mc - 1; i >= 0; --i) { cplx s = g[i]; for (int l = i + 1; l < mc; ++l) s -= H[i][l] * y[l]; y[i] = s / H[i][i]; }
        for (int l = 0; l < mc; ++l) axpy(y[l], V[l], z);
        for (std::size_t i = 0; i < k; ++i) { cplx s = 0; for (int l = 0; l < mc; ++l) s += E[i][l] * y[l]; axpy(-s, U_[i], z); }
        const int rows = std::min<int>(mc + 1, static_cast<int>(V.size()));
        for (int l = 0; l < rows; ++l) { cplx s = 0; for (int q = 0; q < mc; ++q) s += Hraw[l][q] * y[q]; axpy(-s, V[l], r); }
        // 4. Unterraum erweitern: Hraw = Q R (Gram-Schmidt auf den Spalten), C_neu = V Q, U_neu = (V - U E) R^{-1}
        if (C_.size() + mc <= max_recycle_ && mc > 0) {
            std::vector<std::vector<cplx>> Q(mc, std::vector<cplx>(rows, cplx(0))); std::vector<std::vector<cplx>> R(mc, std::vector<cplx>(mc, cplx(0)));
            bool ok = true;
            for (int q = 0; q < mc && ok; ++q) {
                std::vector<cplx> col(rows); for (int l = 0; l < rows; ++l) col[l] = Hraw[l][q];
                for (int pass = 0; pass < 2; ++pass)
                    for (int p = 0; p < q; ++p) { cplx s = 0; for (int l = 0; l < rows; ++l) s += std::conj(Q[p][l]) * col[l]; R[p][q] += s; for (int l = 0; l < rows; ++l) col[l] -= s * Q[p][l]; }
                real cn = 0; for (auto& c : col) cn += std::norm(c); cn = std::sqrt(cn);
                if (!(cn > 1e-14 * bn)) { ok = false; break; }
                R[q][q] = cn; for (int l = 0; l < rows; ++l) Q[q][l] = col[l] / cn;
            }
            if (ok) {
                // W = V_m - U E (Spalten), dann U_neu = W R^{-1} (Rueckwaertseinsetzen spaltenweise: U_neu_q = (W_q - sum_{p<q} U_neu_p R_pq) / R_qq)
                std::vector<std::vector<cplx>> Unew(mc), Cnew(mc);
                for (int q = 0; q < mc; ++q) {
                    std::vector<cplx> Wq = V[q];
                    for (std::size_t i = 0; i < k; ++i) axpy(-E[i][q], U_[i], Wq);
                    for (int p = 0; p < q; ++p) axpy(-R[p][q], Unew[p], Wq);
                    for (auto& v : Wq) v /= R[q][q];
                    Unew[q] = std::move(Wq);
                    Cnew[q].assign(n, cplx(0)); for (int l = 0; l < rows; ++l) axpy(Q[q][l], V[l], Cnew[q]);
                }
                for (int q = 0; q < mc; ++q) { U_.push_back(std::move(Unew[q])); C_.push_back(std::move(Cnew[q])); }
            }
        }
        if (stop) break;
    }
    // x = M z; wahres Residuum
    if (hasM_) M_(z, x); else x = z;
    A_(x, w); real rr = 0; for (std::size_t i = 0; i < n; ++i) rr += std::norm(b[i] - w[i]);
    res.rel_residual = std::sqrt(rr) / bn; res.iterations = total; res.converged = res.rel_residual <= 10 * tol;
    return res;
}

}  // namespace cbem

// ---------------------------------------------------------------------------------------------------------------------------
// GCRO-DR
// ---------------------------------------------------------------------------------------------------------------------------
#include "cbem/linalg/dense.hpp"
#include <algorithm>
#include <numeric>

namespace cbem {

namespace {
cplx dotc2(const std::vector<cplx>& a, const std::vector<cplx>& b) { cplx s = 0; for (std::size_t i = 0; i < a.size(); ++i) s += std::conj(a[i]) * b[i]; return s; }
void axpy2(cplx a, const std::vector<cplx>& x, std::vector<cplx>& y) { for (std::size_t i = 0; i < y.size(); ++i) y[i] += a * x[i]; }
// Kleinste Quadrate min |g - G y| ueber QR (CGS2); gibt das Residuum zurueck
real lsq(const Matrix& G, const std::vector<cplx>& g, std::vector<cplx>& y) {
    Matrix Q, R; qr_cgs2(G, Q, R);
    const std::size_t nc = G.cols; y.assign(nc, cplx(0));
    std::vector<cplx> qg(nc); for (std::size_t c = 0; c < nc; ++c) { cplx s = 0; for (std::size_t r = 0; r < G.rows; ++r) s += std::conj(Q(r, c)) * g[r]; qg[c] = s; }
    for (long i = static_cast<long>(nc) - 1; i >= 0; --i) { cplx s = qg[i]; for (std::size_t j = i + 1; j < nc; ++j) s -= R(i, j) * y[j]; y[i] = s / R(i, i); }
    real res = 0; for (std::size_t r = 0; r < G.rows; ++r) { cplx s = g[r]; for (std::size_t c = 0; c < nc; ++c) s -= G(r, c) * y[c]; res += std::norm(s); }
    return std::sqrt(res);
}
}  // namespace

void GcroDr::applyB(const std::vector<cplx>& v, std::vector<cplx>& out) {
    if (hasM_) { M_(v, tmp_); A_(tmp_, out); } else A_(v, out);
}

GmresResult GcroDr::solve(const std::vector<cplx>& b, std::vector<cplx>& x, real tol, int /*restart*/, int max_iter) {
    const std::size_t n = b.size();
    GmresResult res; const real bn = std::sqrt(std::real(dotc2(b, b)));
    std::vector<cplx> z(n, cplx(0)), r = b, w;
    if (bn == 0) { x.assign(n, cplx(0)); res.converged = true; return res; }
    for (int pass = 0; pass < 2; ++pass)                                        // Projektion auf den Unterraum
        for (std::size_t i = 0; i < C_.size(); ++i) { const cplx c = dotc2(C_[i], r); axpy2(-c, C_[i], r); axpy2(c, U_[i], z); }
    int total = 0;
    while (true) {
        const real beta = std::sqrt(std::real(dotc2(r, r))); res.rel_residual = beta / bn;
        if (res.rel_residual <= tol || total >= max_iter) break;
        const std::size_t k = C_.size();
        const int steps = std::max(10, m_ - static_cast<int>(k));
        // U normieren: U~ = U D, B U~ = C D
        std::vector<std::vector<cplx>> Ut(k); std::vector<real> D(k);
        for (std::size_t i = 0; i < k; ++i) { const real un = std::sqrt(std::real(dotc2(U_[i], U_[i]))); D[i] = 1 / un; Ut[i] = U_[i]; for (auto& v : Ut[i]) v *= D[i]; }
        std::vector<std::vector<cplx>> V; V.reserve(steps + 1); V.push_back(r); for (auto& v : V[0]) v /= beta;
        std::vector<std::vector<cplx>> Hraw(steps + 1, std::vector<cplx>(steps, cplx(0))), E(k, std::vector<cplx>(steps, cplx(0)));
        std::vector<cplx> y; Matrix G; int mc = 0; bool breakdown = false;
        auto build_G = [&](int cols) {                                          // (k + cols + 1) x (k + cols)
            const int rows = static_cast<int>(k) + cols + (breakdown ? 0 : 1);
            Matrix Gm(rows, k + cols);
            for (std::size_t i = 0; i < k; ++i) { Gm(i, i) = D[i]; for (int j = 0; j < cols; ++j) Gm(i, k + j) = E[i][j]; }
            for (int l = 0; l < cols + (breakdown ? 0 : 1); ++l) for (int j = 0; j < cols; ++j) Gm(k + l, k + j) = Hraw[l][j];
            return Gm;
        };
        for (int j = 0; j < steps && total < max_iter; ++j) {
            applyB(V[j], w); ++total;
            for (int pass = 0; pass < 2; ++pass)
                for (std::size_t i = 0; i < k; ++i) { const cplx e = dotc2(C_[i], w); E[i][j] += e; axpy2(-e, C_[i], w); }
            for (int pass = 0; pass < 2; ++pass)
                for (int l = 0; l <= j; ++l) { const cplx h = dotc2(V[l], w); Hraw[l][j] += h; axpy2(-h, V[l], w); }
            const real hn = std::sqrt(std::real(dotc2(w, w))); Hraw[j + 1][j] = hn;
            mc = j + 1;
            if (hn > 1e-14 * bn) { V.emplace_back(w); for (auto& v : V.back()) v /= hn; } else breakdown = true;
            G = build_G(mc);
            std::vector<cplx> g(G.rows, cplx(0)); g[k] = beta;
            res.rel_residual = lsq(G, g, y) / bn;
            if (res.rel_residual <= tol || breakdown) break;
        }
        // z += V^ y, r -= W G y
        for (std::size_t i = 0; i < k; ++i) axpy2(y[i], Ut[i], z);
        for (int j = 0; j < mc; ++j) axpy2(y[k + j], V[j], z);
        std::vector<cplx> Gy(G.rows, cplx(0)); for (std::size_t rr = 0; rr < G.rows; ++rr) for (std::size_t c = 0; c < G.cols; ++c) Gy[rr] += G(rr, c) * y[c];
        for (std::size_t i = 0; i < k; ++i) axpy2(-Gy[i], C_[i], r);
        for (std::size_t l = 0; l + k < G.rows; ++l) axpy2(-Gy[k + l], V[l], r);
        // neuer Unterraum aus harmonischen Ritz-Vektoren
        const std::size_t nc = G.cols, nr = G.rows, kk = std::min<std::size_t>(k_, nc);
        auto Wrow = [&](std::size_t l) -> const std::vector<cplx>& { return l < k ? C_[l] : V[l - k]; };
        auto Vcol = [&](std::size_t l) -> const std::vector<cplx>& { return l < k ? Ut[l] : V[l - k]; };
        Matrix P(nc, kk);
        if (kk == nc) { for (std::size_t i = 0; i < nc; ++i) P(i, i) = 1.0; }
        else {
            Matrix WV(nr, nc);                                                  // W^H V^
            for (std::size_t i = 0; i < k; ++i) for (std::size_t c = 0; c < k; ++c) WV(i, c) = dotc2(C_[i], Ut[c]);
            for (std::size_t l = 0; l + k < nr; ++l) for (std::size_t c = 0; c < k; ++c) WV(k + l, c) = dotc2(V[l], Ut[c]);
            for (std::size_t j = 0; j + k < nc; ++j) WV(k + j, k + j) = 1.0;
            Matrix As(nc, nc), Bs(nc, nc);
            for (std::size_t i = 0; i < nc; ++i) for (std::size_t j = 0; j < nc; ++j) {
                cplx sa = 0, sb = 0; for (std::size_t l = 0; l < nr; ++l) { sa += std::conj(G(l, i)) * G(l, j); sb += std::conj(G(l, i)) * WV(l, j); }
                As(i, j) = sa; Bs(i, j) = sb;
            }
            std::vector<std::size_t> piv; lu_factor(As, piv);
            for (std::size_t j = 0; j < nc; ++j) lu_solve(As, piv, Bs.col(j));   // T = As^{-1} Bs
            std::vector<cplx> mu; Matrix Z;
            if (!eig_complex(Bs, mu, Z)) { for (std::size_t i = 0; i < kk; ++i) P(i, i) = 1.0; }
            else {
                std::vector<std::size_t> ord(nc); std::iota(ord.begin(), ord.end(), 0);
                std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t bq) { return std::abs(mu[a]) > std::abs(mu[bq]); });
                for (std::size_t i = 0; i < kk; ++i) for (std::size_t l = 0; l < nc; ++l) P(l, i) = Z(l, ord[i]);
            }
        }
        Matrix GP(nr, kk);
        for (std::size_t rr = 0; rr < nr; ++rr) for (std::size_t i = 0; i < kk; ++i) { cplx s = 0; for (std::size_t c = 0; c < nc; ++c) s += G(rr, c) * P(c, i); GP(rr, i) = s; }
        Matrix Q, R; qr_cgs2(GP, Q, R);
        std::vector<std::vector<cplx>> Cn(kk, std::vector<cplx>(n, cplx(0))), Un(kk, std::vector<cplx>(n, cplx(0)));
        bool ok = true;
        for (std::size_t i = 0; i < kk && ok; ++i) {
            for (std::size_t l = 0; l < nr; ++l) axpy2(Q(l, i), Wrow(l), Cn[i]);
            for (std::size_t l = 0; l < nc; ++l) axpy2(P(l, i), Vcol(l), Un[i]);
            for (std::size_t p = 0; p < i; ++p) axpy2(-R(p, i), Un[p], Un[i]);
            if (!(std::abs(R(i, i)) > 1e-14)) ok = false; else for (auto& v : Un[i]) v /= R(i, i);
        }
        if (ok) { C_ = std::move(Cn); U_ = std::move(Un); }
        if (breakdown) break;
    }
    if (hasM_) M_(z, x); else x = z;
    A_(x, w); real rr = 0; for (std::size_t i = 0; i < n; ++i) rr += std::norm(b[i] - w[i]);
    res.rel_residual = std::sqrt(rr) / bn; res.iterations = total; res.converged = res.rel_residual <= 10 * tol;
    return res;
}

}  // namespace cbem
