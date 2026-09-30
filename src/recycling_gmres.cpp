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
