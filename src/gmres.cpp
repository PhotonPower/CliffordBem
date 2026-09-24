#include "cbem/solvers/gmres.hpp"
#include <cmath>

namespace cbem {

namespace {
cplx dotc(const std::vector<cplx>& a, const std::vector<cplx>& b) { cplx s = 0; for (std::size_t i = 0; i < a.size(); ++i) s += std::conj(a[i]) * b[i]; return s; }
real nrm(const std::vector<cplx>& a) { return std::sqrt(std::real(dotc(a, a))); }
}  // namespace

GmresResult gmres(const LinOp& A, const std::vector<cplx>& b, std::vector<cplx>& x, const LinOp* M, real tol, int restart, int max_iter) {
    const std::size_t n = b.size();
    if (x.size() != n) x.assign(n, cplx(0));
    GmresResult res; const real bn = nrm(b);
    if (bn == 0) { x.assign(n, cplx(0)); res.converged = true; return res; }
    std::vector<cplx> r(n), w, z;
    auto applyAM = [&](const std::vector<cplx>& v, std::vector<cplx>& out) {
        if (M) { (*M)(v, z); A(z, out); } else A(v, out);
    };
    int total = 0;
    // Startresiduum (x als Loesung in der unvorkonditionierten Variablen)
    while (total < max_iter) {
        A(x, w); for (std::size_t i = 0; i < n; ++i) r[i] = b[i] - w[i];
        real beta = nrm(r); res.rel_residual = beta / bn;
        if (res.rel_residual <= tol) { res.converged = true; break; }
        std::vector<std::vector<cplx>> Vb; Vb.reserve(restart + 1);
        Vb.push_back(r); for (auto& v : Vb[0]) v /= beta;
        std::vector<std::vector<cplx>> H(restart + 1, std::vector<cplx>(restart, cplx(0)));
        std::vector<cplx> cs(restart), sn(restart), g(restart + 1, cplx(0)); g[0] = beta;
        int k = 0;
        for (; k < restart && total < max_iter; ++k, ++total) {
            applyAM(Vb[k], w);
            for (int j = 0; j <= k; ++j) { H[j][k] = dotc(Vb[j], w); for (std::size_t i = 0; i < n; ++i) w[i] -= H[j][k] * Vb[j][i]; }
            real hn = nrm(w); H[k + 1][k] = hn;
            for (int j = 0; j < k; ++j) {                     // bisherige Rotationen
                cplx t = std::conj(cs[j]) * H[j][k] + std::conj(sn[j]) * H[j + 1][k];
                H[j + 1][k] = -sn[j] * H[j][k] + cs[j] * H[j + 1][k]; H[j][k] = t;
            }
            cplx a = H[k][k], bb = H[k + 1][k]; real den = std::sqrt(std::norm(a) + std::norm(bb));
            cs[k] = a / den; sn[k] = bb / den;
            H[k][k] = std::conj(cs[k]) * a + std::conj(sn[k]) * bb; H[k + 1][k] = 0;
            g[k + 1] = -sn[k] * g[k]; g[k] = std::conj(cs[k]) * g[k];
            res.rel_residual = std::abs(g[k + 1]) / bn;
            if (hn > 0) { Vb.emplace_back(w); for (auto& v : Vb.back()) v /= hn; } else { ++k; ++total; break; }
            if (res.rel_residual <= tol) { ++k; ++total; break; }
        }
        // Rueckwaertseinsetzen
        std::vector<cplx> y(k);
        for (int i = k - 1; i >= 0; --i) { cplx s = g[i]; for (int j = i + 1; j < k; ++j) s -= H[i][j] * y[j]; y[i] = s / H[i][i]; }
        std::vector<cplx> u(n, cplx(0));
        for (int j = 0; j < k; ++j) for (std::size_t i = 0; i < n; ++i) u[i] += y[j] * Vb[j][i];
        if (M) (*M)(u, z); else z = u;
        for (std::size_t i = 0; i < n; ++i) x[i] += z[i];
        res.iterations = total;
        if (res.rel_residual <= tol) {
            A(x, w); real rr = 0; for (std::size_t i = 0; i < n; ++i) rr += std::norm(b[i] - w[i]);
            res.rel_residual = std::sqrt(rr) / bn; res.converged = res.rel_residual <= 10 * tol; break;
        }
    }
    res.iterations = total;
    return res;
}

}  // namespace cbem
