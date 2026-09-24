#include "cbem/linalg/dense.hpp"
#include <random>
#include "check.hpp"
using namespace cbem;
int main() {
    std::mt19937 g(1); std::normal_distribution<real> nd;
    for (auto [m, n] : std::vector<std::pair<int, int>>{{10, 6}, {30, 30}, {50, 7}}) {
        Matrix A(m, n); for (auto& v : A.a) v = cplx(nd(g), nd(g));
        Matrix Q, R; qr_cgs2(A, Q, R);
        real err = 0, orth = 0;
        for (int i = 0; i < m; ++i) for (int j = 0; j < n; ++j) { cplx s = 0; for (int k = 0; k < n; ++k) s += Q(i, k) * R(k, j); err = std::max(err, std::abs(s - A(i, j))); }
        for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) { cplx s = 0; for (int i = 0; i < m; ++i) s += std::conj(Q(i, a)) * Q(i, b); orth = std::max(orth, std::abs(s - (a == b ? 1.0 : 0.0))); }
        CHECK(err < 1e-12 && orth < 1e-12, "QR %dx%d: Rekonstruktion %.1e, Orthogonalitaet %.1e", m, n, err, orth);
        Matrix U, V; std::vector<real> s; svd_jacobi(A, U, s, V);
        real e2 = 0;
        for (int i = 0; i < m; ++i) for (int j = 0; j < n; ++j) { cplx t = 0; for (int k = 0; k < n; ++k) t += U(i, k) * s[k] * std::conj(V(j, k)); e2 = std::max(e2, std::abs(t - A(i, j))); }
        bool sorted = true; for (int k = 1; k < n; ++k) sorted &= s[k] <= s[k - 1] + 1e-14;
        CHECK(e2 < 1e-11 && sorted, "SVD %dx%d: Rekonstruktion %.1e", m, n, e2);
    }
    REPORT();
}
