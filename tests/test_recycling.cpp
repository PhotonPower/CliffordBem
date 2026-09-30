// Krylov-Recycling (v0.28): gleiche Loesungen wie GMRES; bei einer Matrix mit wenigen kleinen Ausreisser-Eigenwerten sinken
// die Iterationen folgender rechter Seiten deutlich (die langsamen Richtungen stecken im gespeicherten Unterraum).
#include <cmath>
#include <random>
#include "cbem/solvers/recycling_gmres.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    const std::size_t n = 400;
    std::mt19937 rng(7); std::normal_distribution<real> nd(0, 1);
    // A = Q diag(lambda) Q^H + kleine Stoerung; lambda: 8 Ausreisser nahe 0, Rest in einem Kreis um 1
    std::vector<std::vector<cplx>> Q(n, std::vector<cplx>(n));
    for (auto& r : Q) for (auto& v : r) v = cplx(nd(rng), nd(rng));
    for (std::size_t j = 0; j < n; ++j) {                                 // Gram-Schmidt auf den Spalten
        for (std::size_t p = 0; p < j; ++p) { cplx s = 0; for (std::size_t i = 0; i < n; ++i) s += std::conj(Q[i][p]) * Q[i][j]; for (std::size_t i = 0; i < n; ++i) Q[i][j] -= s * Q[i][p]; }
        real s = 0; for (std::size_t i = 0; i < n; ++i) s += std::norm(Q[i][j]); s = std::sqrt(s); for (std::size_t i = 0; i < n; ++i) Q[i][j] /= s;
    }
    std::vector<cplx> lam(n);
    for (std::size_t i = 0; i < n; ++i) lam[i] = i < 8 ? cplx(1e-3 * (i + 1), 1e-3) : cplx(1.0, 0.0) + 0.5 * std::polar(1.0, 2 * pi * i / n) * std::sqrt(real(i % 7) / 7);
    std::vector<std::vector<cplx>> A(n, std::vector<cplx>(n, cplx(0)));
    for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j) { cplx s = 0; for (std::size_t k = 0; k < n; ++k) s += Q[i][k] * lam[k] * std::conj(Q[j][k]); A[i][j] = s; }
    LinOp op = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { y.assign(n, cplx(0)); for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j) y[i] += A[i][j] * x[j]; };
    RecyclingGmres R(op, nullptr, 200);
    int it_plain = 0, it_rec_later = 0; real maxdiff = 0;
    for (int s = 0; s < 4; ++s) {
        std::vector<cplx> b(n); for (auto& v : b) v = cplx(nd(rng), nd(rng));
        std::vector<cplx> x1, x2;
        const auto g1 = gmres(op, b, x1, nullptr, 1e-10, 300, 2000);
        const auto g2 = R.solve(b, x2, 1e-10, 300, 2000);
        real d = 0, nx = 0; for (std::size_t i = 0; i < n; ++i) { d += std::norm(x1[i] - x2[i]); nx += std::norm(x1[i]); }
        maxdiff = std::max(maxdiff, std::sqrt(d / nx));
        std::printf("  rechte Seite %d: GMRES %d It., Recycling %d It. (Unterraum %zu), rel. Abweichung %.1e\n", s, g1.iterations, g2.iterations, R.recycled(), std::sqrt(d / nx));
        CHECK(g2.converged, "Recycling konvergiert nicht");
        if (s > 0) { it_plain += g1.iterations; it_rec_later += g2.iterations; }
    }
    CHECK(maxdiff < 1e-8, "Recycling liefert andere Loesungen");
    CHECK(it_rec_later < 0.7 * it_plain, "Recycling spart keine Iterationen bei Ausreisser-Eigenwerten");
    REPORT();
}
