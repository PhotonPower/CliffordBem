// Krylov-Recycling (v0.28) und GCRO-DR (v0.29): gleiche Loesungen wie GMRES; bei einer Matrix mit wenigen kleinen
// Ausreisser-Eigenwerten sinken die Iterationen folgender rechter Seiten deutlich (die langsamen Richtungen stecken im
// gespeicherten Unterraum). GCRO-DR mit k = 10 findet die 8 Ausreisser als harmonische Ritz-Vektoren.
#include <algorithm>
#include <cmath>
#include <random>
#include "cbem/linalg/dense.hpp"
#include "cbem/solvers/recycling_gmres.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    {   // Eigenloeser fuer kleine komplexe Matrizen (harmonische Ritz-Werte)
        std::mt19937 g0(3); std::normal_distribution<real> d0(0, 1);
        Matrix A(40, 40); for (auto& v : A.a) v = cplx(d0(g0), d0(g0));
        std::vector<cplx> l; Matrix V; const bool ok = eig_complex(A, l, V); real err = 0;
        for (int i = 0; i < 40; ++i) { real e = 0; for (int r = 0; r < 40; ++r) { cplx s = 0; for (int c = 0; c < 40; ++c) s += A(r, c) * V(c, i); e += std::norm(s - l[i] * V(r, i)); } err = std::max(err, std::sqrt(e)); }
        std::printf("  Eigenloeser (40 x 40): max |A v - lambda v| = %.1e\n", err);
        CHECK(ok && err < 1e-10, "Eigenloeser ungenau");
    }
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
    GcroDr Dr(op, nullptr, 10, 40);
    int it_plain = 0, it_rec_later = 0, it_dr_later = 0; real maxdiff = 0;
    for (int s = 0; s < 4; ++s) {
        std::vector<cplx> b(n); for (auto& v : b) v = cplx(nd(rng), nd(rng));
        std::vector<cplx> x1, x2, x3;
        const auto g1 = gmres(op, b, x1, nullptr, 1e-10, 300, 2000);
        const auto g2 = R.solve(b, x2, 1e-10, 300, 2000);
        const auto g3 = Dr.solve(b, x3, 1e-10, 0, 2000);
        real d = 0, d3 = 0, nx = 0; for (std::size_t i = 0; i < n; ++i) { d += std::norm(x1[i] - x2[i]); d3 += std::norm(x1[i] - x3[i]); nx += std::norm(x1[i]); }
        maxdiff = std::max({maxdiff, std::sqrt(d / nx), std::sqrt(d3 / nx)});
        std::printf("  rechte Seite %d: GMRES %d It., Recycling %d It. (Unterraum %zu), GCRO-DR(10, 40) %d It., rel. Abweichung %.1e / %.1e\n",
                    s, g1.iterations, g2.iterations, R.recycled(), g3.iterations, std::sqrt(d / nx), std::sqrt(d3 / nx));
        CHECK(g2.converged && g3.converged, "Recycling konvergiert nicht");
        if (s > 0) { it_plain += g1.iterations; it_rec_later += g2.iterations; it_dr_later += g3.iterations; }
    }
    CHECK(maxdiff < 1e-8, "Recycling liefert andere Loesungen");
    CHECK(it_rec_later < 0.7 * it_plain, "Recycling spart keine Iterationen bei Ausreisser-Eigenwerten");
    CHECK(it_dr_later < 0.7 * it_plain, "GCRO-DR spart keine Iterationen bei Ausreisser-Eigenwerten");
    REPORT();
}
