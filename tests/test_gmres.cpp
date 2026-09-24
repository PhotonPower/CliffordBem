#include "cbem/solvers/gmres.hpp"
#include <random>
#include "check.hpp"
using namespace cbem;
int main() {
    const int n = 60; std::mt19937 g(5); std::normal_distribution<real> nd;
    std::vector<cplx> A(n * n); for (auto& v : A) v = cplx(nd(g), nd(g)) * 0.1;
    for (int i = 0; i < n; ++i) A[i * n + i] += cplx(2.0, 0.5);
    LinOp op = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { y.assign(n, 0); for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) y[i] += A[i * n + j] * x[j]; };
    LinOp diag = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { y.resize(n); for (int i = 0; i < n; ++i) y[i] = x[i] / A[i * n + i]; };
    std::vector<cplx> xt(n), b; for (auto& v : xt) v = cplx(nd(g), nd(g)); op(xt, b);
    const std::vector<const LinOp*> precs = {nullptr, &diag};
    for (int restart : {10, 100}) for (const LinOp* M : precs) {
        std::vector<cplx> x; GmresResult r = gmres(op, b, x, M, 1e-10, restart, 500);
        real e = 0, s = 0; for (int i = 0; i < n; ++i) { e += std::norm(x[i] - xt[i]); s += std::norm(xt[i]); }
        std::printf("  Neustart %d, %s: %d It., Residuum %.1e, Fehler %.1e\n", restart, M ? "mit Vork." : "ohne     ", r.iterations, r.rel_residual, std::sqrt(e / s));
        CHECK(r.converged && std::sqrt(e / s) < 1e-8, "GMRES nicht konvergiert");
    }
    REPORT();
}
