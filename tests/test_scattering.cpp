// Streuung an der Ikosaeder-Kugel: Vergleich mit dem Python-Prototyp (dichte Loesung, gleiche Diskretisierung).
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/solvers/gmres.hpp"
#include "cbem/sources/fields.hpp"
#include "check.hpp"
using namespace cbem;
static real qext(int n, real om, cplx eps1) {
    TriangleMesh m = make_icosphere(n);
    Medium in{eps1, 1.0}, out{1.0, 1.0};
    cplx k1 = om * std::sqrt(eps1), k2 = om;
    KernelEntries Ein(m, k1), Eout(m, k2); HMatrixParams p; p.eps = 1e-8;
    KernelHMatrix H1(Ein, p), H2(Eout, p); CauchyOperator E1(m, H1), E2(m, H2);
    TransmissionOperator T(m, E1, E2, in, out);
    Vec3 d(0, 0, 1), pol(1, 0, 0);
    auto b = project_plane_wave(m, k2, 1.0, d, pol);
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.precondition(x, y); };
    std::vector<cplx> h; GmresResult r = gmres(A, b, h, &M, 1e-10, 300, 1000);
    std::vector<cplx> hs(h.size()); for (std::size_t i = 0; i < h.size(); ++i) hs[i] = h[i] - b[i];
    real s = extinction_cross_section(m, hs, k2, 1.0, d, pol);
    std::printf("  n=%d, eps1=%g%+gi: Q_ext = %.6f  (%d It., Residuum %.1e)\n", n, eps1.real(), eps1.imag(), s / pi, r.iterations, r.rel_residual);
    return s / pi;
}
int main() {
    // Referenzwerte: Python-Prototyp ap2/scatter3d.py (dichtes LU), gleiche Netze und Quadraturen
    real q1 = qext(4, 1.0, 2.25);          CHECK(std::abs(q1 - 0.20291) < 2e-5, "Glas n=4: %.6f statt 0.20291", q1);
    real q2 = qext(4, 0.5, cplx(-11, 1.2)); CHECK(std::abs(q2 - 0.66915) < 5e-5, "Gold n=4: %.6f statt 0.66915", q2);
    REPORT();
}
