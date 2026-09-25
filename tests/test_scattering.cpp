// Streuung an der Ikosaeder-Kugel: Vergleich mit dem Python-Prototyp (dichte Loesung, gleiche Diskretisierung).
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/solvers/gmres.hpp"
#include "cbem/sources/fields.hpp"
#include "check.hpp"
using namespace cbem;
static real qext(int n, real om, cplx eps1, bool ss) {
    TriangleMesh m = make_icosphere(n);
    Medium in{eps1, 1.0}, out{1.0, 1.0};
    cplx k1 = om * std::sqrt(eps1), k2 = om;
    EntryParams ep; ep.sauter_schwab = ss;
    if (!ss) { ep.adaptive_outer = false; ep.near_subdivision = 4; }   // Nahfeldregel des Python-Prototyps
    KernelEntries Ein(m, k1, ep), Eout(m, k2, ep); HMatrixParams p; p.eps = 1e-8;
    KernelHMatrix H1(Ein, p), H2(Eout, p); CauchyOperator E1(m, H1), E2(m, H2);
    TransmissionOperator T(m, E1, E2, in, out);
    Vec3 d(0, 0, 1), pol(1, 0, 0);
    auto b = project_plane_wave(m, k2, 1.0, d, pol);
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.precondition(x, y); };
    std::vector<cplx> h; GmresResult r = gmres(A, b, h, &M, 1e-10, 300, 1000);
    std::vector<cplx> hs(h.size()); for (std::size_t i = 0; i < h.size(); ++i) hs[i] = h[i] - b[i];
    real s = extinction_cross_section(m, hs, k2, 1.0, d, pol);
    std::printf("  n=%d, eps1=%g%+gi, %s: Q_ext = %.6f  (%d It., Residuum %.1e)\n", n, eps1.real(), eps1.imag(), ss ? "Sauter-Schwab" : "analytisch innen", s / pi, r.iterations, r.rel_residual);
    return s / pi;
}
int main() {
    // Regression: Python-Prototyp ap2/scatter3d.py (dichtes LU) mit derselben Nahfeldmethode
    real q1 = qext(4, 1.0, 2.25, false);          CHECK(std::abs(q1 - 0.20291) < 2e-5, "Glas n=4: %.6f statt 0.20291", q1);
    real q2 = qext(4, 0.5, cplx(-11, 1.2), false); CHECK(std::abs(q2 - 0.66915) < 5e-5, "Gold n=4: %.6f statt 0.66915", q2);
    // Sauter-Schwab: genauere Nahfeldintegrale, Werte naeher an Mie (0.215098 bzw. 0.590018)
    real q3 = qext(4, 1.0, 2.25, true);           CHECK(std::abs(q3 - q1) < 1e-3, "Glas: SS-Aenderung zu gross (Geometriefehler dominiert)");
    real q4 = qext(4, 0.5, cplx(-11, 1.2), true); CHECK(std::abs(q4 - 0.590018) < std::abs(q2 - 0.590018), "Gold: SS nicht naeher an Mie");
    REPORT();
}
