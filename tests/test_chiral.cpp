// Chirale Medien: Transmissionsabbildung gegen den Python-Prototyp (Fresnel-verifiziert, AP 1),
// Reduktion auf den achiralen Fall, Symmetrie sigma_s(chi) = sigma_{-s}(-chi), Vergleich mit chiraler Mie-Loesung.
#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/solvers/gmres.hpp"
#include "cbem/sources/fields.hpp"
#include "check.hpp"
using namespace cbem;
static real sigma(const TriangleMesh& m, real om, const Medium& in, int s, bool chiral_op) {
    const Medium out{1.0, 1.0, 0.0}; const cplx k2 = om;
    HMatrixParams p; p.eps = 1e-8;
    KernelEntries Ep(m, in.k(om, +1)), Em(m, in.k(om, -1)), Eo(m, k2);
    KernelHMatrix Hp(Ep, p), Hm(Em, p), Ho(Eo, p);
    CauchyOperator Cp(m, Hp), Cm(m, Hm), Co(m, Ho); ChiralCauchyOperator Cch(Cp, Cm);
    const BoundaryOperator& E1 = chiral_op ? static_cast<const BoundaryOperator&>(Cch) : static_cast<const BoundaryOperator&>(Cp);
    TransmissionOperator T(m, E1, Co, in, out);
    Vec3 d(0, 0, 1); CVec3 pol = circular_polarization(d, s);
    auto b = project_plane_wave(m, k2, 1.0, d, pol);
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.precondition(x, y); };
    std::vector<cplx> h; gmres(A, b, h, &M, 1e-10, 300, 1000);
    for (std::size_t i = 0; i < h.size(); ++i) h[i] -= b[i];
    return extinction_cross_section(m, h, k2, 1.0, d, pol);
}
int main() {
    // (1) J gegen Python (chiral_J.py), n = (0.36, -0.48, 0.8), innen (2.25, 1.3, 0.3+0.05i), aussen Vakuum
    struct E { int r, c; cplx v; } ref[] = {
        {1,1,cplx(1.394654362531231e+00,9.415439210350452e-04)}, {1,6,cplx(3.642447156055447e-03,-2.051479351024435e-02)},
        {6,2,cplx(3.691581097383703e-03,-2.079152303232346e-02)}, {4,4,cplx(9.797746297838583e-01,4.649599610049607e-03)},
        {3,4,cplx(-1.367252258290261e-02,7.700564086045728e-02)}, {4,3,cplx(1.798739336323678e-02,-1.013076222728116e-01)},
        {0,0,cplx(1.140175425099138e+00,0.0)}, {7,7,cplx(1.5,0.0)}, {2,5,cplx(-6.475461610765239e-03,3.647074401821217e-02)}};
    Mat8 J = transmission_map(Vec3(0.36, -0.48, 0.8), Medium{2.25, 1.3, cplx(0.3, 0.05)}, Medium{1.0, 1.0, 0.0});
    real ej = 0; for (auto& e : ref) ej = std::max(ej, std::abs(J[e.r * 8 + e.c] - e.v));
    std::printf("  J chiral gegen Python: max. Abw. %.1e\n", ej); CHECK(ej < 1e-12, "J chiral weicht ab");
    // (2) chiraler Operator mit chi = 0 = achiraler Operator; (3) Symmetrie; (4) Mie
    TriangleMesh m = make_icosphere(4); const real om = 1.0;
    real s_ach = sigma(m, om, Medium{2.25, 1.0, 0.0}, +1, false), s_ch0 = sigma(m, om, Medium{2.25, 1.0, 0.0}, +1, true);
    std::printf("  chi = 0: achiral %.8f, chiraler Operator %.8f\n", s_ach, s_ch0);
    CHECK(std::abs(s_ach - s_ch0) < 1e-9 * s_ach, "chi = 0 stimmt nicht mit achiral ueberein");
    real sp = sigma(m, om, Medium{2.25, 1.0, 0.2}, +1, true), smm = sigma(m, om, Medium{2.25, 1.0, -0.2}, -1, true);
    real sm = sigma(m, om, Medium{2.25, 1.0, 0.2}, -1, true);
    std::printf("  chi = 0.2: sigma_+ = %.6f, sigma_- = %.6f;  sigma_-(chi=-0.2) = %.6f\n", sp, sm, smm);
    CHECK(std::abs(sp - smm) < 1e-7 * sp, "Symmetrie sigma_s(chi) = sigma_-s(-chi) verletzt");
    // chirale Mie (tools/mie_chiral.py): Q_+ = 0.35557, Q_- = 0.13045; hier grobes Netz (320 Dreiecke), Konvergenz in apps/scatter_chiral
    CHECK(sp > sm, "falsches Vorzeichen des Zirkulardichroismus");
    REPORT();
}
