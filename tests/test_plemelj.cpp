// Spurtrennung: E h = +h fuer Spuren innerer, E h = -h fuer Spuren aeusserer Dirac-Loesungen.
#include "cbem/operators/cauchy_operator.hpp"
#include "cbem/kernel/dirac_kernel.hpp"
#include "check.hpp"
using namespace cbem;
static Multivector dirac_solution(const Vec3& x, const Vec3& xs, cplx k, const Multivector& c) {
    Vec3 z = x - xs; KernelValue kv = dirac_kernel_full(z, k);
    Multivector Phi = Multivector::blade(0, kv.s) + Multivector::vector(CVec3{kv.vcoef * z.x, kv.vcoef * z.y, kv.vcoef * z.z});
    return Phi * c;
}
int main() {
    TriangleMesh m = make_icosphere(4);
    MeshQuadrature q(m, QuadRule::subdivided(2));
    Multivector c; for (int b = 0; b < 8; ++b) c.c[b] = cplx(0.3 + 0.1 * b, -0.2 + 0.05 * b);
    for (cplx k : {cplx(0), cplx(1.5), cplx(1.2, 0.4)}) {
        KernelEntries E(m, k); HMatrixParams p; p.eps = 1e-8; KernelHMatrix H(E, p); CauchyOperator op(m, H);
        for (int side : {+1, -1}) {
            Vec3 xs = side > 0 ? Vec3(2.3, -0.8, 0.9) : Vec3(0.1, 0.25, -0.2);
            std::vector<cplx> h(8 * m.size(), 0.0);
            for (std::size_t t = 0; t < m.size(); ++t) {
                for (int a = 0; a < q.q; ++a) { Multivector v = dirac_solution(q.points(t)[a], xs, k, c); for (int b = 0; b < 8; ++b) h[8 * t + b] += q.weights(t)[a] * v.c[b]; }
                for (int b = 0; b < 8; ++b) h[8 * t + b] /= std::sqrt(m.area[t]);
            }
            std::vector<cplx> y; op.apply(h, y);
            real num = 0, den = 0; for (std::size_t i = 0; i < h.size(); ++i) { num += std::norm(y[i] - real(side) * h[i]); den += std::norm(h[i]); }
            real err = std::sqrt(num / den);
            std::printf("  k = %g%+gi, %s: |E h %s h|/|h| = %.2e\n", k.real(), k.imag(), side > 0 ? "innen" : "aussen", side > 0 ? "-" : "+", err);
            CHECK(err < 0.03, "Spurtrennung verletzt (%.2e)", err);
        }
    }
    REPORT();
}
