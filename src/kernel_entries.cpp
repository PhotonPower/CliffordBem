#include "cbem/assembly/kernel_entries.hpp"
#include "cbem/kernel/dirac_kernel.hpp"
#include "cbem/kernel/triangle_integrals.hpp"

namespace cbem {

KernelEntries::KernelEntries(const TriangleMesh& mesh, cplx k, EntryParams prm)
    : m_(mesh), k_(k), prm_(prm), q7_(mesh, QuadRule::dunavant7()), qn_(mesh, QuadRule::subdivided(prm.near_subdivision)) {}

bool KernelEntries::is_near(std::size_t i, std::size_t j) const {
    real d = norm(m_.centroid[i] - m_.centroid[j]);
    return d < prm_.near_factor * std::max(m_.hmax[i], m_.hmax[j]);
}

KernelComp KernelEntries::far(std::size_t i, std::size_t j) const {
    KernelComp K{0, 0, 0, 0};
    const Vec3* xi = q7_.points(i); const real* wi = q7_.weights(i);
    const Vec3* yj = q7_.points(j); const real* wj = q7_.weights(j);
    for (int a = 0; a < q7_.q; ++a)
        for (int b = 0; b < q7_.q; ++b) {
            Vec3 z = xi[a] - yj[b]; KernelValue kv = dirac_kernel_full(z, k_); real ww = wi[a] * wj[b];
            K[0] += ww * kv.s; K[1] += ww * kv.vcoef * z.x; K[2] += ww * kv.vcoef * z.y; K[3] += ww * kv.vcoef * z.z;
        }
    return K;
}

KernelComp KernelEntries::near(std::size_t i, std::size_t j) const {
    KernelComp K{0, 0, 0, 0};
    const bool swap = (i != j) && (m_.hmax[i] > 1.5 * m_.hmax[j]);
    // aeussere Punkte auf dem kleineren Dreieck, inneres Integral analytisch ueber das andere
    std::size_t outer = swap ? j : i, inner = swap ? i : j;
    const Vec3* xo = qn_.points(outer); const real* wo = qn_.weights(outer);
    auto tri = m_.vertices(inner); const Vec3& n = m_.normal[inner];
    const cplx ik = cplx(0, 1) * k_;
    Vec3 g{0, 0, 0}; real inv = 0;
    for (int a = 0; a < qn_.q; ++a) {
        Vec3 Ig; real Ii; triangle_integrals(xo[a], tri, n, Ig, Ii);
        g += Ig * wo[a]; inv += Ii * wo[a];
    }
    // Phi_0-Anteil: ungerade in (x - y); Selbstterm exakt 0
    real sgn = swap ? -1.0 : 1.0;
    if (i != j) { K[1] += sgn * g.x / (4 * pi); K[2] += sgn * g.y / (4 * pi); K[3] += sgn * g.z / (4 * pi); }
    K[0] += -ik * inv / (4 * pi);
    // Rest mit Gauss: aeussere Punkte x (7 Punkte auf i bzw. Nahregel), z = x_i - y_j
    const Vec3* yi = swap ? q7_.points(i) : xo;  const real* wyi = swap ? q7_.weights(i) : wo;  int ni = swap ? q7_.q : qn_.q;
    const Vec3* yj = swap ? xo : q7_.points(j); const real* wyj = swap ? wo : q7_.weights(j); int nj = swap ? qn_.q : q7_.q;
    for (int a = 0; a < ni; ++a)
        for (int b = 0; b < nj; ++b) {
            Vec3 z = yi[a] - yj[b]; KernelValue kv = dirac_kernel_remainder(z, k_); real ww = wyi[a] * wyj[b];
            K[0] += ww * kv.s; K[1] += ww * kv.vcoef * z.x; K[2] += ww * kv.vcoef * z.y; K[3] += ww * kv.vcoef * z.z;
        }
    return K;
}

}  // namespace cbem
