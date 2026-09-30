#include "cbem/sources/chiral_incidence.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cbem {

namespace {
Multivector projector(int s) { return Multivector::blade(0, 0.5) + Multivector::blade(7, cplx(0, 0.5 * s)); }   // P_s = (1 + s iI)/2
real proj_norm2(const Multivector& a) { real s = 0; for (const auto& c : a.c) s += std::norm(c); return s; }
}  // namespace

PlaneWaveIncidence plane_wave_incidence(const Medium& m, real omega, const Vec3& d, const CVec3& p) {
    if (!(std::abs(m.chi) > 0)) return {m.k(omega), 0};
    const Vec3 dn = d / norm(d);
    const CVec3 dxp{dn.y * p[2] - dn.z * p[1], dn.z * p[0] - dn.x * p[2], dn.x * p[1] - dn.y * p[0]};
    const Multivector F = Multivector::vector(p) + Multivector::blade(7) * Multivector::vector(dxp);   // p + I d x p
    const real np = proj_norm2(projector(+1) * F), nm = proj_norm2(projector(-1) * F);
    if (std::min(np, nm) > 1e-16 * std::max(np, nm))
        throw std::invalid_argument("chirales Aussenmedium: nur Helizitaetswellen (zirkulare Polarisation) sind Eigenmoden");
    const int s = np > nm ? +1 : -1;
    return {m.k(omega, s), s};
}

std::vector<cplx> helicity_part(const std::vector<cplx>& h, int sigma) {
    if (sigma == 0) return h;
    const Multivector P = projector(sigma);
    std::vector<cplx> out(h.size());
    for (std::size_t t = 0; t < h.size() / 8; ++t) {
        Multivector F; for (int c = 0; c < 8; ++c) F.c[c] = h[8 * t + c];
        const Multivector G = P * F; for (int c = 0; c < 8; ++c) out[8 * t + c] = G.c[c];
    }
    return out;
}

real extinction_in_medium(const TriangleMesh& mesh, const std::vector<cplx>& hs, const Medium& m, const PlaneWaveIncidence& inc,
                          const Vec3& d, const CVec3& p) {
    return extinction_cross_section(mesh, helicity_part(hs, inc.proj), inc.k, m.eps, d, p);
}

cplx forward_amplitude_in_medium(const TriangleMesh& mesh, const std::vector<cplx>& hs, const Medium& m, const PlaneWaveIncidence& inc,
                                 const Vec3& d, const CVec3& p) {
    return forward_amplitude(mesh, helicity_part(hs, inc.proj), inc.k, m.eps, d, p);
}

}  // namespace cbem
