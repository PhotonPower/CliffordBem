#include "cbem/kernel/triangle_integrals.hpp"
#include <cmath>

namespace cbem {

void triangle_integrals(const Vec3& x, const std::array<Vec3, 3>& p, const Vec3& n, Vec3& Igrad, real& Iinv) {
    const real w = dot(x - p[0], n);
    const Vec3 rho = x - n * w;
    const real aw = std::abs(w);
    Vec3 gpar{0, 0, 0}; real inv = 0, bsum = 0;
    for (int e = 0; e < 3; ++e) {
        const Vec3& q1 = p[e]; const Vec3& q2 = p[(e + 1) % 3];
        Vec3 d = q2 - q1; real L = norm(d); Vec3 t = d / L; Vec3 m = cross(t, n);
        real lm = dot(q1 - rho, t), lp = dot(q2 - rho, t), P0 = dot(q1 - rho, m);
        real Rm = norm(x - q1), Rp = norm(x - q2), R0sq = P0 * P0 + w * w;
        real f = (lm + lp >= 0) ? std::log((Rp + lp) / (Rm + lm)) : std::log((Rm - lm) / (Rp - lp));
        real beta = std::atan2(P0 * lp, R0sq + aw * Rp) - std::atan2(P0 * lm, R0sq + aw * Rm);
        gpar += m * f; inv += P0 * f; bsum += beta;
    }
    inv -= aw * bsum;
    const real sgn = (w > 0) ? 1.0 : (w < 0 ? -1.0 : 0.0);
    Igrad = gpar + n * (sgn * bsum);
    Iinv = inv;
}

}  // namespace cbem
