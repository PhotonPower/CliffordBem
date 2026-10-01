#include "cbem/kernel/triangle_integrals.hpp"
#include <algorithm>
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
        // f = int_e dl / R = asinh(l+/R0) - asinh(l-/R0); R0 nach unten begrenzt (Punkte auf der Kantengeraden)
        const real R0 = std::max(std::sqrt(R0sq), 1e-13 * L);
        real f = std::asinh(lp / R0) - std::asinh(lm / R0);
        (void)Rm; (void)Rp;
        real beta = std::atan2(P0 * lp, R0sq + aw * Rp) - std::atan2(P0 * lm, R0sq + aw * Rm);
        gpar += m * f; inv += P0 * f; bsum += beta;
    }
    inv -= aw * bsum;
    const real sgn = (w > 0) ? 1.0 : (w < 0 ? -1.0 : 0.0);
    Igrad = gpar + n * (sgn * bsum);
    Iinv = inv;
}

void triangle_integrals_linear(const Vec3& x, const std::array<Vec3, 3>& p, const Vec3& n, std::array<Vec3, 3>& Igrad,
                               std::array<real, 3>& Iinv) {
    const real w = dot(x - p[0], n);
    const Vec3 xs = x - n * w;                                          // Fusspunkt in der Ebene
    const real aw = std::abs(w);
    Vec3 gpar{0, 0, 0}, V1{0, 0, 0}; real inv = 0, bsum = 0;
    // M3 = inv (1 - n n^T) - sum_e m_e (x) (P0 f m_e + (R+ - R-) t_e): als Summe von Dyaden gespeichert
    std::array<Vec3, 3> me, ce;
    for (int e = 0; e < 3; ++e) {
        const Vec3& q1 = p[e]; const Vec3& q2 = p[(e + 1) % 3];
        Vec3 d = q2 - q1; real L = norm(d); Vec3 t = d / L; Vec3 m = cross(t, n);
        real lm = dot(q1 - xs, t), lp = dot(q2 - xs, t), P0 = dot(q1 - xs, m);
        real Rm = norm(x - q1), Rp = norm(x - q2), R0sq = P0 * P0 + w * w;
        const real R0 = std::max(std::sqrt(R0sq), 1e-13 * L);
        real f = std::asinh(lp / R0) - std::asinh(lm / R0);              // int_e dl / R
        real beta = std::atan2(P0 * lp, R0sq + aw * Rp) - std::atan2(P0 * lm, R0sq + aw * Rm);
        gpar += m * f; inv += P0 * f; bsum += beta;
        V1 += m * (0.5 * (lp * Rp - lm * Rm + R0sq * f));                 // int_e R dl = [l R + R0^2 asinh(l/R0)] / 2
        me[e] = m; ce[e] = m * (P0 * f) + t * (Rp - Rm);                   // int_e rho / R dl, rho = P0 m + l t
    }
    inv -= aw * bsum;
    const real sgn = (w > 0) ? 1.0 : (w < 0 ? -1.0 : 0.0);
    const Vec3 Ig0 = gpar + n * (sgn * bsum);                            // = triangle_integrals
    const Vec3 V3 = gpar * -1.0;
    // baryzentrische Koordinaten des Fusspunkts und Gradienten g_b = n x (p[b+2] - p[b+1]) / (2 A)
    const real A2 = dot(cross(p[1] - p[0], p[2] - p[0]), n);
    for (int b = 0; b < 3; ++b) {
        const Vec3& a = p[(b + 1) % 3]; const Vec3& c = p[(b + 2) % 3];
        const Vec3 g = cross(n, c - a) / A2;
        const real lam = dot(cross(c - a, xs - a), n) / A2;               // lambda_b(x*)
        Vec3 M3g = (g - n * dot(n, g)) * inv;
        for (int e = 0; e < 3; ++e) M3g = M3g - me[e] * dot(ce[e], g);
        Iinv[b] = lam * inv + dot(g, V1);
        Igrad[b] = Ig0 * lam + n * (w * dot(g, V3)) - M3g;
    }
}

}  // namespace cbem
