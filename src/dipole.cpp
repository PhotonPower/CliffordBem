#include "cbem/sources/dipole.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include "cbem/geometry/quadrature.hpp"

namespace cbem {

namespace {
CVec3 ccross(const CVec3& a, const CVec3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
CVec3 cvec(const Vec3& v) { return {v.x, v.y, v.z}; }
// Gauss-Legendre-Knoten und -Gewichte auf [-1, 1] (Newton-Iteration)
void gauss_legendre(int n, std::vector<real>& x, std::vector<real>& w) {
    x.resize(n); w.resize(n);
    for (int i = 0; i < n; ++i) {
        real z = std::cos(pi * (i + 0.75) / (n + 0.5)), pp = 0;
        for (int it = 0; it < 100; ++it) {
            real p1 = 1, p2 = 0;
            for (int j = 1; j <= n; ++j) { const real p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
            pp = n * (z * p1 - p2) / (z * z - 1);
            const real dz = p1 / pp; z -= dz; if (std::abs(dz) < 1e-15) break;
        }
        x[i] = z; w[i] = 2 / ((1 - z * z) * pp * pp);
    }
}
}  // namespace

void dipole_field(const Medium& m, real omega, const Vec3& r0, const CVec3& p, const Vec3& x, CVec3& E, CVec3& H, const CVec3& md) {
    if (std::abs(m.chi) > 0) throw std::invalid_argument("dipole_field: chirales Medium nicht unterstuetzt");
    const cplx k = m.k(omega), ik = cplx(0, 1) * k;
    const Vec3 z = x - r0; const real r = norm(z); const Vec3 n = z / r;
    const cplx G = std::exp(ik * r) / (4 * pi * r);
    const CVec3 nc = cvec(n), nxp = ccross(nc, p), nxpxn = ccross(nxp, nc);
    const cplx np = nc[0] * p[0] + nc[1] * p[1] + nc[2] * p[2];
    const cplx f = (1.0 / (r * r) - ik / r);
    for (int a = 0; a < 3; ++a) {
        E[a] = (k * k * nxpxn[a] + (3.0 * nc[a] * np - p[a]) * f) * G / m.eps;
        H[a] = omega * k * nxp[a] * G * (1.0 - 1.0 / (ik * r));
    }
    if (md[0] != cplx(0) || md[1] != cplx(0) || md[2] != cplx(0)) {      // magnetischer Dipol (dual)
        const CVec3 nxm = ccross(nc, md), nxmxn = ccross(nxm, nc);
        const cplx nm = nc[0] * md[0] + nc[1] * md[1] + nc[2] * md[2];
        for (int a = 0; a < 3; ++a) {
            H[a] += (k * k * nxmxn[a] + (3.0 * nc[a] * nm - md[a]) * f) * G / m.mu;
            E[a] += -omega * k * nxm[a] * G * (1.0 - 1.0 / (ik * r));
        }
    }
}

std::vector<cplx> project_dipole(const TriangleMesh& mesh, const Medium& m, real omega, const Vec3& r0, const CVec3& p, int max_sub, const CVec3& md) {
    const cplx se = std::sqrt(m.eps), sm = std::sqrt(m.mu);
    std::map<int, QuadRule> rules;
    std::vector<cplx> h(8 * mesh.size(), cplx(0));
    std::vector<int> subs(mesh.size());
    for (std::size_t t = 0; t < mesh.size(); ++t) {                    // Unterteilung je Dreieck, Regeln vorab (seriell)
        const auto v = mesh.vertices(t);
        real rad = 0; for (const Vec3& q : v) rad = std::max(rad, norm(q - mesh.centroid[t]));
        const real dist = std::max(norm(r0 - mesh.centroid[t]) - rad, 1e-3 * mesh.hmax[t]);
        subs[t] = std::min(max_sub, std::max(2, static_cast<int>(std::ceil(4 * mesh.hmax[t] / dist))));
        if (!rules.count(subs[t])) rules.emplace(subs[t], QuadRule::subdivided(subs[t]));
    }
    CBEM_OMP(omp parallel for schedule(dynamic, 32))
    for (std::size_t t = 0; t < mesh.size(); ++t) {
        const auto v = mesh.vertices(t);
        const QuadRule& R = rules.at(subs[t]);
        for (std::size_t a = 0; a < R.w.size(); ++a) {
            const Vec3 x = v[0] * R.bary[a][0] + v[1] * R.bary[a][1] + v[2] * R.bary[a][2];
            CVec3 E, H; dipole_field(m, omega, r0, p, x, E, H, md);
            const Multivector F = Multivector::vector(E) * se + Multivector::blade(7) * Multivector::vector(H) * sm;
            const real w = R.w[a] * mesh.area[t];
            for (int b = 0; b < 8; ++b) h[8 * t + b] += w * F.c[b];
        }
        for (int b = 0; b < 8; ++b) h[8 * t + b] /= std::sqrt(mesh.area[t]);
    }
    return h;
}

int helicity_projector_sign(int s) {
    // F = p + I d x p der zirkularen Welle p = circular_polarization(d, s) liegt ganz im Bild von P_+ oder P_-
    const Vec3 d(0, 0, 1); const CVec3 p = circular_polarization(d, s);
    const CVec3 dxp{d.y * p[2] - d.z * p[1], d.z * p[0] - d.x * p[2], d.x * p[1] - d.y * p[0]};
    const Multivector F = Multivector::vector(p) + Multivector::blade(7) * Multivector::vector(dxp);
    real np = 0, nm = 0;
    const Multivector A = (Multivector::blade(0, 0.5) + Multivector::blade(7, cplx(0, 0.5))) * F, B = (Multivector::blade(0, 0.5) + Multivector::blade(7, cplx(0, -0.5))) * F;
    for (int c = 0; c < 8; ++c) { np += std::norm(A.c[c]); nm += std::norm(B.c[c]); }
    return np > nm ? +1 : -1;
}

DipoleRates dipole_rates(const TriangleMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b, const Medium& m,
                         real omega, const Vec3& r0, const CVec3& p, int ntheta, const CVec3& md) {
    if (std::abs(m.chi) > 0) throw std::invalid_argument("dipole_rates: chirales Aussenmedium nicht unterstuetzt");
    DipoleRates R;
    const cplx k = m.k(omega), se = std::sqrt(m.eps);
    std::vector<cplx> hs(h.size()); for (std::size_t i = 0; i < h.size(); ++i) hs[i] = h[i] - b[i];
    real hm = 0; for (real x : outer.hmax) hm += x; hm /= std::max<std::size_t>(1, outer.hmax.size());
    R.distance = distance_to_surface(outer, {r0}, 10 * hm)[0];
    R.too_close = R.distance < 0.02 * hm;
    // Gesamtrate aus dem Streufeld am Dipolort
    const Multivector Fs = scattered_field(outer, hs, k, {r0})[0];
    const CVec3 Es{Fs.c[1] / se, Fs.c[2] / se, Fs.c[4] / se};
    const cplx sm = std::sqrt(m.mu);
    const int BIV[3] = {6, 5, 3}; const real BS[3] = {1, -1, 1};
    CVec3 Hs; for (int a = 0; a < 3; ++a) Hs[a] = Fs.c[BIV[a]] * BS[a] / sm;
    const real p2 = std::norm(p[0]) + std::norm(p[1]) + std::norm(p[2]), m2 = std::norm(md[0]) + std::norm(md[1]) + std::norm(md[2]);
    const cplx pe = std::conj(p[0]) * Es[0] + std::conj(p[1]) * Es[1] + std::conj(p[2]) * Es[2];
    const cplx mh = std::conj(md[0]) * Hs[0] + std::conj(md[1]) * Hs[1] + std::conj(md[2]) * Hs[2];
    const real k3 = std::real(k * k * k);
    R.total = 1 + (std::imag(pe) + std::imag(mh)) / (k3 * p2 / (6 * pi * std::real(m.eps)) + k3 * m2 / (6 * pi * std::real(m.mu)));
    // Strahlende Rate: Fernfeld von Dipol + Streufeld ueber die Einheitskugel
    std::vector<real> ct, wt; gauss_legendre(ntheta, ct, wt);
    const int nphi = 2 * ntheta;
    real P = 0, P0 = 0, Pp = 0, Pm = 0, P0p = 0, P0m = 0;
    const int sp = helicity_projector_sign(+1);                         // Projektor der Helizitaet s = +1
    const Multivector Pplus = Multivector::blade(0, 0.5) + Multivector::blade(7, cplx(0, 0.5 * sp));
    const Multivector Pminus = Multivector::blade(0, 0.5) + Multivector::blade(7, cplx(0, -0.5 * sp));
    auto split = [&](const CVec3& E, const Vec3& xh, real& a, real& b2) {  // Helizitaetsanteile des transversalen Feldes E
        const CVec3 xc = cvec(xh), xe = ccross(xc, E);
        const Multivector F = Multivector::vector(E) + Multivector::blade(7) * Multivector::vector(xe);
        const Multivector A = Pplus * F, B = Pminus * F; a = 0; b2 = 0;
        for (int c = 0; c < 8; ++c) { a += std::norm(A.c[c]); b2 += std::norm(B.c[c]); }
    };
    CBEM_OMP(omp parallel for schedule(dynamic) reduction(+ : P, P0, Pp, Pm, P0p, P0m))
    for (int i = 0; i < ntheta; ++i) {
        const real st = std::sqrt(std::max(0.0, 1 - ct[i] * ct[i]));
        for (int j = 0; j < nphi; ++j) {
            const real ph = 2 * pi * (j + 0.5) / nphi, w = wt[i] * 2 * pi / nphi;
            const Vec3 xh(st * std::cos(ph), st * std::sin(ph), ct[i]); const CVec3 xc = cvec(xh);
            const CVec3 e0 = ccross(ccross(xc, p), xc), em = ccross(xc, md);
            const cplx ex = std::exp(-cplx(0, 1) * k * dot(xh, r0));
            const cplx ph0 = k * k / (4 * pi * m.eps) * ex, phm = -omega * k / (4 * pi) * ex;
            const Multivector Ff = far_field(outer, hs, k, xh);
            real a = 0, a0 = 0; CVec3 E0v, Etv;
            for (int c = 0; c < 3; ++c) {
                const int VEC[3] = {1, 2, 4};
                const cplx E0 = ph0 * e0[c] + phm * em[c], Et = E0 + Ff.c[VEC[c]] / se;
                E0v[c] = E0; Etv[c] = Et;
                a += std::norm(Et); a0 += std::norm(E0);
            }
            P += w * a; P0 += w * a0;
            real hp, hm2, hp0, hm0; split(Etv, xh, hp, hm2); split(E0v, xh, hp0, hm0);
            Pp += w * hp; Pm += w * hm2; P0p += w * hp0; P0m += w * hm0;
        }
    }
    R.radiative = P / P0;
    R.nonradiative = R.total - R.radiative;
    // Aufteilung nach Helizitaet: |P+ F|^2 + |P- F|^2 = |F|^2 = 2 |E|^2 fuer transversale Felder
    R.rad_plus = R.radiative * Pp / (Pp + Pm); R.rad_minus = R.radiative * Pm / (Pp + Pm);
    R.glum = 2 * (Pp - Pm) / (Pp + Pm); R.glum_free = 2 * (P0p - P0m) / (P0p + P0m);
    return R;
}

real fluorescence_enhancement(const real exc[3], const DipoleRates rates[3], real q0) {
    real F = 0;
    for (int a = 0; a < 3; ++a) F += exc[a] * rates[a].radiative / (rates[a].total + (1 - q0) / q0);
    return F / q0;
}

}  // namespace cbem
