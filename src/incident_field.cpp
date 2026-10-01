#include "cbem/sources/incident_field.hpp"
#include <cmath>
#include <stdexcept>
#include "cbem/geometry/quadrature.hpp"
#include "cbem/sources/dipole.hpp"

namespace cbem {

namespace {
CVec3 cr(const CVec3& a, const CVec3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
void gauss_legendre_ab(int n, real a, real b, std::vector<real>& x, std::vector<real>& w) {   // Knoten auf [a, b]
    x.resize(n); w.resize(n);
    for (int i = 0; i < n; ++i) {
        real z = std::cos(pi * (i + 0.75) / (n + 0.5)), pp = 0;
        for (int it = 0; it < 100; ++it) {
            real p1 = 1, p2 = 0;
            for (int j = 1; j <= n; ++j) { const real p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
            pp = n * (z * p1 - p2) / (z * z - 1);
            const real dz = p1 / pp; z -= dz; if (std::abs(dz) < 1e-15) break;
        }
        x[i] = 0.5 * (a + b) + 0.5 * (b - a) * z; w[i] = (b - a) / ((1 - z * z) * pp * pp);
    }
}
}  // namespace

std::vector<cplx> IncidentField::project(const TriangleMesh& mesh, const Medium& m, int sub) const {
    MeshQuadrature q(mesh, QuadRule::subdivided(sub));
    const cplx se = std::sqrt(m.eps), sm = std::sqrt(m.mu);
    std::vector<cplx> h(8 * mesh.size(), cplx(0));
    CBEM_OMP(omp parallel for schedule(dynamic, 32))
    for (std::size_t t = 0; t < mesh.size(); ++t) {
        for (int a = 0; a < q.q; ++a) {
            CVec3 E, H; eval(q.points(t)[a], E, H);
            const Multivector F = Multivector::vector(E) * se + Multivector::blade(7) * Multivector::vector(H) * sm;
            const real w = q.weights(t)[a];
            for (int b = 0; b < 8; ++b) h[8 * t + b] += w * F.c[b];
        }
        for (int b = 0; b < 8; ++b) h[8 * t + b] /= std::sqrt(mesh.area[t]);
    }
    return h;
}

// --- ebene Welle ---------------------------------------------------------------------------------------------------------
PlaneWaveField::PlaneWaveField(const Medium& m, real omega, const Vec3& d, const CVec3& p) : m_(m), d_(d / norm(d)), p_(p) {
    k_ = plane_wave_incidence(m, omega, d_, p).k;                        // chiral: Helizitaetswelle mit k_sigma
    dxp_ = cr(CVec3{d_.x, d_.y, d_.z}, p);
    p2_ = std::norm(p[0]) + std::norm(p[1]) + std::norm(p[2]);
    C0_ = std::abs(std::real(std::sqrt(m.eps) / std::sqrt(m.mu))) * p2_;
}
void PlaneWaveField::eval(const Vec3& x, CVec3& E, CVec3& H) const {
    const cplx ph = std::exp(cplx(0, 1) * k_ * dot(d_, x)), r = std::sqrt(m_.eps) / std::sqrt(m_.mu);
    for (int a = 0; a < 3; ++a) { E[a] = p_[a] * ph; H[a] = r * dxp_[a] * ph; }
}

// --- Dipol ---------------------------------------------------------------------------------------------------------------
DipoleField::DipoleField(const Medium& m, real omega, const Vec3& r0, const CVec3& p, const CVec3& md) : m_(m), omega_(omega), r0_(r0), p_(p), md_(md) {
    if (std::abs(m.chi) > 0) throw std::invalid_argument("DipoleField: chirales Medium nicht unterstuetzt");
}
void DipoleField::eval(const Vec3& x, CVec3& E, CVec3& H) const { dipole_field(m_, omega_, r0_, p_, x, E, H, md_); }

// --- Strahl --------------------------------------------------------------------------------------------------------------
BeamField BeamField::focused(const Medium& m, real omega, const Vec3& focus, real NA, real f0, const CVec3& pp, int nt, int np) {
    BeamField B(m, omega, focus); B.k_ = m.k(omega);
    const real n = std::real(std::sqrt(m.eps * m.mu));
    if (!(NA > 0 && NA < n)) throw std::invalid_argument("BeamField: 0 < NA < n verlangt");
    const real tmax = std::asin(NA / n), smax = std::sin(tmax);
    std::vector<real> th, wth; gauss_legendre_ab(nt, 0.0, tmax, th, wth);
    for (int i = 0; i < nt; ++i) {
        const real st = std::sin(th[i]), ct = std::cos(th[i]);
        const real f = std::exp(-(st * st) / (f0 * f0 * smax * smax)) * std::sqrt(ct);
        for (int j = 0; j < np; ++j) {
            const real ph = 2 * pi * (j + 0.5) / np, cp = std::cos(ph), sp = std::sin(ph);
            const CVec3 phi{-sp, cp, 0.0}, the{ct * cp, ct * sp, -st};
            const cplx er = pp[0] * cp + pp[1] * sp, ep = -pp[0] * sp + pp[1] * cp;
            const real w = wth[i] * st * 2 * pi / np;                           // dOmega = sin theta dtheta dphi
            CVec3 a; for (int c = 0; c < 3; ++c) a[c] = w * f * (ep * phi[c] + er * the[c]);
            B.add(Vec3(st * cp, st * sp, ct), a, w);
        }
    }
    B.finish(); return B;
}

BeamField BeamField::gaussian(const Medium& m, real omega, const Vec3& focus, real w0, const CVec3& p, int nt, int np) {
    BeamField B(m, omega, focus); B.k_ = m.k(omega);
    const real kr = std::real(B.k_);
    // Winkelbereich: bis exp(-(k w0 sin theta/2)^2) < 1e-12, hoechstens pi/2
    const real smax = std::min(1.0, 2 * std::sqrt(12 * std::log(10.0)) / (kr * w0)), tmax = std::asin(smax);
    std::vector<real> th, wth; gauss_legendre_ab(nt, 0.0, tmax, th, wth);
    for (int i = 0; i < nt; ++i) {
        const real st = std::sin(th[i]), ct = std::cos(th[i]), g = std::exp(-std::pow(kr * w0 * st / 2, 2));
        for (int j = 0; j < np; ++j) {
            const real ph = 2 * pi * (j + 0.5) / np, cp = std::cos(ph), sp = std::sin(ph), w = wth[i] * st * 2 * pi / np;
            const CVec3 a{w * g * p[0] * ct, w * g * p[1] * ct, -w * g * st * (p[0] * cp + p[1] * sp)};
            B.add(Vec3(st * cp, st * sp, ct), a, w);
        }
    }
    B.finish(); return B;
}

void BeamField::add(const Vec3& dir, const CVec3& a, real w) {
    if (!(std::abs(m_.chi) > 0)) { dir_.push_back(dir); amp_.push_back(a); dw_.push_back(w); kc_.push_back(k_); return; }
    for (int s : {+1, -1}) {                                                    // Helizitaetsanteile mit eigener Wellenzahl
        CVec3 e = circular_polarization(dir, s); for (auto& c : e) c /= std::sqrt(2.0);
        const cplx c = std::conj(e[0]) * a[0] + std::conj(e[1]) * a[1] + std::conj(e[2]) * a[2];
        if (std::abs(c) == 0) continue;
        dir_.push_back(dir); amp_.push_back(CVec3{c * e[0], c * e[1], c * e[2]}); dw_.push_back(w);
        kc_.push_back(plane_wave_incidence(m_, omega_, dir, circular_polarization(dir, s)).k);
    }
}

void BeamField::finish() {
    // Leistung nach Parseval: P = 1/2 sqrt(eps/mu) (2 pi/k)^2 int |a|^2 dOmega = 1/2 sqrt(eps/mu) (2 pi/k)^2 sum |amp_i|^2 / dOmega_i
    const real z = std::real(std::sqrt(m_.eps / m_.mu));
    real S = 0;
    for (std::size_t i = 0; i < amp_.size(); ++i) { real a2 = 0; for (auto c : amp_[i]) a2 += std::norm(c); S += a2 / dw_[i] * std::pow(2 * pi / std::real(kc_[i]), 2); }
    P_ = 0.5 * z * S;
    const real s = 1 / std::sqrt(P_);
    for (auto& a : amp_) for (auto& c : a) c *= s;
    P_ = 1.0;
    CVec3 E, H; eval(focus_, E, H); I0_ = std::norm(E[0]) + std::norm(E[1]) + std::norm(E[2]);
}

real BeamField::poynting_flux(real dz, real L, int ng) const {
    real S = 0;
    CBEM_OMP(omp parallel for reduction(+ : S) schedule(dynamic))
    for (int i = 0; i < ng; ++i) for (int j = 0; j < ng; ++j) {
        const Vec3 x(focus_.x - L + 2 * L * (i + 0.5) / ng, focus_.y - L + 2 * L * (j + 0.5) / ng, focus_.z + dz);
        CVec3 E, H; eval(x, E, H); const CVec3 Hc{std::conj(H[0]), std::conj(H[1]), std::conj(H[2])};
        S += 0.5 * std::real(cr(E, Hc)[2]) * (2 * L / ng) * (2 * L / ng);
    }
    return S;
}

void BeamField::eval(const Vec3& x, CVec3& E, CVec3& H) const {
    E = CVec3{}; H = CVec3{};
    const Vec3 r = x - focus_; const cplx z = std::sqrt(m_.eps / m_.mu);
    for (std::size_t i = 0; i < dir_.size(); ++i) {
        const cplx ph = std::exp(cplx(0, 1) * kc_[i] * dot(dir_[i], r));
        const CVec3 kxa = cr(CVec3{dir_[i].x, dir_[i].y, dir_[i].z}, amp_[i]);
        for (int c = 0; c < 3; ++c) { E[c] += amp_[i][c] * ph; H[c] += z * kxa[c] * ph; }
    }
}

}  // namespace cbem
