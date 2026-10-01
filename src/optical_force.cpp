#include "cbem/sources/optical_force.hpp"
#include <cmath>
#include <stdexcept>
#include "cbem/geometry/quadrature.hpp"

namespace cbem {

namespace {
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
Vec3 integrate(const std::vector<NearFieldPoint>& f, const std::vector<Vec3>& n, const std::vector<real>& w, const Medium& m) {
    Vec3 F(0, 0, 0);
    for (std::size_t i = 0; i < f.size(); ++i) F = F + stress_dot_normal(f[i].E, f[i].H, n[i], m) * w[i];
    return F;
}
}  // namespace

Vec3 stress_dot_normal(const CVec3& E, const CVec3& H, const Vec3& n, const Medium& m) {
    if (std::abs(m.chi) > 0) throw std::invalid_argument("Spannungstensor: chirales Medium nicht unterstuetzt");
    const cplx En = E[0] * n.x + E[1] * n.y + E[2] * n.z, Hn = H[0] * n.x + H[1] * n.y + H[2] * n.z;
    real e2 = 0, h2 = 0; for (int a = 0; a < 3; ++a) { e2 += std::norm(E[a]); h2 += std::norm(H[a]); }
    const real eps = std::real(m.eps), mu = std::real(m.mu), u = 0.5 * (eps * e2 + mu * h2);
    const real nn[3] = {n.x, n.y, n.z}; real t[3];
    for (int a = 0; a < 3; ++a) t[a] = 0.5 * (std::real(eps * E[a] * std::conj(En) + mu * H[a] * std::conj(Hn)) - u * nn[a]);
    return Vec3(t[0], t[1], t[2]);
}

std::vector<Vec3> force_from_traces(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, const std::vector<std::size_t>& body_begin) {
    const cplx se = std::sqrt(m.eps), sm = std::sqrt(m.mu);
    const int VEC[3] = {1, 2, 4}, BIV[3] = {6, 5, 3}; const real BS[3] = {1, -1, 1};
    std::vector<Vec3> F(body_begin.size() - 1, Vec3(0, 0, 0));
    for (std::size_t b = 0; b + 1 < body_begin.size(); ++b)
        for (std::size_t t = body_begin[b]; t < body_begin[b + 1]; ++t) {
            const real s = 1 / std::sqrt(outer.area[t]); CVec3 E, H;
            for (int a = 0; a < 3; ++a) { E[a] = h[8 * t + VEC[a]] * s / se; H[a] = h[8 * t + BIV[a]] * BS[a] * s / sm; }
            F[b] = F[b] + stress_dot_normal(E, H, outer.normal[t], m) * outer.area[t];
        }
    return F;
}

Vec3 force_on_sphere(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d, const CVec3& p,
                     const Vec3& c, real R, int ntheta, const NearFieldOptions& o) {
    std::vector<real> ct, wt; gauss_legendre(ntheta, ct, wt);
    const int nphi = 2 * ntheta;
    std::vector<Vec3> pts, nrm; std::vector<real> w;
    for (int i = 0; i < ntheta; ++i) {
        const real st = std::sqrt(std::max(0.0, 1 - ct[i] * ct[i]));
        for (int j = 0; j < nphi; ++j) {
            const real ph = 2 * pi * (j + 0.5) / nphi; const Vec3 nv(st * std::cos(ph), st * std::sin(ph), ct[i]);
            pts.push_back(c + nv * R); nrm.push_back(nv); w.push_back(wt[i] * 2 * pi / nphi * R * R);
        }
    }
    const auto f = exterior_near_field(outer, h, m, omega, d, p, pts, o);
    for (const auto& q : f) if (q.inside || q.too_close) throw std::invalid_argument("force_on_sphere: Kugel schneidet einen Koerper");
    return integrate(f, nrm, w, m);
}

Vec3 force_on_offset(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d, const CVec3& p,
                     const TriangleMesh& body, real delta, const NearFieldOptions& o) {
    TriangleMesh S = offset_surface(body, delta); S.compute_geometry();
    MeshQuadrature q(S, QuadRule::dunavant7());
    std::vector<Vec3> pts, nrm; std::vector<real> w;
    for (std::size_t t = 0; t < S.size(); ++t) {
        const Vec3* qp = q.points(t); const real* qw = q.weights(t);
        for (int a = 0; a < q.q; ++a) { pts.push_back(qp[a]); nrm.push_back(S.normal[t]); w.push_back(qw[a]); }
    }
    const auto f = exterior_near_field(outer, h, m, omega, d, p, pts, o);
    for (const auto& x : f) if (x.inside) throw std::invalid_argument("force_on_offset: Parallelflaeche schneidet einen Koerper");
    return integrate(f, nrm, w, m);
}

std::vector<FieldGradient> fields_with_gradients(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega,
                                                 const Vec3& d, const CVec3& p, const std::vector<Vec3>& x, real delta, const NearFieldOptions& o) {
    std::vector<Vec3> pts;
    const Vec3 ex[3] = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    for (const Vec3& q : x) { pts.push_back(q); for (int i = 0; i < 3; ++i) { pts.push_back(q + ex[i] * delta); pts.push_back(q - ex[i] * delta); } }
    const auto f = exterior_near_field(outer, h, m, omega, d, p, pts, o);
    std::vector<FieldGradient> g(x.size());
    for (std::size_t k = 0; k < x.size(); ++k) {
        const auto& c = f[7 * k];
        if (c.inside) throw std::invalid_argument("fields_with_gradients: Punkt innerhalb eines Koerpers");
        g[k].E = c.E; g[k].H = c.H;
        for (int i = 0; i < 3; ++i) {
            const auto& a = f[7 * k + 1 + 2 * i]; const auto& b = f[7 * k + 2 + 2 * i];
            for (int j = 0; j < 3; ++j) { g[k].dE[i][j] = (a.E[j] - b.E[j]) / (2 * delta); g[k].dH[i][j] = (a.H[j] - b.H[j]) / (2 * delta); }
        }
    }
    return g;
}

Vec3 dipole_particle_force(const FieldGradient& g, const DipolePolarizability& a, const Medium& m, real omega) {
    if (std::abs(m.chi) > 0) throw std::invalid_argument("dipole_particle_force: chirales Medium nicht unterstuetzt");
    const cplx k = m.k(omega), se = std::sqrt(m.eps), sm = std::sqrt(m.mu), I(0, 1);
    CVec3 pv, mv;
    for (int j = 0; j < 3; ++j) {
        const cplx e = se * g.E[j], hh = sm * g.H[j];
        pv[j] = se * (a.Ae * e + I * a.Ac * hh);
        mv[j] = sm * (a.Am * hh - I * a.Ac * e);
    }
    real F[3];
    for (int i = 0; i < 3; ++i) {
        cplx s = 0; for (int j = 0; j < 3; ++j) s += pv[j] * std::conj(g.dE[i][j]) + mv[j] * std::conj(g.dH[i][j]);
        F[i] = 0.5 * std::real(s);
    }
    const CVec3 pm{pv[1] * std::conj(mv[2]) - pv[2] * std::conj(mv[1]), pv[2] * std::conj(mv[0]) - pv[0] * std::conj(mv[2]), pv[0] * std::conj(mv[1]) - pv[1] * std::conj(mv[0])};
    const real c = std::real(omega * k * k * k) / (12 * pi);
    for (int i = 0; i < 3; ++i) F[i] -= c * std::real(pm[i]);
    return Vec3(F[0], F[1], F[2]);
}

}  // namespace cbem
