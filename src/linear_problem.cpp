#include "cbem/problems/linear_problem.hpp"

#include <cmath>
#include <stdexcept>

#include "cbem/geometry/quadrature.hpp"
#include "cbem/solvers/gmres.hpp"
#include "cbem/sources/chiral_incidence.hpp"

namespace cbem {

namespace {

// lambda_k an einem Punkt des Dreiecks t (baryzentrisch)
std::array<real, 3> barycentric(const TriangleMesh& m, std::size_t t, const Vec3& x) {
    const auto v = m.vertices(t);
    const Vec3& n = m.normal[t];
    const real A2 = dot(cross(v[1] - v[0], v[2] - v[0]), n);
    std::array<real, 3> l{};
    for (int k = 0; k < 3; ++k) {
        const Vec3& a = v[(k + 1) % 3]; const Vec3& c = v[(k + 2) % 3];
        l[k] = dot(cross(c - a, x - a), n) / A2;
    }
    return l;
}

}  // namespace

std::vector<cplx> project_plane_wave_linear(const TriangleMesh& m, cplx k, cplx eps, const Vec3& d0, const CVec3& p, int sub) {
    const Vec3 d = d0 / norm(d0);
    const QuadRule R = QuadRule::subdivided(sub);
    MeshQuadrature q(m, R);
    const cplx se = std::sqrt(eps);
    const CVec3 dp = {d.y * p[2] - d.z * p[1], d.z * p[0] - d.x * p[2], d.x * p[1] - d.y * p[0]};
    const Multivector A = Multivector::vector(p) * se + Multivector::blade(7) * Multivector::vector(dp) * se;
    std::vector<cplx> h(24 * m.size(), cplx(0));
    CBEM_OMP(omp parallel for schedule(static))
    for (std::size_t t = 0; t < m.size(); ++t) {
        cplx Lk[3] = {0, 0, 0};                                          // int lambda_k e^{ik d.y}
        for (int a = 0; a < q.q; ++a) {
            const cplx ph = std::exp(cplx(0, 1) * k * dot(d, q.points(t)[a])) * q.weights(t)[a];
            for (int kk = 0; kk < 3; ++kk) Lk[kk] += ph * R.bary[a][kk];
        }
        const auto S = psi_matrix(m.area[t]);
        for (int a = 0; a < 3; ++a) {
            cplx c = 0; for (int kk = 0; kk < 3; ++kk) c += S[a * 3 + kk] * Lk[kk];
            for (int b = 0; b < 8; ++b) h[8 * (3 * t + a) + b] = c * A.c[b];
        }
    }
    return h;
}

Multivector linear_trace_value(const TriangleMesh& m, const std::vector<cplx>& h, std::size_t t, const Vec3& x) {
    const auto l = barycentric(m, t, x);
    const auto S = psi_matrix(m.area[t]);
    Multivector v;
    for (int a = 0; a < 3; ++a) {
        real psi = 0; for (int kk = 0; kk < 3; ++kk) psi += S[a * 3 + kk] * l[kk];
        for (int b = 0; b < 8; ++b) v.c[b] += psi * h[8 * (3 * t + a) + b];
    }
    return v;
}

std::vector<cplx> linear_to_constant(const TriangleMesh& m, const std::vector<cplx>& h) {
    // Mittelwert: int psi_a = sqrt(A/3); konstante Koeffizienten = int h / sqrt(A) = sum_a h_a / sqrt(3)
    std::vector<cplx> c(8 * m.size(), cplx(0));
    for (std::size_t t = 0; t < m.size(); ++t)
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 8; ++b) c[8 * t + b] += h[8 * (3 * t + a) + b] / std::sqrt(3.0);
    return c;
}

Multivector far_field_linear(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, const Vec3& xh, int sub) {
    const QuadRule R = QuadRule::subdivided(sub);
    MeshQuadrature q(m, R);
    const Multivector onepx = Multivector::blade(0) + Multivector::vector(xh);
    Multivector F;
    for (std::size_t t = 0; t < m.size(); ++t) {
        const auto S = psi_matrix(m.area[t]);
        cplx Lk[3] = {0, 0, 0};                                          // int lambda_k e^{-ik xh.y}
        for (int a = 0; a < q.q; ++a) {
            const cplx ph = std::exp(cplx(0, -1) * k * dot(xh, q.points(t)[a])) * q.weights(t)[a];
            for (int kk = 0; kk < 3; ++kk) Lk[kk] += ph * R.bary[a][kk];
        }
        Multivector hm;                                                  // int h(y) e^{-ik xh.y} dS
        for (int a = 0; a < 3; ++a) {
            cplx c = 0; for (int kk = 0; kk < 3; ++kk) c += S[a * 3 + kk] * Lk[kk];
            for (int b = 0; b < 8; ++b) hm.c[b] += c * hs[8 * (3 * t + a) + b];
        }
        F = F + onepx * (Multivector::vector(m.normal[t]) * hm);
    }
    return F * (cplx(0, -1) * k / (4 * pi));
}

real extinction_cross_section_linear(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
    const Multivector F = far_field_linear(m, hs, k, d / norm(d));
    const cplx se = std::sqrt(eps);
    const cplx pe = (std::conj(p[0]) * F.c[1] + std::conj(p[1]) * F.c[2] + std::conj(p[2]) * F.c[4]) / se;
    const real pn = std::norm(p[0]) + std::norm(p[1]) + std::norm(p[2]);
    return std::real(4 * pi / k * std::imag(pe)) / pn;
}

cplx forward_amplitude_linear(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
    const Multivector F = far_field_linear(m, hs, k, d / norm(d));
    const cplx se = std::sqrt(eps);
    const cplx pe = (std::conj(p[0]) * F.c[1] + std::conj(p[1]) * F.c[2] + std::conj(p[2]) * F.c[4]) / se;
    const real pn = std::norm(p[0]) + std::norm(p[1]) + std::norm(p[2]);
    return cplx(0, -1) * k * pe / pn;
}

LinearScatteringProblem::LinearScatteringProblem(const std::vector<TriangleMesh>& bodies, const std::vector<Medium>& media, real omega,
                                                 Medium outer, HMatrixParams hp, EntryParams ep)
    : mb_(make_multibody(bodies)), omega_(omega), outer_(outer) {
    if (media.size() != bodies.size()) throw std::invalid_argument("LinearScatteringProblem: ein Medium je Koerper");
    { std::vector<const TriangleMesh*> ps; for (auto& b : bodies) ps.push_back(&b); require_separated_all(ps, {}, "LinearScatteringProblem"); }
    const TriangleMesh& m = mb_.all;
    auto add = [&](const TriangleMesh& mesh, cplx k) -> LinearCauchyOperator* {
        ents_.push_back(std::make_unique<LinearKernelEntries>(mesh, k, ep));
        hms_.push_back(std::make_unique<KernelHMatrix>(*ents_.back(), hp));
        cops_.push_back(std::make_unique<LinearCauchyOperator>(mesh, *hms_.back()));
        return cops_.back().get();
    };
    if (std::abs(outer.chi) > 0) {                                      // chirales Aussenmedium (v0.59): P+ E_{k+} + P- E_{k-}
        LinearCauchyOperator* op = add(m, outer.k(omega, +1));
        LinearCauchyOperator* om = add(m, outer.k(omega, -1));
        chops_.push_back(std::make_unique<ChiralCauchyOperator>(*op, *om));
        outer_op_ = chops_.back().get();
    } else outer_op_ = add(m, outer.k(omega));
    std::vector<const BoundaryOperator*> blocks;
    std::vector<std::size_t> begin3;                                     // Basisfunktionen je Koerper: 3 x Dreiecke
    for (std::size_t b = 0; b < bodies.size(); ++b) {
        const Medium& md = media[b]; const TriangleMesh& mp = mb_.parts[b];
        if (std::abs(md.chi) > 0) {
            LinearCauchyOperator* pp = add(mp, md.k(omega, +1));
            LinearCauchyOperator* pm = add(mp, md.k(omega, -1));
            chops_.push_back(std::make_unique<ChiralCauchyOperator>(*pp, *pm));
            blocks.push_back(chops_.back().get());
        } else blocks.push_back(add(mp, md.k(omega)));
    }
    for (std::size_t v : mb_.body_begin) begin3.push_back(3 * v);
    if (bodies.size() == 1) inner_ptr_ = blocks[0];
    else { inner_ = std::make_unique<BlockDiagonalOperator>(blocks, begin3); inner_ptr_ = inner_.get(); }
    std::vector<Vec3> normals(3 * m.size());
    std::vector<Medium> per_index(3 * m.size());
    for (std::size_t b = 0; b < bodies.size(); ++b)
        for (std::size_t t = mb_.body_begin[b]; t < mb_.body_begin[b + 1]; ++t)
            for (int a = 0; a < 3; ++a) { normals[3 * t + a] = m.normal[t]; per_index[3 * t + a] = media[b]; }
    T_ = std::make_unique<TransmissionOperator>(normals, *inner_ptr_, *outer_op_, per_index, outer);
}

double LinearScatteringProblem::hmatrix_bytes() const { double b = 0; for (auto& h : hms_) b += h->stats().bytes(); return b; }

std::size_t LinearScatteringProblem::near_pairs() const { std::size_t n = 0; for (auto& e : ents_) n += e->near_pairs(); return n; }

PlaneWaveResult LinearScatteringProblem::solve_rhs(const std::vector<cplx>& b, const SolveOptions& o) const {
    if (b.size() != T_->size()) throw std::invalid_argument("LinearScatteringProblem::solve_rhs: rechte Seite der Laenge 24 N erwartet");
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->precondition(x, y); };
    PlaneWaveResult r; GmresResult g = gmres(A, b, r.h, &M, o.tol, o.restart, o.max_iter);
    r.iterations = g.iterations; r.residual = g.rel_residual;
    return r;
}

PlaneWaveResult LinearScatteringProblem::solve_plane_wave(const Vec3& d0, const CVec3& p, const SolveOptions& o) const {
    const Vec3 d = d0 / norm(d0);
    const TriangleMesh& m = mb_.all;
    // chirales Aussenmedium: Helizitaetswelle mit k_sigma, optisches Theorem im Kanal sigma (wie extinction_in_medium)
    const PlaneWaveIncidence inc = plane_wave_incidence(outer_, omega_, d, p);
    const auto b = project_plane_wave_linear(m, inc.k, outer_.eps, d, p);
    PlaneWaveResult r = solve_rhs(b, o);
    std::vector<cplx> hs(r.h.size()); for (std::size_t i = 0; i < hs.size(); ++i) hs[i] = r.h[i] - b[i];
    const std::vector<cplx> hc = helicity_part(hs, inc.proj);
    r.sigma_ext = extinction_cross_section_linear(m, hc, inc.k, outer_.eps, d, p);
    r.forward = forward_amplitude_linear(m, hc, inc.k, outer_.eps, d, p);
    return r;
}

}  // namespace cbem
