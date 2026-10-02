#include "cbem/problems/curved_problem.hpp"

#include <cmath>
#include <stdexcept>

#include "cbem/solvers/gmres.hpp"
#include "cbem/sources/chiral_incidence.hpp"

namespace cbem {

std::vector<cplx> project_plane_wave_curved(const QuadraticMesh& m, cplx k, cplx eps, const Vec3& d0, const CVec3& p, int sub) {
    const Vec3 d = d0 / norm(d0);
    CurvedQuadrature Q(m, QuadRule::subdivided(sub));
    const auto S = curved_psi_matrices(m);
    const cplx se = std::sqrt(eps);
    const CVec3 dp = {d.y * p[2] - d.z * p[1], d.z * p[0] - d.x * p[2], d.x * p[1] - d.y * p[0]};
    const Multivector A = Multivector::vector(p) * se + Multivector::blade(7) * Multivector::vector(dp) * se;
    std::vector<cplx> h(24 * m.size(), cplx(0));
    CBEM_OMP(omp parallel for schedule(static))
    for (std::size_t t = 0; t < m.size(); ++t) {
        cplx c[3] = {0, 0, 0};
        for (int q = 0; q < Q.q; ++q) {
            const cplx ph = std::exp(cplx(0, 1) * k * dot(d, Q.points(t)[q])) * Q.weights(t)[q];
            for (int a = 0; a < 3; ++a) { real psi = 0; for (int kk = 0; kk < 3; ++kk) psi += S[t][a * 3 + kk] * Q.lam[q][kk]; c[a] += ph * psi; }
        }
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 8; ++b) h[8 * (3 * t + a) + b] = c[a] * A.c[b];
    }
    return h;
}

Multivector far_field_curved(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, const Vec3& xh, int sub) {
    CurvedQuadrature Q(m, QuadRule::subdivided(sub));
    const auto S = curved_psi_matrices(m);
    Multivector acc;                                                     // int e^{-ik xh.y} n(y) h(y) dS
    for (std::size_t t = 0; t < m.size(); ++t)
        for (int q = 0; q < Q.q; ++q) {
            Multivector h;
            for (int a = 0; a < 3; ++a) {
                real psi = 0; for (int kk = 0; kk < 3; ++kk) psi += S[t][a * 3 + kk] * Q.lam[q][kk];
                for (int b = 0; b < 8; ++b) h.c[b] += psi * hs[8 * (3 * t + a) + b];
            }
            const cplx ph = std::exp(cplx(0, -1) * k * dot(xh, Q.points(t)[q])) * Q.weights(t)[q];
            acc = acc + (Multivector::vector(Q.normals(t)[q]) * h) * ph;
        }
    const Multivector onepx = Multivector::blade(0) + Multivector::vector(xh);
    return (onepx * acc) * (cplx(0, -1) * k / (4 * pi));
}

real extinction_cross_section_curved(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
    const Multivector F = far_field_curved(m, hs, k, d / norm(d));
    const cplx pe = (std::conj(p[0]) * F.c[1] + std::conj(p[1]) * F.c[2] + std::conj(p[2]) * F.c[4]) / std::sqrt(eps);
    const real pn = std::norm(p[0]) + std::norm(p[1]) + std::norm(p[2]);
    return std::real(4 * pi / k * std::imag(pe)) / pn;
}

cplx forward_amplitude_curved(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
    const Multivector F = far_field_curved(m, hs, k, d / norm(d));
    const cplx pe = (std::conj(p[0]) * F.c[1] + std::conj(p[1]) * F.c[2] + std::conj(p[2]) * F.c[4]) / std::sqrt(eps);
    const real pn = std::norm(p[0]) + std::norm(p[1]) + std::norm(p[2]);
    return cplx(0, -1) * k * pe / pn;
}

CurvedScatteringProblem::CurvedScatteringProblem(const std::vector<QuadraticMesh>& bodies, const std::vector<Medium>& media, real omega,
                                                 Medium outer, HMatrixParams hp, EntryParams ep, CurvedNearParams np)
    : parts_(bodies), omega_(omega), outer_(outer) {
    if (media.size() != bodies.size() || bodies.empty()) throw std::invalid_argument("CurvedScatteringProblem: ein Medium je Koerper");
    { std::vector<const TriangleMesh*> ps; for (auto& b : parts_) ps.push_back(&b.flat); require_separated_all(ps, {}, "CurvedScatteringProblem"); }
    all_ = merge_quadratic(parts_, &begin_);
    auto add = [&](const QuadraticMesh& mesh, cplx k) -> CurvedCauchyOperator* {
        ents_.push_back(std::make_unique<CurvedKernelEntries>(mesh, k, ep, np));
        hms_.push_back(std::make_unique<CurvedHMatrix>(*ents_.back(), hp));
        cops_.push_back(std::make_unique<CurvedCauchyOperator>(mesh, *hms_.back()));
        return cops_.back().get();
    };
    if (std::abs(outer.chi) > 0) {                                      // chirales Aussenmedium (v0.59): P+ E_{k+} + P- E_{k-}
        CurvedCauchyOperator* op = add(all_, outer.k(omega, +1));
        CurvedCauchyOperator* om = add(all_, outer.k(omega, -1));
        chops_.push_back(std::make_unique<ChiralCauchyOperator>(*op, *om));
        outer_op_ = chops_.back().get();
    } else outer_op_ = add(all_, outer.k(omega));
    std::vector<const BoundaryOperator*> blocks;
    for (std::size_t b = 0; b < parts_.size(); ++b) {
        const Medium& md = media[b];
        if (std::abs(md.chi) > 0) {
            CurvedCauchyOperator* pp = add(parts_[b], md.k(omega, +1));
            CurvedCauchyOperator* pm = add(parts_[b], md.k(omega, -1));
            chops_.push_back(std::make_unique<ChiralCauchyOperator>(*pp, *pm));
            blocks.push_back(chops_.back().get());
        } else blocks.push_back(add(parts_[b], md.k(omega)));
    }
    std::vector<std::size_t> begin3; for (std::size_t v : begin_) begin3.push_back(3 * v);
    if (parts_.size() == 1) inner_ptr_ = blocks[0];
    else { inner_ = std::make_unique<BlockDiagonalOperator>(blocks, begin3); inner_ptr_ = inner_.get(); }
    std::vector<Medium> per_el(all_.size());
    for (std::size_t b = 0; b < parts_.size(); ++b) for (std::size_t t = begin_[b]; t < begin_[b + 1]; ++t) per_el[t] = media[b];
    T_ = std::make_unique<CurvedTransmissionOperator>(all_, curved_psi_matrices(all_), *inner_ptr_, *outer_op_, per_el, outer);
}

double CurvedScatteringProblem::hmatrix_bytes() const { double b = 0; for (auto& h : hms_) b += h->stats().bytes(); return b; }
double CurvedScatteringProblem::near_seconds() const { double s = 0; for (auto& e : ents_) s += e->near_seconds(); return s; }

PlaneWaveResult CurvedScatteringProblem::solve_rhs(const std::vector<cplx>& b, const SolveOptions& o) const {
    if (b.size() != T_->size()) throw std::invalid_argument("CurvedScatteringProblem::solve_rhs: rechte Seite der Laenge 24 N erwartet");
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->precondition(x, y); };
    PlaneWaveResult r; GmresResult g = gmres(A, b, r.h, &M, o.tol, o.restart, o.max_iter);
    r.iterations = g.iterations; r.residual = g.rel_residual;
    return r;
}

PlaneWaveResult CurvedScatteringProblem::solve_plane_wave(const Vec3& d0, const CVec3& p, const SolveOptions& o) const {
    const Vec3 d = d0 / norm(d0);
    // chirales Aussenmedium: Helizitaetswelle mit k_sigma, optisches Theorem im Kanal sigma (wie extinction_in_medium)
    const PlaneWaveIncidence inc = plane_wave_incidence(outer_, omega_, d, p);
    const auto b = project_plane_wave_curved(all_, inc.k, outer_.eps, d, p);
    PlaneWaveResult r = solve_rhs(b, o);
    std::vector<cplx> hs(r.h.size()); for (std::size_t i = 0; i < hs.size(); ++i) hs[i] = r.h[i] - b[i];
    const std::vector<cplx> hc = helicity_part(hs, inc.proj);
    r.sigma_ext = extinction_cross_section_curved(all_, hc, inc.k, outer_.eps, d, p);
    r.forward = forward_amplitude_curved(all_, hc, inc.k, outer_.eps, d, p);
    return r;
}

}  // namespace cbem
