#include "cbem/problems/scattering_problem.hpp"
#include <stdexcept>
#include "cbem/sources/fields.hpp"
#include "cbem/operators/dense_blocks.hpp"

namespace cbem {

ScatteringProblem::ScatteringProblem(const std::vector<TriangleMesh>& bodies, const std::vector<Medium>& media, real omega,
                                     Medium outer, HMatrixParams hp, EntryParams ep, bool union_interior)
    : mb_(make_multibody(bodies)), omega_(omega), outer_(outer) {
    if (media.size() != bodies.size()) throw std::invalid_argument("ScatteringProblem: ein Medium je Koerper");
    if (std::abs(outer.chi) > 0) throw std::invalid_argument("ScatteringProblem: chirales Aussenmedium nicht unterstuetzt");
    const TriangleMesh& m = mb_.all;
    auto add = [&](const TriangleMesh& mesh, cplx k) -> CauchyOperator* {
        ents_.push_back(std::make_unique<KernelEntries>(mesh, k, ep));
        hms_.push_back(std::make_unique<KernelHMatrix>(*ents_.back(), hp));
        cops_.push_back(std::make_unique<CauchyOperator>(mesh, *hms_.back()));
        return cops_.back().get();
    };
    outer_op_ = add(m, outer.k(omega)); outer_entries_ = ents_.back().get();
    std::vector<Medium> per_tri(m.size());
    for (std::size_t b = 0; b < bodies.size(); ++b)
        for (std::size_t t = mb_.body_begin[b]; t < mb_.body_begin[b + 1]; ++t) per_tri[t] = media[b];
    if (union_interior) {                        // ein Innenoperator auf der Vereinigung (gleiche achirale Medien)
        for (const auto& md : media)
            if (std::abs(md.chi) > 0 || md.eps != media[0].eps || md.mu != media[0].mu) throw std::invalid_argument("union_interior: gleiche achirale Medien noetig");
        inner_ptr_ = add(m, media[0].k(omega));
        inner_parts_ = {{{ents_.back().get(), 0}}}; inner_begin_ = {0, m.size()};
    } else {
        inner_begin_ = mb_.body_begin;
        std::vector<const BoundaryOperator*> blocks;
        for (std::size_t b = 0; b < bodies.size(); ++b) {
            const Medium& md = media[b]; const TriangleMesh& mp = mb_.parts[b];
            if (std::abs(md.chi) > 0) {
                CauchyOperator* p = add(mp, md.k(omega, +1)); const KernelEntries* ep = ents_.back().get();
                CauchyOperator* q = add(mp, md.k(omega, -1)); const KernelEntries* eq = ents_.back().get();
                chops_.push_back(std::make_unique<ChiralCauchyOperator>(*p, *q)); blocks.push_back(chops_.back().get());
                inner_parts_.push_back({{ep, +1}, {eq, -1}});
            } else { blocks.push_back(add(mp, md.k(omega))); inner_parts_.push_back({{ents_.back().get(), 0}}); }
        }
        if (bodies.size() == 1) inner_ptr_ = blocks[0];
        else { inner_ = std::make_unique<BlockDiagonalOperator>(blocks, mb_.body_begin); inner_ptr_ = inner_.get(); }
    }
    T_ = std::make_unique<TransmissionOperator>(m, *inner_ptr_, *outer_op_, per_tri, outer);
    // Hilfsgroessen fuer Einzeleintraege
    Lcn_.resize(m.size()); body_of_.assign(m.size(), 0);
    const Multivector e[4] = {Multivector::blade(0), Multivector::blade(1), Multivector::blade(2), Multivector::blade(4)};
    for (std::size_t t = 0; t < m.size(); ++t) { Multivector n = Multivector::vector(m.normal[t]); for (int c = 0; c < 4; ++c) Lcn_[t][c] = (e[c] * n).left_matrix(); }
    for (std::size_t b = 0; b + 1 < inner_begin_.size(); ++b) for (std::size_t t = inner_begin_[b]; t < inner_begin_[b + 1]; ++t) body_of_[t] = b;
}

Mat8 ScatteringProblem::system_entry(std::size_t i, std::size_t j) const {
    const TriangleMesh& m = mb_.all; const real s = -2.0 / std::sqrt(m.area[i] * m.area[j]);
    auto cauchy = [&](const KernelEntries& E, std::size_t a, std::size_t b) {
        KernelComp K = E.exact(a, b); Mat8 R{};
        for (int c = 0; c < 4; ++c) for (int q = 0; q < 64; ++q) R[q] += K[c] * Lcn_[j][c][q];
        for (auto& v : R) v *= s;
        return R;
    };
    Mat8 E2 = cauchy(*outer_entries_, i, j), E1{};
    const std::size_t b = body_of_[i];
    if (body_of_[j] == b) {
        const std::size_t li = i - inner_begin_[b], lj = j - inner_begin_[b];
        for (const auto& part : inner_parts_[b]) {
            Mat8 C = cauchy(*part.E, li, lj);
            if (part.helicity != 0) { Mat8 P = helicity_projector(part.helicity), PC{}; for (int r = 0; r < 8; ++r) for (int q = 0; q < 8; ++q) { cplx v = 0; for (int t = 0; t < 8; ++t) v += P[r * 8 + t] * C[t * 8 + q]; PC[r * 8 + q] = v; } C = PC; }
            for (int q = 0; q < 64; ++q) E1[q] += C[q];
        }
    }
    const Mat8& Jj = T_->J()[j]; Mat8 T{};
    for (int r = 0; r < 8; ++r) for (int q = 0; q < 8; ++q) {
        cplx v = 0; for (int t = 0; t < 8; ++t) v += E1[r * 8 + t] * Jj[t * 8 + q];
        T[r * 8 + q] = 0.5 * E2[r * 8 + q] - 0.5 * v;
    }
    if (i == j) for (int r = 0; r < 8; ++r) for (int q = 0; q < 8; ++q) T[r * 8 + q] += 0.5 * ((r == q ? 1.0 : 0.0) + Jj[r * 8 + q]);
    return T;
}

void ScatteringProblem::use_hodlr_preconditioner(HodlrParams p) {
    prec_.reset();
    hodlr_ = std::make_unique<HodlrSolver>(mb_.all, [this](std::size_t i, std::size_t j) { return system_entry(i, j); }, p);
}

Matrix ScatteringProblem::inner_block(const std::vector<std::size_t>& B) const {
    const std::size_t nb = B.size(); Matrix A(8 * nb, 8 * nb);
    const Mat8 Pp = helicity_projector(+1), Pm = helicity_projector(-1);
    for (std::size_t body = 0; body + 1 < inner_begin_.size(); ++body) {
        std::vector<std::size_t> pos, loc;                       // Positionen in B, lokale Indizes im Koerper
        for (std::size_t a = 0; a < nb; ++a)
            if (B[a] >= inner_begin_[body] && B[a] < inner_begin_[body + 1]) { pos.push_back(a); loc.push_back(B[a] - inner_begin_[body]); }
        if (pos.empty()) continue;
        for (const auto& part : inner_parts_[body]) {
            Matrix C = cauchy_block(*part.E, loc, loc);
            for (std::size_t a = 0; a < pos.size(); ++a)
                for (std::size_t b = 0; b < pos.size(); ++b)
                    for (int r = 0; r < 8; ++r)
                        for (int q = 0; q < 8; ++q) {
                            cplx v;
                            if (part.helicity == 0) v = C(8 * a + r, 8 * b + q);
                            else { const Mat8& P = part.helicity > 0 ? Pp : Pm; v = 0; for (int t = 0; t < 8; ++t) v += P[r * 8 + t] * C(8 * a + t, 8 * b + q); }
                            A(8 * pos[a] + r, 8 * pos[b] + q) += v;
                        }
        }
    }
    return A;
}

void ScatteringProblem::use_block_preconditioner(const std::vector<std::vector<std::size_t>>& groups) {
    hodlr_.reset();
    prec_ = std::make_unique<BlockPreconditioner>(mb_.all, [this](const std::vector<std::size_t>& B) { return inner_block(B); },
                                                  *outer_entries_, *T_, groups);
}

double ScatteringProblem::hmatrix_bytes() const { double b = 0; for (auto& h : hms_) b += h->stats().bytes(); return b; }

PlaneWaveResult ScatteringProblem::solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o) const {
    const TriangleMesh& m = mb_.all; const cplx k2 = outer_.k(omega_);
    auto b = project_plane_wave(m, k2, outer_.eps, d, p);
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) {
        if (hodlr_) hodlr_->apply(x, y); else if (prec_) prec_->apply(x, y); else T_->precondition(x, y); };
    PlaneWaveResult r; GmresResult g = gmres(A, b, r.h, &M, o.tol, o.restart, o.max_iter);
    r.iterations = g.iterations; r.residual = g.rel_residual;
    std::vector<cplx> hs(r.h.size()); for (std::size_t i = 0; i < hs.size(); ++i) hs[i] = r.h[i] - b[i];
    r.sigma_ext = extinction_cross_section(m, hs, k2, outer_.eps, d, p);
    r.forward = forward_amplitude(m, hs, k2, outer_.eps, d, p);
    return r;
}

}  // namespace cbem
