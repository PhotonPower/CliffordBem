#include "cbem/problems/scattering_problem.hpp"
#include <stdexcept>
#include "cbem/sources/fields.hpp"

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
    outer_op_ = add(m, outer.k(omega));
    std::vector<Medium> per_tri(m.size());
    for (std::size_t b = 0; b < bodies.size(); ++b)
        for (std::size_t t = mb_.body_begin[b]; t < mb_.body_begin[b + 1]; ++t) per_tri[t] = media[b];
    if (union_interior) {                        // ein Innenoperator auf der Vereinigung (gleiche achirale Medien)
        for (const auto& md : media)
            if (std::abs(md.chi) > 0 || md.eps != media[0].eps || md.mu != media[0].mu) throw std::invalid_argument("union_interior: gleiche achirale Medien noetig");
        inner_ptr_ = add(m, media[0].k(omega));
    } else {
        std::vector<const BoundaryOperator*> blocks;
        for (std::size_t b = 0; b < bodies.size(); ++b) {
            const Medium& md = media[b]; const TriangleMesh& mp = mb_.parts[b];
            if (std::abs(md.chi) > 0) {
                CauchyOperator* p = add(mp, md.k(omega, +1)); CauchyOperator* q = add(mp, md.k(omega, -1));
                chops_.push_back(std::make_unique<ChiralCauchyOperator>(*p, *q)); blocks.push_back(chops_.back().get());
            } else blocks.push_back(add(mp, md.k(omega)));
        }
        if (bodies.size() == 1) inner_ptr_ = blocks[0];
        else { inner_ = std::make_unique<BlockDiagonalOperator>(blocks, mb_.body_begin); inner_ptr_ = inner_.get(); }
    }
    T_ = std::make_unique<TransmissionOperator>(m, *inner_ptr_, *outer_op_, per_tri, outer);
}

double ScatteringProblem::hmatrix_bytes() const { double b = 0; for (auto& h : hms_) b += h->stats().bytes(); return b; }

PlaneWaveResult ScatteringProblem::solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o) const {
    const TriangleMesh& m = mb_.all; const cplx k2 = outer_.k(omega_);
    auto b = project_plane_wave(m, k2, outer_.eps, d, p);
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->precondition(x, y); };
    PlaneWaveResult r; GmresResult g = gmres(A, b, r.h, &M, o.tol, o.restart, o.max_iter);
    r.iterations = g.iterations; r.residual = g.rel_residual;
    std::vector<cplx> hs(r.h.size()); for (std::size_t i = 0; i < hs.size(); ++i) hs[i] = r.h[i] - b[i];
    r.sigma_ext = extinction_cross_section(m, hs, k2, outer_.eps, d, p);
    return r;
}

}  // namespace cbem
