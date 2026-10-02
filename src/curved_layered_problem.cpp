#include "cbem/problems/curved_layered_problem.hpp"

#include <stdexcept>
#include <string>

#include "cbem/problems/curved_problem.hpp"
#include "cbem/solvers/gmres.hpp"
#include "cbem/sources/chiral_incidence.hpp"

namespace cbem {

int add_body(CurvedLayeredGeometry& g, const QuadraticMesh& surface, const Medium& inner, int parent) {
    const int r = g.add_region(inner);
    g.add_surface(surface, r, parent);
    return r;
}

std::vector<int> add_layered_body(CurvedLayeredGeometry& g, const std::vector<QuadraticMesh>& surf, const std::vector<Medium>& media, int parent) {
    if (surf.empty() || surf.size() != media.size()) throw std::invalid_argument("add_layered_body: ein Medium je Flaeche");
    std::vector<int> regions;
    int out = parent;
    for (std::size_t l = 0; l < surf.size(); ++l) {
        const int r = g.add_region(media[l]);
        g.add_surface(surf[l], r, out);
        regions.push_back(r); out = r;
    }
    return regions;
}

std::vector<int> add_coated_body(CurvedLayeredGeometry& g, const QuadraticMesh& surface, const Medium& core,
                                 const std::vector<Coating>& coatings, bool outward, int parent) {
    const std::size_t L = coatings.size();
    std::vector<real> t(L + 1, 0.0);                           // t_j = d_0 + ... + d_{j-1}
    for (std::size_t j = 0; j < L; ++j) {
        if (!(coatings[j].thickness > 0)) throw std::invalid_argument("add_coated_body: Schichtdicke muss positiv sein");
        t[j + 1] = t[j] + coatings[j].thickness;
    }
    std::vector<QuadraticMesh> surf; std::vector<Medium> media;
    for (std::size_t l = 0; l <= L; ++l) {                     // von aussen nach innen
        const real off = t[L - l] - (outward ? 0.0 : t[L]);
        surf.push_back(off == 0.0 ? surface : offset_surface(surface, off));
        media.push_back(l < L ? coatings[L - 1 - l].medium : core);
    }
    return add_layered_body(g, surf, media, parent);
}

CurvedLayeredTransmissionOperator::CurvedLayeredTransmissionOperator(std::vector<Region> regions, std::vector<std::vector<cplx>> JG,
                                                                     std::vector<std::vector<cplx>> P)
    : R_(std::move(regions)), JG_(std::move(JG)), P_(std::move(P)) {}

void CurvedLayeredTransmissionOperator::apply_J(const std::vector<cplx>& x, std::vector<cplx>& y) const { curved_block_apply(JG_, x, y); }

void CurvedLayeredTransmissionOperator::precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const { curved_block_apply(P_, x, y); }

void CurvedLayeredTransmissionOperator::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    std::vector<cplx> Jx;
    curved_block_apply(JG_, x, Jx);
    y.assign(x.size(), cplx(0));
    std::vector<cplx> w, Ew;
    for (const Region& R : R_) {
        const std::size_t n = R.el.size();
        w.resize(24 * n);
        for (std::size_t l = 0; l < n; ++l) {
            const cplx* u = R.is_inner[l] ? &Jx[24 * R.el[l]] : &x[24 * R.el[l]];
            const real sg = R.is_inner[l] ? 1.0 : -1.0;       // Normale aus dem Gebiet heraus
            for (int q = 0; q < 24; ++q) w[24 * l + q] = sg * u[q];
        }
        R.E->apply(w, Ew);
        for (std::size_t l = 0; l < n; ++l) {
            const std::size_t g = R.el[l];
            const cplx* u = R.is_inner[l] ? &Jx[24 * g] : &x[24 * g];
            for (int q = 0; q < 24; ++q) y[24 * g + q] += 0.5 * (u[q] - Ew[24 * l + q]);
        }
    }
}

CurvedLayeredScatteringProblem::CurvedLayeredScatteringProblem(const CurvedLayeredGeometry& g, real omega, HMatrixParams hp, EntryParams ep,
                                                               CurvedNearParams np)
    : g_(g), omega_(omega) {
    const std::size_t S = g_.surfaces.size(), NR = g_.region_medium.size();
    if (S == 0) throw std::invalid_argument("CurvedLayeredScatteringProblem: keine Flaechen");
    for (std::size_t r = 0; r < NR; ++r) {                        // Flaechen, die an dasselbe Gebiet grenzen, duerfen sich nicht beruehren
        std::vector<const TriangleMesh*> ps;
        for (std::size_t s = 0; s < S; ++s) if (g_.outside[s] == static_cast<int>(r)) ps.push_back(&g_.surfaces[s].flat);
        if (ps.size() > 1) require_separated_all(ps, {}, "CurvedLayeredScatteringProblem, Gebiet " + std::to_string(r));
    }
    for (std::size_t s = 0; s < S; ++s) {
        const int a = g_.inside[s], b = g_.outside[s];
        if (a <= 0 || a >= static_cast<int>(NR) || b < 0 || b >= static_cast<int>(NR) || a == b)
            throw std::invalid_argument("CurvedLayeredScatteringProblem: ungueltige Gebiete an Flaeche " + std::to_string(s));
    }
    all_ = merge_quadratic(g_.surfaces, &begin_);
    std::vector<Medium> in(all_.size()), out(all_.size());
    for (std::size_t s = 0; s < S; ++s)
        for (std::size_t t = begin_[s]; t < begin_[s + 1]; ++t) { in[t] = g_.region_medium[g_.inside[s]]; out[t] = g_.region_medium[g_.outside[s]]; }
    std::vector<std::vector<cplx>> JG, P;
    curved_transmission_blocks(all_, curved_psi_matrices(all_), in, out, JG, P);

    auto add = [&](const QuadraticMesh& mesh, cplx k) -> CurvedCauchyOperator* {
        ents_.push_back(std::make_unique<CurvedKernelEntries>(mesh, k, ep, np));
        hms_.push_back(std::make_unique<CurvedHMatrix>(*ents_.back(), hp));
        cops_.push_back(std::make_unique<CurvedCauchyOperator>(mesh, *hms_.back()));
        return cops_.back().get();
    };
    std::vector<CurvedLayeredTransmissionOperator::Region> regions;
    for (std::size_t r = 0; r < NR; ++r) {
        std::vector<QuadraticMesh> parts; CurvedLayeredTransmissionOperator::Region R;
        for (std::size_t s = 0; s < S; ++s) {
            const bool is_in = g_.inside[s] == static_cast<int>(r), is_out = g_.outside[s] == static_cast<int>(r);
            if (!is_in && !is_out) continue;
            parts.push_back(g_.surfaces[s]);
            for (std::size_t t = begin_[s]; t < begin_[s + 1]; ++t) { R.el.push_back(t); R.is_inner.push_back(is_in ? 1 : 0); }
        }
        if (parts.empty()) {
            if (r == 0) throw std::invalid_argument("CurvedLayeredScatteringProblem: keine Flaeche am Aussenraum");
            continue;                                              // unbenutztes Gebiet
        }
        region_mesh_.push_back(std::make_unique<QuadraticMesh>(merge_quadratic(parts)));
        const QuadraticMesh& rm = *region_mesh_.back();
        const Medium& md = g_.region_medium[r];
        if (std::abs(md.chi) > 0) {
            CurvedCauchyOperator* p = add(rm, md.k(omega, +1));
            CurvedCauchyOperator* q = add(rm, md.k(omega, -1));
            chops_.push_back(std::make_unique<ChiralCauchyOperator>(*p, *q));
            R.E = chops_.back().get();
        } else R.E = add(rm, md.k(omega));
        if (r == 0) ext_el_ = R.el;
        regions.push_back(std::move(R));
    }
    T_ = std::make_unique<CurvedLayeredTransmissionOperator>(std::move(regions), std::move(JG), std::move(P));
}

LayeredResult CurvedLayeredScatteringProblem::solve_plane_wave(const Vec3& d0, const CVec3& p, const SolveOptions& o) const {
    const Vec3 d = d0 / norm(d0);
    const Medium& ext = g_.region_medium[0];
    const PlaneWaveIncidence inc = plane_wave_incidence(ext, omega_, d, p);   // chiral: Helizitaetswelle mit k_sigma
    const QuadraticMesh& me = *region_mesh_[0];                   // Rand des Aussenraums
    const std::vector<cplx> be = project_plane_wave_curved(me, inc.k, ext.eps, d, p);
    std::vector<cplx> b(T_->size(), cplx(0));
    for (std::size_t l = 0; l < ext_el_.size(); ++l) for (int q = 0; q < 24; ++q) b[24 * ext_el_[l] + q] = be[24 * l + q];
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->precondition(x, y); };
    LayeredResult r; const GmresResult gr = gmres(A, b, r.h, &M, o.tol, o.restart, o.max_iter);
    r.iterations = gr.iterations; r.residual = gr.rel_residual;
    std::vector<cplx> hs(be.size());
    for (std::size_t l = 0; l < ext_el_.size(); ++l) for (int q = 0; q < 24; ++q) hs[24 * l + q] = r.h[24 * ext_el_[l] + q] - be[24 * l + q];
    const std::vector<cplx> hc = helicity_part(hs, inc.proj);
    r.sigma_ext = extinction_cross_section_curved(me, hc, inc.k, ext.eps, d, p);
    r.forward = forward_amplitude_curved(me, hc, inc.k, ext.eps, d, p);
    return r;
}

double CurvedLayeredScatteringProblem::hmatrix_bytes() const { double b = 0; for (auto& h : hms_) b += h->stats().bytes(); return b; }
double CurvedLayeredScatteringProblem::near_seconds() const { double t = 0; for (auto& e : ents_) t += e->near_seconds(); return t; }

}  // namespace cbem
