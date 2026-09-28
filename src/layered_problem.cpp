#include "cbem/problems/layered_problem.hpp"
#include <stdexcept>
#include <string>
#include "cbem/sources/fields.hpp"

namespace cbem {

int add_body(LayeredGeometry& g, const TriangleMesh& surface, const Medium& inner, int parent) {
    const int r = g.add_region(inner);
    g.add_surface(surface, r, parent);
    return r;
}

std::vector<int> add_layered_body(LayeredGeometry& g, const std::vector<TriangleMesh>& surf, const std::vector<Medium>& media, int parent) {
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

std::vector<int> add_coated_body(LayeredGeometry& g, const TriangleMesh& surface, const Medium& core,
                                 const std::vector<Coating>& coatings, bool outward, int parent) {
    const std::size_t L = coatings.size();
    std::vector<real> t(L + 1, 0.0);                           // t_j = d_0 + ... + d_{j-1}
    for (std::size_t j = 0; j < L; ++j) {
        if (!(coatings[j].thickness > 0)) throw std::invalid_argument("add_coated_body: Schichtdicke muss positiv sein");
        t[j + 1] = t[j] + coatings[j].thickness;
    }
    std::vector<TriangleMesh> surf; std::vector<Medium> media;
    for (std::size_t l = 0; l <= L; ++l) {                     // von aussen nach innen
        const real off = t[L - l] - (outward ? 0.0 : t[L]);
        surf.push_back(off == 0.0 ? surface : offset_surface(surface, off));
        media.push_back(l < L ? coatings[L - 1 - l].medium : core);
    }
    return add_layered_body(g, surf, media, parent);
}

LayeredTransmissionOperator::LayeredTransmissionOperator(std::vector<Region> regions, std::vector<Mat8> J)
    : R_(std::move(regions)), J_(std::move(J)) {
    P_.resize(J_.size());
    for (std::size_t t = 0; t < J_.size(); ++t) {
        Mat8 A = J_[t]; for (int i = 0; i < 8; ++i) A[i * 8 + i] += 1.0;
        P_[t] = inverse8(A); for (auto& v : P_[t]) v *= 2.0;
    }
}

void LayeredTransmissionOperator::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    const std::size_t N = J_.size();
    std::vector<cplx> Jx(8 * N);
    for (std::size_t t = 0; t < N; ++t) cbem::apply(J_[t], &x[8 * t], &Jx[8 * t]);
    y.assign(8 * N, cplx(0));
    std::vector<cplx> w, Ew;
    for (const Region& R : R_) {
        const std::size_t n = R.tri.size();
        w.resize(8 * n);
        for (std::size_t l = 0; l < n; ++l) {
            const cplx* u = R.is_inner[l] ? &Jx[8 * R.tri[l]] : &x[8 * R.tri[l]];
            const real sg = R.is_inner[l] ? 1.0 : -1.0;       // Normale aus dem Gebiet heraus
            for (int q = 0; q < 8; ++q) w[8 * l + q] = sg * u[q];
        }
        R.E->apply(w, Ew);
        for (std::size_t l = 0; l < n; ++l) {
            const std::size_t g = R.tri[l];
            const cplx* u = R.is_inner[l] ? &Jx[8 * g] : &x[8 * g];
            for (int q = 0; q < 8; ++q) y[8 * g + q] += 0.5 * (u[q] - Ew[8 * l + q]);
        }
    }
}

void LayeredTransmissionOperator::precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    y.resize(x.size());
    for (std::size_t t = 0; t < J_.size(); ++t) cbem::apply(P_[t], &x[8 * t], &y[8 * t]);
}

LayeredScatteringProblem::LayeredScatteringProblem(const LayeredGeometry& g, real omega, HMatrixParams hp, EntryParams ep)
    : g_(g), omega_(omega), all_(make_multibody(g.surfaces)) {
    const std::size_t S = g_.surfaces.size(), NR = g_.region_medium.size();
    if (S == 0) throw std::invalid_argument("LayeredScatteringProblem: keine Flaechen");
    if (std::abs(g_.region_medium[0].chi) > 0) throw std::invalid_argument("LayeredScatteringProblem: chirales Aussenmedium nicht unterstuetzt");
    for (std::size_t s = 0; s < S; ++s) {
        const int a = g_.inside[s], b = g_.outside[s];
        if (a <= 0 || a >= static_cast<int>(NR) || b < 0 || b >= static_cast<int>(NR) || a == b)
            throw std::invalid_argument("LayeredScatteringProblem: ungueltige Gebiete an Flaeche " + std::to_string(s));
    }
    const TriangleMesh& m = all_.all;
    std::vector<Mat8> J(m.size());
    for (std::size_t s = 0; s < S; ++s)
        for (std::size_t t = all_.body_begin[s]; t < all_.body_begin[s + 1]; ++t)
            J[t] = transmission_map(m.normal[t], g_.region_medium[g_.inside[s]], g_.region_medium[g_.outside[s]]);

    auto add = [&](const TriangleMesh& mesh, cplx k) -> CauchyOperator* {
        ents_.push_back(std::make_unique<KernelEntries>(mesh, k, ep));
        hms_.push_back(std::make_unique<KernelHMatrix>(*ents_.back(), hp));
        cops_.push_back(std::make_unique<CauchyOperator>(mesh, *hms_.back()));
        return cops_.back().get();
    };
    std::vector<LayeredTransmissionOperator::Region> regions;
    for (std::size_t r = 0; r < NR; ++r) {
        std::vector<TriangleMesh> parts; LayeredTransmissionOperator::Region R;
        for (std::size_t s = 0; s < S; ++s) {
            const bool in = g_.inside[s] == static_cast<int>(r), out = g_.outside[s] == static_cast<int>(r);
            if (!in && !out) continue;
            parts.push_back(g_.surfaces[s]);
            for (std::size_t t = all_.body_begin[s]; t < all_.body_begin[s + 1]; ++t) { R.tri.push_back(t); R.is_inner.push_back(in ? 1 : 0); }
        }
        if (parts.empty()) {
            if (r == 0) throw std::invalid_argument("LayeredScatteringProblem: keine Flaeche am Aussenraum");
            continue;                                              // unbenutztes Gebiet
        }
        region_mesh_.push_back(std::make_unique<TriangleMesh>(make_multibody(parts).all));
        const TriangleMesh& rm = *region_mesh_.back();
        const Medium& md = g_.region_medium[r];
        if (std::abs(md.chi) > 0) {
            CauchyOperator* p = add(rm, md.k(omega, +1));
            CauchyOperator* q = add(rm, md.k(omega, -1));
            chops_.push_back(std::make_unique<ChiralCauchyOperator>(*p, *q));
            R.E = chops_.back().get();
        } else R.E = add(rm, md.k(omega));
        if (r == 0) ext_tri_ = R.tri;
        regions.push_back(std::move(R));
    }
    T_ = std::make_unique<LayeredTransmissionOperator>(std::move(regions), std::move(J));
}

LayeredResult LayeredScatteringProblem::solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o) const {
    const Medium& ext = g_.region_medium[0]; const cplx k = ext.k(omega_);
    const TriangleMesh& me = *region_mesh_[0];                  // Rand des Aussenraums
    const std::vector<cplx> be = project_plane_wave(me, k, ext.eps, d, p);
    std::vector<cplx> b(T_->size(), cplx(0));
    for (std::size_t l = 0; l < ext_tri_.size(); ++l) for (int q = 0; q < 8; ++q) b[8 * ext_tri_[l] + q] = be[8 * l + q];
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->precondition(x, y); };
    LayeredResult r; GmresResult gr = gmres(A, b, r.h, &M, o.tol, o.restart, o.max_iter);
    r.iterations = gr.iterations; r.residual = gr.rel_residual;
    std::vector<cplx> hs(be.size());
    for (std::size_t l = 0; l < ext_tri_.size(); ++l) for (int q = 0; q < 8; ++q) hs[8 * l + q] = r.h[8 * ext_tri_[l] + q] - be[8 * l + q];
    r.sigma_ext = extinction_cross_section(me, hs, k, ext.eps, d, p);
    r.forward = forward_amplitude(me, hs, k, ext.eps, d, p);
    return r;
}

double LayeredScatteringProblem::hmatrix_bytes() const { double b = 0; for (auto& h : hms_) b += h->stats().bytes(); return b; }
std::size_t LayeredScatteringProblem::near_pairs() const { std::size_t n = 0; for (auto& e : ents_) n += e->near_pairs(); return n; }
double LayeredScatteringProblem::near_seconds() const { double t = 0; for (auto& e : ents_) t += e->near_seconds(); return t; }

}  // namespace cbem
