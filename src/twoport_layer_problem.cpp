#include "cbem/problems/twoport_layer_problem.hpp"
#include <cmath>
#include <stdexcept>
#include "cbem/sources/chiral_incidence.hpp"
#include "cbem/sources/fields.hpp"

namespace cbem {

namespace {
// Inverse einer komplexen n x n-Matrix (Gauss-Jordan mit Pivotsuche), zeilenweise
std::vector<cplx> inverse_n(std::vector<cplx> A, int n) {
    std::vector<cplx> I(n * n, cplx(0)); for (int i = 0; i < n; ++i) I[i * n + i] = 1.0;
    for (int c = 0; c < n; ++c) {
        int piv = c; for (int r = c + 1; r < n; ++r) if (std::abs(A[r * n + c]) > std::abs(A[piv * n + c])) piv = r;
        if (std::abs(A[piv * n + c]) < 1e-300) throw std::runtime_error("TwoPort: singulaerer Vorkonditionierer");
        for (int q = 0; q < n; ++q) { std::swap(A[c * n + q], A[piv * n + q]); std::swap(I[c * n + q], I[piv * n + q]); }
        const cplx dd = A[c * n + c]; for (int q = 0; q < n; ++q) { A[c * n + q] /= dd; I[c * n + q] /= dd; }
        for (int r = 0; r < n; ++r) if (r != c) { const cplx f = A[r * n + c]; if (f == cplx(0)) continue;
            for (int q = 0; q < n; ++q) { A[r * n + q] -= f * A[c * n + q]; I[r * n + q] -= f * I[c * n + q]; } }
    }
    return I;
}
Multivector apply8(const Mat8& M, const Multivector& x) { Multivector y; cbem::apply(M, x.c.data(), y.c.data()); return y; }
Multivector central(cplx a, cplx b) { return Multivector::blade(0, 0.5 * (a + b)) + Multivector::blade(7, cplx(0, 0.5) * (a - b)); }  // a P+ + b P-
}  // namespace

std::vector<TriangleMesh> TwoPortLayerProblem::layer_surfaces(const TriangleMesh& core_surface, const std::vector<Coating>& layers) {
    std::vector<TriangleMesh> S{core_surface}; real t = 0;
    for (const Coating& c : layers) { t += c.thickness; S.push_back(offset_surface(core_surface, t)); }
    return S;
}

TwoPortLayerProblem::TwoPortLayerProblem(const TriangleMesh& core_surface, const TriangleMesh& outer_surface, const Medium& core,
                                         const Medium& layer, real d, real omega, Medium outer, HMatrixParams hp, EntryParams ep,
                                         TwoPortOptions opt)
    : outer_m_(outer), omega_(omega), opt_(opt) {
    build({TwoPortBody{{core_surface, outer_surface}, core, {Coating{d, layer}}}}, hp, ep);
}

TwoPortLayerProblem::TwoPortLayerProblem(const std::vector<TriangleMesh>& surfaces, const Medium& core, const std::vector<Coating>& layers,
                                         real omega, Medium outer, HMatrixParams hp, EntryParams ep, TwoPortOptions opt)
    : outer_m_(outer), omega_(omega), opt_(opt) {
    build({TwoPortBody{surfaces, core, layers}}, hp, ep);
}

TwoPortLayerProblem::TwoPortLayerProblem(const std::vector<TwoPortBody>& bodies, real omega, Medium outer, HMatrixParams hp,
                                         EntryParams ep, TwoPortOptions opt)
    : outer_m_(outer), omega_(omega), opt_(opt) {
    build(bodies, hp, ep);
}

void TwoPortLayerProblem::setup_body(Body& Bd) {
    Bd.L_in = Bd.lay.size();
    for (auto& s : Bd.S) if (s.normal.size() != s.size()) s.compute_geometry();
    Bd.N = Bd.S[0].size();
    for (auto& s : Bd.S) if (s.size() != Bd.N || s.P.size() != Bd.S[0].P.size()) throw std::invalid_argument("TwoPort: Flaechen mit verschiedener Topologie");
    // Unterteilung dicker Schichten nach der groessten Hauptkruemmung der Kernflaeche
    if (opt_.split > 0) {
        SurfaceFV f0(Bd.S[0]); real kmax = 0;
        for (const auto& Sm : f0.shape) {                               // groesster Eigenwert-Betrag von S (Potenziteration)
            real v[3] = {1, 0.7, 0.3}, lam = 0;
            for (int it = 0; it < 30; ++it) {
                real w[3] = {0, 0, 0}; for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) w[a] += Sm[3 * a + b] * v[b];
                const real nw = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]); if (!(nw > 0)) break;
                lam = nw; for (int a = 0; a < 3; ++a) v[a] = w[a] / nw;
            }
            kmax = std::max(kmax, lam);
        }
        if (kmax > 0) {
            const real dmax = opt_.split / kmax;
            std::vector<TriangleMesh> S2{Bd.S[0]}; std::vector<Layer> L2;
            for (std::size_t l = 0; l < Bd.lay.size(); ++l) {
                const int parts = std::max(1, static_cast<int>(std::ceil(Bd.lay[l].d / dmax - 1e-9)));
                for (int q = 1; q <= parts; ++q) {
                    Layer Ly = Bd.lay[l]; Ly.d = Bd.lay[l].d / parts; L2.push_back(Ly);
                    if (q == parts) { S2.push_back(Bd.S[l + 1]); continue; }
                    TriangleMesh Mq = Bd.S[l + 1]; const real a = real(q) / parts;   // lineare Interpolation der Knoten
                    for (std::size_t v = 0; v < Mq.P.size(); ++v) Mq.P[v] = Bd.S[l].P[v] * (1 - a) + Bd.S[l + 1].P[v] * a;
                    Mq.compute_geometry(); S2.push_back(Mq);
                }
            }
            Bd.S = S2; Bd.lay = L2;
        }
    }
    Bd.L = Bd.lay.size();
    const std::size_t N = Bd.N, L = Bd.L;
    for (std::size_t j = 0; j < L; ++j) Bd.fv.push_back(std::make_unique<SurfaceFV>(Bd.S[j]));
    Bd.sq.assign(L + 1, std::vector<real>(N));
    for (std::size_t j = 0; j <= L; ++j) for (std::size_t t = 0; t < N; ++t) Bd.sq[j][t] = std::sqrt(Bd.S[j].area[t]);
    const TriangleMesh& mb = Bd.S[0];
    Bd.Jb.resize(N); Bd.Jbs.resize(N);
    for (std::size_t t = 0; t < N; ++t) { Bd.Jb[t] = transmission_map(mb.normal[t], Bd.lay[0].m, Bd.core); Bd.Jbs[t] = transmission_map(Bd.fv[0]->nsm[t], Bd.lay[0].m, Bd.core); }
    for (std::size_t l = 0; l < L; ++l) {
        Layer& Ly = Bd.lay[l];
        if (!(Ly.d > 0)) throw std::invalid_argument("TwoPort: Schichtdicke muss positiv sein");
        real smax_lap = 0;
        for (std::size_t t = 0; t < N; ++t) { real s = 0; for (real w : Bd.fv[l]->ring[t].lw) s += std::abs(w); smax_lap = std::max(smax_lap, 2 * s); }
        const Medium& above = l + 1 < L ? Bd.lay[l + 1].m : outer_m_;
        Ly.Jtop.resize(N); Ly.Jtops.resize(N);
        for (std::size_t t = 0; t < N; ++t) { Ly.Jtop[t] = transmission_map(Bd.S[l + 1].normal[t], Ly.m, above); Ly.Jtops[t] = transmission_map(Bd.fv[l]->nsm[t], Ly.m, above); }
        Ly.numid = 0.5 * Ly.d;
        // Pole von g_l: sigma_m = ((2m+1) pi / d)^2 bis sigma_M >= 10 s_max, Rest linear in s
        const real smax = smax_lap + std::max(std::norm(Ly.m.k(omega_, +1)), std::norm(Ly.m.k(omega_, -1)));
        if (opt_.poles >= 0) Ly.M = opt_.poles;
        else { const real q = Ly.d * std::sqrt(10 * smax) / pi; Ly.M = std::min(200, std::max(0, static_cast<int>(std::ceil((q - 1) / 2)))); }
        real s2 = pi * pi / 8, s4 = std::pow(pi, 4) / 96;
        for (int m = 0; m < Ly.M; ++m) { const real o = 2 * m + 1; Ly.sigma.push_back(std::pow(o * pi / Ly.d, 2)); s2 -= 1 / (o * o); s4 -= 1 / std::pow(o, 4); }
        Ly.tail0 = (4 / Ly.d) * (Ly.d * Ly.d / (pi * pi)) * s2;
        Ly.tail1 = (4 / Ly.d) * std::pow(Ly.d / pi, 4) * s4;
    }
}

void TwoPortLayerProblem::build(std::vector<TwoPortBody> bodies, HMatrixParams hp, EntryParams ep) {
    if (bodies.empty()) throw std::invalid_argument("TwoPort: keine Koerper");
    for (auto& bd : bodies) {
        if (bd.layers.empty() || bd.surfaces.size() != bd.layers.size() + 1) throw std::invalid_argument("TwoPort: je Schicht eine Flaeche mehr als Schichten");
        auto Bd = std::make_unique<Body>(); Bd->S = bd.surfaces; Bd->core = bd.core;
        for (const Coating& c : bd.layers) { Layer Ly; Ly.m = c.medium; Ly.d = c.thickness; Bd->lay.push_back(Ly); }
        setup_body(*Bd);
        B_.push_back(std::move(Bd));
    }
    // Positionen: [X_L aller Koerper | v aller Koerper | Zwischenspuren je Koerper]
    Ntot_ = 0; for (auto& b : B_) { b->pre = Ntot_; Ntot_ += b->N; }
    std::size_t off = 16 * Ntot_;
    for (auto& b : B_) {
        b->offXL = 8 * b->pre; b->offV = 8 * Ntot_ + 8 * b->pre;
        b->offInt.clear(); for (std::size_t j = 1; j < b->L; ++j) { b->offInt.push_back(off); off += 8 * b->N; }
    }
    size_ = off;
    // Operatoren: E_2 auf der Vereinigung der Aussenflaechen, E_1 blockdiagonal auf den Kernen
    std::vector<TriangleMesh> outs; for (auto& b : B_) outs.push_back(b->S.back());
    { std::vector<const TriangleMesh*> ps; for (auto& o : outs) ps.push_back(&o); require_separated_all(ps, {}, "TwoPortLayerProblem (Aussenflaechen)"); }
    outer_ = make_multibody(outs);
    auto add = [&](const TriangleMesh& mesh, cplx k) {
        ents_.push_back(std::make_unique<KernelEntries>(mesh, k, ep));
        hms_.push_back(std::make_unique<KernelHMatrix>(*ents_.back(), hp));
        cops_.push_back(std::make_unique<CauchyOperator>(mesh, *hms_.back()));
        return cops_.back().get();
    };
    if (std::abs(outer_m_.chi) > 0) {
        const CauchyOperator* p = add(outer_.all, outer_m_.k(omega_, +1)); const CauchyOperator* q = add(outer_.all, outer_m_.k(omega_, -1));
        chops_.push_back(std::make_unique<ChiralCauchyOperator>(*p, *q)); E2op_ = chops_.back().get();
    } else E2op_ = add(outer_.all, outer_m_.k(omega_));
    std::vector<const BoundaryOperator*> inner; std::vector<std::size_t> begin{0};
    for (auto& b : B_) {
        const TriangleMesh& mb = b->S[0];
        if (std::abs(b->core.chi) > 0) {
            const CauchyOperator* p = add(mb, b->core.k(omega_, +1)); const CauchyOperator* q = add(mb, b->core.k(omega_, -1));
            chops_.push_back(std::make_unique<ChiralCauchyOperator>(*p, *q)); inner.push_back(chops_.back().get());
        } else inner.push_back(add(mb, b->core.k(omega_)));
        begin.push_back(begin.back() + b->N);
    }
    if (B_.size() == 1) E1_ = inner[0];
    else { E1bd_ = std::make_unique<BlockDiagonalOperator>(inner, begin); E1_ = E1bd_.get(); }
    // Vorkonditionierer: lokaler Teil je Dreieck (Zeile I: 1/2 X_L + 1/2 rho v; Zeile l: Jtop X_l - u_bot, skaliert mit sqrt|tau_0|)
    for (auto& bp : B_) {
        Body& Bd = *bp; const std::size_t L = Bd.L; const int nb = static_cast<int>(8 * (L + 1));
        auto lpos = [&](std::size_t j) -> std::size_t { return j == L ? 0 : j == 0 ? 1 : j + 1; };   // lokale Blockposition
        Bd.Pinv.resize(Bd.N);
        for (std::size_t t = 0; t < Bd.N; ++t) {
            std::vector<cplx> A(nb * nb, cplx(0));
            auto put = [&](std::size_t row, std::size_t col, const Mat8* M, cplx s) {
                for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
                    A[(8 * row + i) * nb + 8 * col + j] += s * (M ? (*M)[i * 8 + j] : cplx(i == j ? 1.0 : 0.0));
            };
            const real rho = Bd.sq[L][t] / Bd.sq[0][t];
            put(0, lpos(L), nullptr, 0.5); put(0, lpos(0), nullptr, 0.5 * rho);
            for (std::size_t l = 1; l <= L; ++l) {
                put(l, lpos(l), &Bd.lay[l - 1].Jtop[t], Bd.sq[0][t] / Bd.sq[l][t]);
                if (l == 1) put(l, lpos(0), &Bd.Jb[t], -1.0);
                else put(l, lpos(l - 1), nullptr, -Bd.sq[0][t] / Bd.sq[l - 1][t]);
            }
            Bd.Pinv[t] = inverse_n(A, nb);
        }
    }
}

void TwoPortLayerProblem::resolvent(const SurfaceFV& fv, const std::vector<Multivector>& z, cplx shift, std::vector<Multivector>& x) const {
    const std::size_t N = z.size();
    const auto& R = fv.ring;
    std::vector<cplx> b(8 * N), sol;
    for (std::size_t t = 0; t < N; ++t) for (int c = 0; c < 8; ++c) b[8 * t + c] = z[t].c[c];
    LinOp A = [&](const std::vector<cplx>& v, std::vector<cplx>& y) {
        y.resize(v.size());
        for (std::size_t t = 0; t < N; ++t)
            for (int c = 0; c < 8; ++c) {
                cplx lap = 0; for (std::size_t j = 0; j < R[t].nb.size(); ++j) lap += R[t].lw[j] * (v[8 * R[t].nb[j] + c] - v[8 * t + c]);
                y[8 * t + c] = -lap + shift * v[8 * t + c];
            }
    };
    std::vector<cplx> diag(N);
    for (std::size_t t = 0; t < N; ++t) { real s = 0; for (real w : R[t].lw) s += w; diag[t] = s + shift; }
    LinOp P = [&](const std::vector<cplx>& v, std::vector<cplx>& y) { y.resize(v.size()); for (std::size_t i = 0; i < v.size(); ++i) y[i] = v[i] / diag[i / 8]; };
    GmresResult g = gmres(A, b, sol, &P, opt_.inner_tol, 60, 2000);
    inner_its_ += g.iterations; ++inner_calls_;
    x.resize(N);
    for (std::size_t t = 0; t < N; ++t) for (int c = 0; c < 8; ++c) x[t].c[c] = sol[8 * t + c];
}

void TwoPortLayerProblem::apply_g(const Body& Bd, std::size_t l, const std::vector<Multivector>& z, std::vector<Multivector>& gz) const {
    const std::size_t N = Bd.N; const Layer& Ly = Bd.lay[l]; const SurfaceFV& fv = *Bd.fv[l];
    const cplx kp = Ly.m.k(omega_, +1), km = Ly.m.k(omega_, -1);
    const bool chiral = std::abs(Ly.m.chi) > 0;
    const Multivector Pp = central(1.0, 0.0), Pm = central(0.0, 1.0), K2 = central(kp * kp, km * km);
    const auto& R = fv.ring;
    gz.resize(N);
    for (std::size_t t = 0; t < N; ++t) {                                   // Rest: tail0 z - tail1 s z, s z = -K^2 z - Delta z
        Multivector lap;
        for (std::size_t j = 0; j < R[t].nb.size(); ++j) lap = lap + (z[R[t].nb[j]] - z[t]) * R[t].lw[j];
        const Multivector sz = (K2 * z[t]) * (-1.0) - lap;
        gz[t] = z[t] * Ly.tail0 - sz * Ly.tail1;
    }
    std::vector<Multivector> zp, zm, x;
    if (chiral) { zp.resize(N); zm.resize(N); for (std::size_t t = 0; t < N; ++t) { zp[t] = Pp * z[t]; zm[t] = Pm * z[t]; } }
    const real w = 4 / Ly.d;
    for (int m = 0; m < Ly.M; ++m) {
        if (!chiral) { resolvent(fv, z, Ly.sigma[m] - kp * kp, x); for (std::size_t t = 0; t < N; ++t) gz[t] = gz[t] + x[t] * w; }
        else {
            resolvent(fv, zp, Ly.sigma[m] - kp * kp, x); for (std::size_t t = 0; t < N; ++t) gz[t] = gz[t] + x[t] * w;
            resolvent(fv, zm, Ly.sigma[m] - km * km, x); for (std::size_t t = 0; t < N; ++t) gz[t] = gz[t] + x[t] * w;
        }
    }
}

void TwoPortLayerProblem::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    const std::size_t Nt = Ntot_;
    std::vector<cplx> XL(x.begin(), x.begin() + 8 * Nt), V(x.begin() + 8 * Nt, x.begin() + 16 * Nt), E2X, E1V;
    E2op_->apply(XL, E2X); E1_->apply(V, E1V);
    y.assign(size_, cplx(0));
    for (const auto& bp : B_) {                                             // Zeile I
        const Body& Bd = *bp;
        for (std::size_t t = 0; t < Bd.N; ++t) {
            const real rho = Bd.sq[Bd.L][t] / Bd.sq[0][t];
            for (int c = 0; c < 8; ++c) { const std::size_t i = 8 * (Bd.pre + t) + c; y[i] = 0.5 * (XL[i] + E2X[i]) + rho * 0.5 * (V[i] - E1V[i]); }
        }
    }
    for (const auto& bp : B_) {                                             // Zweitor-Zeilen je Koerper
        const Body& Bd = *bp; const std::size_t N = Bd.N, L = Bd.L;
        std::vector<std::vector<Multivector>> X(L + 1, std::vector<Multivector>(N));
        for (std::size_t j = 0; j <= L; ++j) {
            const std::size_t o = Bd.pos(j);
            for (std::size_t t = 0; t < N; ++t) { const real s = 1 / Bd.sq[j][t]; for (int c = 0; c < 8; ++c) X[j][t].c[c] = x[o + 8 * t + c] * s; }
        }
        std::vector<Multivector> ut(N), ubot(N), z(N), DF, DSF, DDF, Bz(N), gBz;
        for (std::size_t l = 1; l <= L; ++l) {
            const Layer& Ly = Bd.lay[l - 1];
            for (std::size_t t = 0; t < N; ++t) {
                ut[t] = apply8(Ly.Jtop[t], X[l][t]);
                ubot[t] = l == 1 ? apply8(Bd.Jb[t], X[0][t]) : X[l - 1][t];
                z[t] = apply8(Ly.Jtops[t], X[l][t]) + (l == 1 ? apply8(Bd.Jbs[t], X[0][t]) : X[l - 1][t]);
            }
            const SurfaceFV& fv = *Bd.fv[l - 1];                              // untere Flaeche der Schicht l
            fv.dirac_fit(z, DF, DSF, DDF);
            const Multivector iK = central(Ly.m.k(omega_, +1), Ly.m.k(omega_, -1)) * cplx(0, 1);
            for (std::size_t t = 0; t < N; ++t) {
                Multivector Xv = iK * z[t] - DF[t];
                if (opt_.curvature) Xv = Xv + DSF[t] * Ly.numid;
                Bz[t] = Multivector::vector(fv.nsm[t]) * Xv;                   // B in der Mitte der Schicht l
            }
            apply_g(Bd, l - 1, Bz, gBz);
            const std::size_t o = Bd.row(l);
            for (std::size_t t = 0; t < N; ++t) {
                const Multivector r = ut[t] - ubot[t] - gBz[t];
                for (int c = 0; c < 8; ++c) y[o + 8 * t + c] = r.c[c] * Bd.sq[0][t];
            }
        }
    }
}

void TwoPortLayerProblem::precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    y.assign(x.size(), cplx(0));
    for (const auto& bp : B_) {
        const Body& Bd = *bp; const std::size_t L = Bd.L; const int nb = static_cast<int>(8 * (L + 1));
        std::vector<cplx> in(nb), out(nb);
        std::vector<std::size_t> blk(L + 1);                                // Blockposition s -> Offset (s = 0: X_L, 1: v, s >= 2: X_{s-1})
        blk[0] = Bd.offXL; blk[1] = Bd.offV; for (std::size_t s = 2; s <= L; ++s) blk[s] = Bd.offInt[s - 2];
        for (std::size_t t = 0; t < Bd.N; ++t) {
            for (std::size_t r = 0; r <= L; ++r) for (int c = 0; c < 8; ++c) in[8 * r + c] = x[blk[r] + 8 * t + c];   // Zeile r liegt bei blk[r]
            for (int r = 0; r < nb; ++r) { out[r] = 0; for (int c = 0; c < nb; ++c) out[r] += Bd.Pinv[t][r * nb + c] * in[c]; }
            for (std::size_t s = 0; s <= L; ++s) for (int c = 0; c < 8; ++c) y[blk[s] + 8 * t + c] = out[8 * s + c];
        }
    }
}

LayeredResult TwoPortLayerProblem::solve_plane_wave(const Vec3& dir, const CVec3& p, const SolveOptions& o) const {
    const PlaneWaveIncidence inc = plane_wave_incidence(outer_m_, omega_, dir, p);
    const TriangleMesh& ma = outer_.all;
    const std::vector<cplx> ba = project_plane_wave(ma, inc.k, outer_m_.eps, dir, p);
    LayeredResult r = solve_rhs(ba, o);
    std::vector<cplx> hs(8 * Ntot_); for (std::size_t i = 0; i < 8 * Ntot_; ++i) hs[i] = r.h[i] - ba[i];
    r.sigma_ext = extinction_in_medium(ma, hs, outer_m_, inc, dir, p);
    r.forward = forward_amplitude_in_medium(ma, hs, outer_m_, inc, dir, p);
    return r;
}

LayeredResult TwoPortLayerProblem::solve_rhs(const std::vector<cplx>& ba, const SolveOptions& o) const {
    std::vector<cplx> b(size_, cplx(0)); std::copy(ba.begin(), ba.end(), b.begin());
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { precondition(x, y); };
    LayeredResult r; GmresResult g = gmres(A, b, r.h, &M, o.tol, o.restart, o.max_iter);
    r.iterations = g.iterations; r.residual = g.rel_residual;
    return r;
}

}  // namespace cbem
