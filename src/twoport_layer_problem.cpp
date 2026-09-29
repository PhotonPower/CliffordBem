#include "cbem/problems/twoport_layer_problem.hpp"
#include <cmath>
#include <stdexcept>
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
    : S_{core_surface, outer_surface}, core_(core), outer_(outer), omega_(omega), opt_(opt) {
    lay_.resize(1); lay_[0].m = layer; lay_[0].d = d;
    build(hp, ep);
}

TwoPortLayerProblem::TwoPortLayerProblem(const std::vector<TriangleMesh>& surfaces, const Medium& core, const std::vector<Coating>& layers,
                                         real omega, Medium outer, HMatrixParams hp, EntryParams ep, TwoPortOptions opt)
    : S_(surfaces), core_(core), outer_(outer), omega_(omega), opt_(opt) {
    if (layers.empty() || surfaces.size() != layers.size() + 1) throw std::invalid_argument("TwoPort: je Schicht eine Flaeche mehr als Schichten");
    for (const Coating& c : layers) { Layer Ly; Ly.m = c.medium; Ly.d = c.thickness; lay_.push_back(Ly); }
    build(hp, ep);
}

void TwoPortLayerProblem::build(HMatrixParams hp, EntryParams ep) {
    L_in_ = lay_.size();
    for (auto& s : S_) if (s.normal.size() != s.size()) s.compute_geometry();
    N_ = S_[0].size();
    for (auto& s : S_) if (s.size() != N_ || s.P.size() != S_[0].P.size()) throw std::invalid_argument("TwoPort: Flaechen mit verschiedener Topologie");
    // Unterteilung dicker Schichten nach der groessten Hauptkruemmung der Kernflaeche
    if (opt_.split > 0) {
        SurfaceFV f0(S_[0]); real kmax = 0;
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
            std::vector<TriangleMesh> S2{S_[0]}; std::vector<Layer> L2;
            for (std::size_t l = 0; l < lay_.size(); ++l) {
                const int parts = std::max(1, static_cast<int>(std::ceil(lay_[l].d / dmax - 1e-9)));
                for (int q = 1; q <= parts; ++q) {
                    Layer Ly = lay_[l]; Ly.d = lay_[l].d / parts; L2.push_back(Ly);
                    if (q == parts) { S2.push_back(S_[l + 1]); continue; }
                    TriangleMesh Mq = S_[l + 1]; const real a = real(q) / parts;   // lineare Interpolation der Knoten
                    for (std::size_t v = 0; v < Mq.P.size(); ++v) Mq.P[v] = S_[l].P[v] * (1 - a) + S_[l + 1].P[v] * a;
                    Mq.compute_geometry(); S2.push_back(Mq);
                }
            }
            S_ = S2; lay_ = L2;
        }
    }
    L_ = lay_.size();
    if (std::abs(outer_.chi) > 0) throw std::invalid_argument("TwoPort: chirales Aussenmedium nicht unterstuetzt");
    const TriangleMesh& mb = S_[0]; const TriangleMesh& ma = S_[L_];
    for (std::size_t j = 0; j < L_; ++j) fv_.push_back(std::make_unique<SurfaceFV>(S_[j]));
    K2_ = std::make_unique<KernelEntries>(ma, outer_.k(omega_), ep); H2_ = std::make_unique<KernelHMatrix>(*K2_, hp);
    E2_ = std::make_unique<CauchyOperator>(ma, *H2_);
    if (std::abs(core_.chi) > 0) {
        K1p_ = std::make_unique<KernelEntries>(mb, core_.k(omega_, +1), ep); H1p_ = std::make_unique<KernelHMatrix>(*K1p_, hp);
        K1m_ = std::make_unique<KernelEntries>(mb, core_.k(omega_, -1), ep); H1m_ = std::make_unique<KernelHMatrix>(*K1m_, hp);
        E1p_ = std::make_unique<CauchyOperator>(mb, *H1p_); E1m_ = std::make_unique<CauchyOperator>(mb, *H1m_);
        E1ch_ = std::make_unique<ChiralCauchyOperator>(*E1p_, *E1m_); E1_ = E1ch_.get();
    } else {
        K1p_ = std::make_unique<KernelEntries>(mb, core_.k(omega_), ep); H1p_ = std::make_unique<KernelHMatrix>(*K1p_, hp);
        E1p_ = std::make_unique<CauchyOperator>(mb, *H1p_); E1_ = E1p_.get();
    }
    sq_.assign(L_ + 1, std::vector<real>(N_));
    for (std::size_t j = 0; j <= L_; ++j) for (std::size_t t = 0; t < N_; ++t) sq_[j][t] = std::sqrt(S_[j].area[t]);
    // Transmissionsabbildungen: Kern -> m_1 an Gamma_0; X_l (m_{l+1}) -> m_l an Gamma_l
    Jb_.resize(N_); Jbs_.resize(N_);
    for (std::size_t t = 0; t < N_; ++t) { Jb_[t] = transmission_map(mb.normal[t], lay_[0].m, core_); Jbs_[t] = transmission_map(fv_[0]->nsm[t], lay_[0].m, core_); }
    for (std::size_t l = 0; l < L_; ++l) {
        real smax_lap = 0;
        for (std::size_t t = 0; t < N_; ++t) { real s = 0; for (real w : fv_[l]->ring[t].lw) s += std::abs(w); smax_lap = std::max(smax_lap, 2 * s); }
        Layer& Ly = lay_[l];
        if (!(Ly.d > 0)) throw std::invalid_argument("TwoPort: Schichtdicke muss positiv sein");
        const Medium& above = l + 1 < L_ ? lay_[l + 1].m : outer_;
        Ly.Jtop.resize(N_); Ly.Jtops.resize(N_);
        for (std::size_t t = 0; t < N_; ++t) { Ly.Jtop[t] = transmission_map(S_[l + 1].normal[t], Ly.m, above); Ly.Jtops[t] = transmission_map(fv_[l]->nsm[t], Ly.m, above); }
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
    // Vorkonditionierer: lokaler Teil je Dreieck (Zeile I: 1/2 X_L + 1/2 rho v; Zeile l: Jtop X_l - u_bot, skaliert mit sqrt|tau_0|)
    const int nb = static_cast<int>(8 * (L_ + 1));
    Pinv_.resize(N_);
    for (std::size_t t = 0; t < N_; ++t) {
        std::vector<cplx> A(nb * nb, cplx(0));
        auto put = [&](std::size_t row, std::size_t col, const Mat8* M, cplx s) {   // Block (row, col) += s M (M = nullptr: s I)
            for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
                A[(8 * row + i) * nb + 8 * col + j] += s * (M ? (*M)[i * 8 + j] : cplx(i == j ? 1.0 : 0.0));
        };
        const real rho = sq_[L_][t] / sq_[0][t];
        put(0, slot(L_), nullptr, 0.5); put(0, slot(0), nullptr, 0.5 * rho);
        for (std::size_t l = 1; l <= L_; ++l) {
            put(l, slot(l), &lay_[l - 1].Jtop[t], sq_[0][t] / sq_[l][t]);
            if (l == 1) put(l, slot(0), &Jb_[t], -1.0);
            else put(l, slot(l - 1), nullptr, -sq_[0][t] / sq_[l - 1][t]);
        }
        Pinv_[t] = inverse_n(A, nb);
    }
}

void TwoPortLayerProblem::resolvent(const SurfaceFV& fv, const std::vector<Multivector>& z, cplx shift, std::vector<Multivector>& x) const {
    const std::size_t N = N_;
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

void TwoPortLayerProblem::apply_g(std::size_t l, const std::vector<Multivector>& z, std::vector<Multivector>& gz) const {
    const std::size_t N = N_; const Layer& Ly = lay_[l]; const SurfaceFV& fv = *fv_[l];
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
    const std::size_t N = N_, L = L_;
    auto seg = [&](std::size_t j) { return x.begin() + 8 * N * slot(j); };
    std::vector<cplx> XL(seg(L), seg(L) + 8 * N), V(seg(0), seg(0) + 8 * N), E2X, E1V;
    E2_->apply(XL, E2X); E1_->apply(V, E1V);
    y.assign(size(), cplx(0));
    for (std::size_t i = 0; i < 8 * N; ++i) y[i] = 0.5 * (XL[i] + E2X[i]) + (sq_[L][i / 8] / sq_[0][i / 8]) * 0.5 * (V[i] - E1V[i]);   // Zeile I
    // Werte aller Spuren
    std::vector<std::vector<Multivector>> X(L + 1, std::vector<Multivector>(N));
    for (std::size_t j = 0; j <= L; ++j) {
        const auto b = x.begin() + 8 * N * slot(j);
        for (std::size_t t = 0; t < N; ++t) { const real s = 1 / sq_[j][t]; for (int c = 0; c < 8; ++c) X[j][t].c[c] = b[8 * t + c] * s; }
    }
    std::vector<Multivector> ut(N), ubot(N), z(N), DF, DSF, DDF, Bz(N), gBz;
    for (std::size_t l = 1; l <= L; ++l) {
        const Layer& Ly = lay_[l - 1];
        for (std::size_t t = 0; t < N; ++t) {
            ut[t] = apply8(Ly.Jtop[t], X[l][t]);
            ubot[t] = l == 1 ? apply8(Jb_[t], X[0][t]) : X[l - 1][t];
            z[t] = apply8(Ly.Jtops[t], X[l][t]) + (l == 1 ? apply8(Jbs_[t], X[0][t]) : X[l - 1][t]);
        }
        const SurfaceFV& fv = *fv_[l - 1];                                  // untere Flaeche der Schicht l
        fv.dirac_fit(z, DF, DSF, DDF);
        const Multivector iK = central(Ly.m.k(omega_, +1), Ly.m.k(omega_, -1)) * cplx(0, 1);
        for (std::size_t t = 0; t < N; ++t) {
            Multivector Xv = iK * z[t] - DF[t];
            if (opt_.curvature) Xv = Xv + DSF[t] * Ly.numid;
            Bz[t] = Multivector::vector(fv.nsm[t]) * Xv;                   // B in der Mitte der Schicht l
        }
        apply_g(l - 1, Bz, gBz);
        for (std::size_t t = 0; t < N; ++t) {
            const Multivector r = ut[t] - ubot[t] - gBz[t];
            for (int c = 0; c < 8; ++c) y[8 * N * l + 8 * t + c] = r.c[c] * sq_[0][t];
        }
    }
}

void TwoPortLayerProblem::precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    const std::size_t N = N_; const int nb = static_cast<int>(8 * (L_ + 1)); y.resize(x.size());
    std::vector<cplx> in(nb), out(nb);
    for (std::size_t t = 0; t < N; ++t) {
        for (std::size_t r = 0; r <= L_; ++r) for (int c = 0; c < 8; ++c) in[8 * r + c] = x[8 * N * r + 8 * t + c];   // Zeilen 0..L
        for (int r = 0; r < nb; ++r) { out[r] = 0; for (int c = 0; c < nb; ++c) out[r] += Pinv_[t][r * nb + c] * in[c]; }
        for (std::size_t s = 0; s <= L_; ++s) for (int c = 0; c < 8; ++c) y[8 * N * s + 8 * t + c] = out[8 * s + c];   // Positionen
    }
}

LayeredResult TwoPortLayerProblem::solve_plane_wave(const Vec3& dir, const CVec3& p, const SolveOptions& o) const {
    const cplx k = outer_.k(omega_);
    const TriangleMesh& ma = S_[L_];
    const std::vector<cplx> ba = project_plane_wave(ma, k, outer_.eps, dir, p);
    std::vector<cplx> b(size(), cplx(0)); std::copy(ba.begin(), ba.end(), b.begin());
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { precondition(x, y); };
    LayeredResult r; GmresResult g = gmres(A, b, r.h, &M, o.tol, o.restart, o.max_iter);
    r.iterations = g.iterations; r.residual = g.rel_residual;
    std::vector<cplx> hs(8 * N_); for (std::size_t i = 0; i < 8 * N_; ++i) hs[i] = r.h[8 * N_ * slot(L_) + i] - ba[i];
    r.sigma_ext = extinction_cross_section(ma, hs, k, outer_.eps, dir, p);
    r.forward = forward_amplitude(ma, hs, k, outer_.eps, dir, p);
    return r;
}

}  // namespace cbem
