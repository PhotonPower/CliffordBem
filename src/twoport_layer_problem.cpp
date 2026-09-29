#include "cbem/problems/twoport_layer_problem.hpp"
#include <cmath>
#include <stdexcept>
#include "cbem/sources/fields.hpp"

namespace cbem {

namespace {
// Inverse einer komplexen 16x16-Matrix (Gauss-Jordan mit Pivotsuche), zeilenweise
std::array<cplx, 256> inverse16(std::array<cplx, 256> A) {
    std::array<cplx, 256> I{}; for (int i = 0; i < 16; ++i) I[i * 16 + i] = 1.0;
    for (int c = 0; c < 16; ++c) {
        int piv = c; for (int r = c + 1; r < 16; ++r) if (std::abs(A[r * 16 + c]) > std::abs(A[piv * 16 + c])) piv = r;
        if (std::abs(A[piv * 16 + c]) < 1e-300) throw std::runtime_error("TwoPort: singulaerer Vorkonditionierer");
        for (int q = 0; q < 16; ++q) { std::swap(A[c * 16 + q], A[piv * 16 + q]); std::swap(I[c * 16 + q], I[piv * 16 + q]); }
        const cplx dd = A[c * 16 + c]; for (int q = 0; q < 16; ++q) { A[c * 16 + q] /= dd; I[c * 16 + q] /= dd; }
        for (int r = 0; r < 16; ++r) if (r != c) { const cplx f = A[r * 16 + c]; if (f == cplx(0)) continue;
            for (int q = 0; q < 16; ++q) { A[r * 16 + q] -= f * A[c * 16 + q]; I[r * 16 + q] -= f * I[c * 16 + q]; } }
    }
    return I;
}
Multivector apply8(const Mat8& M, const Multivector& x) { Multivector y; cbem::apply(M, x.c.data(), y.c.data()); return y; }
Multivector central(cplx a, cplx b) { return Multivector::blade(0, 0.5 * (a + b)) + Multivector::blade(7, cplx(0, 0.5) * (a - b)); }  // a P+ + b P-
}  // namespace

TwoPortLayerProblem::TwoPortLayerProblem(const TriangleMesh& core_surface, const TriangleMesh& outer_surface, const Medium& core,
                                         const Medium& layer, real d, real omega, Medium outer, HMatrixParams hp, EntryParams ep,
                                         TwoPortOptions opt)
    : mb_(core_surface), ma_(outer_surface), core_(core), layer_(layer), outer_(outer), d_(d), omega_(omega), opt_(opt) {
    if (mb_.normal.size() != mb_.size()) mb_.compute_geometry();
    if (ma_.normal.size() != ma_.size()) ma_.compute_geometry();
    if (mb_.size() != ma_.size() || mb_.P.size() != ma_.P.size()) throw std::invalid_argument("TwoPort: Flaechen mit verschiedener Topologie");
    if (!(d > 0)) throw std::invalid_argument("TwoPort: Schichtdicke muss positiv sein");
    if (std::abs(outer.chi) > 0) throw std::invalid_argument("TwoPort: chirales Aussenmedium nicht unterstuetzt");
    N_ = mb_.size();
    fv_ = std::make_unique<SurfaceFV>(mb_);
    K2_ = std::make_unique<KernelEntries>(ma_, outer.k(omega), ep); H2_ = std::make_unique<KernelHMatrix>(*K2_, hp);
    E2_ = std::make_unique<CauchyOperator>(ma_, *H2_);
    if (std::abs(core.chi) > 0) {
        K1p_ = std::make_unique<KernelEntries>(mb_, core.k(omega, +1), ep); H1p_ = std::make_unique<KernelHMatrix>(*K1p_, hp);
        K1m_ = std::make_unique<KernelEntries>(mb_, core.k(omega, -1), ep); H1m_ = std::make_unique<KernelHMatrix>(*K1m_, hp);
        E1p_ = std::make_unique<CauchyOperator>(mb_, *H1p_); E1m_ = std::make_unique<CauchyOperator>(mb_, *H1m_);
        E1ch_ = std::make_unique<ChiralCauchyOperator>(*E1p_, *E1m_); E1_ = E1ch_.get();
    } else {
        K1p_ = std::make_unique<KernelEntries>(mb_, core.k(omega), ep); H1p_ = std::make_unique<KernelHMatrix>(*K1p_, hp);
        E1p_ = std::make_unique<CauchyOperator>(mb_, *H1p_); E1_ = E1p_.get();
    }
    Jc2_.resize(N_); Jc1_.resize(N_); Jc2s_.resize(N_); Jc1s_.resize(N_); rho_.resize(N_); Pinv_.resize(N_);
    for (std::size_t t = 0; t < N_; ++t) {
        Jc2_[t] = transmission_map(ma_.normal[t], layer, outer); Jc1_[t] = transmission_map(mb_.normal[t], layer, core);
        Jc2s_[t] = transmission_map(fv_->nsm[t], layer, outer); Jc1s_[t] = transmission_map(fv_->nsm[t], layer, core);
        rho_[t] = std::sqrt(ma_.area[t] / mb_.area[t]);
        std::array<cplx, 256> A{};
        for (int i = 0; i < 8; ++i) { A[i * 16 + i] = 0.5; A[i * 16 + 8 + i] = 0.5 * rho_[t]; }
        for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) { A[(8 + i) * 16 + j] = Jc2_[t][i * 8 + j] / rho_[t]; A[(8 + i) * 16 + 8 + j] = -Jc1_[t][i * 8 + j]; }
        Pinv_[t] = inverse16(A);
    }
    // Pole von g: sigma_m = ((2m+1) pi / d)^2 bis sigma_M >= 10 s_max (Gershgorin-Schranke fuer -Delta), Rest linear in s
    real smax = 0;
    for (std::size_t t = 0; t < N_; ++t) { real s = 0; for (real w : fv_->ring[t].lw) s += std::abs(w); smax = std::max(smax, 2 * s); }
    smax += std::max(std::norm(layer.k(omega, +1)), std::norm(layer.k(omega, -1)));
    if (opt.poles >= 0) M_ = opt.poles;
    else { const real q = d * std::sqrt(10 * smax) / pi; M_ = std::max(0, static_cast<int>(std::ceil((q - 1) / 2))); M_ = std::min(M_, 200); }
    real s2 = pi * pi / 8, s4 = std::pow(pi, 4) / 96;
    for (int m = 0; m < M_; ++m) { const real o = 2 * m + 1; sigma_.push_back(std::pow(o * pi / d, 2)); s2 -= 1 / (o * o); s4 -= 1 / std::pow(o, 4); }
    tail0_ = (4 / d) * (d * d / (pi * pi)) * s2;
    tail1_ = (4 / d) * std::pow(d / pi, 4) * s4;
}

void TwoPortLayerProblem::resolvent(const std::vector<Multivector>& z, cplx shift, std::vector<Multivector>& x) const {
    const std::size_t N = N_;
    const auto& R = fv_->ring;
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

void TwoPortLayerProblem::apply_g(const std::vector<Multivector>& z, std::vector<Multivector>& gz) const {
    const std::size_t N = N_;
    const cplx kp = layer_.k(omega_, +1), km = layer_.k(omega_, -1);
    const bool chiral = std::abs(layer_.chi) > 0;
    const Multivector Pp = central(1.0, 0.0), Pm = central(0.0, 1.0), K2 = central(kp * kp, km * km);
    // Rest: tail0 z - tail1 s z mit s z = -K^2 z - Delta z
    const auto& R = fv_->ring;
    gz.resize(N);
    for (std::size_t t = 0; t < N; ++t) {
        Multivector lap;
        for (std::size_t j = 0; j < R[t].nb.size(); ++j) lap = lap + (z[R[t].nb[j]] - z[t]) * R[t].lw[j];
        const Multivector sz = (K2 * z[t]) * (-1.0) - lap;
        gz[t] = z[t] * tail0_ - sz * tail1_;
    }
    std::vector<Multivector> zp(N), zm(N), x;
    if (chiral) for (std::size_t t = 0; t < N; ++t) { zp[t] = Pp * z[t]; zm[t] = Pm * z[t]; }
    for (int m = 0; m < M_; ++m) {
        const real w = 4 / d_;
        if (!chiral) {
            resolvent(z, sigma_[m] - kp * kp, x);
            for (std::size_t t = 0; t < N; ++t) gz[t] = gz[t] + x[t] * w;
        } else {
            resolvent(zp, sigma_[m] - kp * kp, x); for (std::size_t t = 0; t < N; ++t) gz[t] = gz[t] + x[t] * w;
            resolvent(zm, sigma_[m] - km * km, x); for (std::size_t t = 0; t < N; ++t) gz[t] = gz[t] + x[t] * w;
        }
    }
}

void TwoPortLayerProblem::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    const std::size_t N = N_;
    std::vector<cplx> H(x.begin(), x.begin() + 8 * N), V(x.begin() + 8 * N, x.end()), E2H, E1V;
    E2_->apply(H, E2H); E1_->apply(V, E1V);
    y.assign(16 * N, cplx(0));
    for (std::size_t i = 0; i < 8 * N; ++i) y[i] = 0.5 * (H[i] + E2H[i]) + rho_[i / 8] * 0.5 * (V[i] - E1V[i]);   // Zeile I
    // Zeile II (Werte auf Gamma_b)
    std::vector<Multivector> ua(N), ub(N), z(N);
    for (std::size_t t = 0; t < N; ++t) {
        Multivector ha, vb;
        const real sa = 1 / std::sqrt(ma_.area[t]), sb = 1 / std::sqrt(mb_.area[t]);
        for (int c = 0; c < 8; ++c) { ha.c[c] = H[8 * t + c] * sa; vb.c[c] = V[8 * t + c] * sb; }
        ua[t] = apply8(Jc2_[t], ha); ub[t] = apply8(Jc1_[t], vb);
        z[t] = apply8(Jc2s_[t], ha) + apply8(Jc1s_[t], vb);
    }
    std::vector<Multivector> DF, DSF, DDF, Bz(N), gBz;
    fv_->dirac_fit(z, DF, DSF, DDF);
    const Multivector iK = central(layer_.k(omega_, +1), layer_.k(omega_, -1)) * cplx(0, 1);
    for (std::size_t t = 0; t < N; ++t) {
        Multivector X = iK * z[t] - DF[t];
        if (opt_.curvature) X = X + DSF[t] * (0.5 * d_);
        Bz[t] = Multivector::vector(fv_->nsm[t]) * X;                     // B in der Schichtmitte
    }
    apply_g(Bz, gBz);
    for (std::size_t t = 0; t < N; ++t) {
        const real sb = std::sqrt(mb_.area[t]);
        const Multivector r = ua[t] - ub[t] - gBz[t];
        for (int c = 0; c < 8; ++c) y[8 * N + 8 * t + c] = r.c[c] * sb;
    }
}

void TwoPortLayerProblem::precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    const std::size_t N = N_; y.resize(x.size());
    for (std::size_t t = 0; t < N; ++t) {
        cplx in[16], out[16];
        for (int c = 0; c < 8; ++c) { in[c] = x[8 * t + c]; in[8 + c] = x[8 * N + 8 * t + c]; }
        for (int r = 0; r < 16; ++r) { out[r] = 0; for (int c = 0; c < 16; ++c) out[r] += Pinv_[t][r * 16 + c] * in[c]; }
        for (int c = 0; c < 8; ++c) { y[8 * t + c] = out[c]; y[8 * N + 8 * t + c] = out[8 + c]; }
    }
}

LayeredResult TwoPortLayerProblem::solve_plane_wave(const Vec3& dir, const CVec3& p, const SolveOptions& o) const {
    const cplx k = outer_.k(omega_);
    const std::vector<cplx> ba = project_plane_wave(ma_, k, outer_.eps, dir, p);
    std::vector<cplx> b(16 * N_, cplx(0)); std::copy(ba.begin(), ba.end(), b.begin());
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { precondition(x, y); };
    LayeredResult r; GmresResult g = gmres(A, b, r.h, &M, o.tol, o.restart, o.max_iter);
    r.iterations = g.iterations; r.residual = g.rel_residual;
    std::vector<cplx> hs(8 * N_); for (std::size_t i = 0; i < 8 * N_; ++i) hs[i] = r.h[i] - ba[i];
    r.sigma_ext = extinction_cross_section(ma_, hs, k, outer_.eps, dir, p);
    r.forward = forward_amplitude(ma_, hs, k, outer_.eps, dir, p);
    return r;
}

}  // namespace cbem
