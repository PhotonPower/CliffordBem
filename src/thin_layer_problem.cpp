#include "cbem/problems/thin_layer_problem.hpp"
#include <map>
#include <stdexcept>
#include "cbem/sources/fields.hpp"

namespace cbem {

namespace {
const int VEC[3] = {1, 2, 4};
const int BIV[3] = {6, 5, 3};      // I e1 = e23, I e2 = -e13, I e3 = e12
const real BS[3] = {1, -1, 1};
using C3 = std::array<cplx, 3>;
C3 cross_rc(const Vec3& a, const C3& b) { return {a.y * b[2] - a.z * b[1], a.z * b[0] - a.x * b[2], a.x * b[1] - a.y * b[0]}; }
cplx dot_rc(const Vec3& a, const C3& b) { return a.x * b[0] + a.y * b[1] + a.z * b[2]; }
}  // namespace

SurfaceFV::SurfaceFV(const TriangleMesh& m) : mesh(m) {
    const std::size_t N = m.size();
    std::map<std::pair<int, int>, std::vector<std::pair<std::size_t, int>>> em;
    for (std::size_t t = 0; t < N; ++t)
        for (int a = 0; a < 3; ++a) {
            int p = m.T[t][a], q = m.T[t][(a + 1) % 3];
            em[{std::min(p, q), std::max(p, q)}].push_back({t, a});
        }
    edges.resize(N);
    auto conormal = [&](std::size_t t, int a) {
        const Vec3 A = m.P[m.T[t][a]], B = m.P[m.T[t][(a + 1) % 3]];
        Vec3 c = cross(B - A, m.normal[t]); c = c / norm(c);
        if (dot(c, (A + B) * 0.5 - m.centroid[t]) < 0) c = c * -1.0;   // nach aussen (vom Dreieck weg)
        return c;
    };
    for (const auto& kv : em) {
        if (kv.second.size() != 2) throw std::runtime_error("SurfaceFV: Netz nicht geschlossen/konform (Kante mit " + std::to_string(kv.second.size()) + " Dreiecken)");
        const auto [t0, a0] = kv.second[0]; const auto [t1, a1] = kv.second[1];
        const real len = norm(m.P[kv.first.first] - m.P[kv.first.second]);
        const Vec3 c0 = conormal(t0, a0), c1 = conormal(t1, a1);
        edges[t0][a0] = Edge{t1, len, c0, c1, Vec3{}};
        edges[t1][a1] = Edge{t0, len, c1, c0, Vec3{}};
    }
    // glatte Normalen: winkelgewichtete Knotennormalen, am Schwerpunkt gemittelt
    {
        std::vector<Vec3> vn(m.P.size(), Vec3{});
        for (std::size_t t = 0; t < N; ++t)
            for (int a = 0; a < 3; ++a) {
                const Vec3 e1 = m.P[m.T[t][(a + 1) % 3]] - m.P[m.T[t][a]], e2 = m.P[m.T[t][(a + 2) % 3]] - m.P[m.T[t][a]];
                const real w = std::acos(std::max(-1.0, std::min(1.0, dot(e1, e2) / (norm(e1) * norm(e2)))));
                vn[m.T[t][a]] += m.normal[t] * w;
            }
        for (auto& v : vn) v = v / norm(v);
        nsm.resize(N);
        for (std::size_t t = 0; t < N; ++t) { Vec3 s = vn[m.T[t][0]] + vn[m.T[t][1]] + vn[m.T[t][2]]; nsm[t] = s / norm(s); }
    }
    // Kleinste-Quadrate-Gradient: g = sum_j a_j (phi_j - phi_t),  a_j = M^{-1} r_j,  r_j = P_t (c_j - c_t),  M = sum r r^T
    for (std::size_t t = 0; t < N; ++t) {
        const Vec3& n = m.normal[t];
        Vec3 u = cross(n, std::abs(n.x) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)); u = u / norm(u); const Vec3 v = cross(n, u);
        real M[2][2] = {{0, 0}, {0, 0}}; real r2[3][2];
        for (int a = 0; a < 3; ++a) {
            const Vec3 r = m.centroid[edges[t][a].nb] - m.centroid[t];
            r2[a][0] = dot(r, u); r2[a][1] = dot(r, v);
            for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) M[i][j] += r2[a][i] * r2[a][j];
        }
        const real det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
        if (!(std::abs(det) > 1e-14 * (M[0][0] + M[1][1]) * (M[0][0] + M[1][1]))) throw std::runtime_error("SurfaceFV: entartete Nachbarschaft");
        for (int a = 0; a < 3; ++a) {
            const real x = ( M[1][1] * r2[a][0] - M[0][1] * r2[a][1]) / det;
            const real y = (-M[1][0] * r2[a][0] + M[0][0] * r2[a][1]) / det;
            edges[t][a].lsq = u * x + v * y;
        }
    }
    // Formoperator: Kleinste-Quadrate-Gradient der Normalenkomponenten, symmetrisiert und tangential projiziert
    shape.assign(N, std::array<real, 9>{});
    for (std::size_t t = 0; t < N; ++t) {
        real M[3][3] = {};
        for (const Edge& e : edges[t]) {
            const Vec3 dn = nsm[e.nb] - nsm[t];
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) M[a][b] += e.lsq[a] * dn[b];
        }
        const Vec3& n = nsm[t]; real P[3][3], Sy[3][3], R[3][3];
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) { P[a][b] = (a == b) - n[a] * n[b]; Sy[a][b] = 0.5 * (M[a][b] + M[b][a]); }
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) { R[a][b] = 0; for (int c = 0; c < 3; ++c) R[a][b] += P[a][c] * Sy[c][b]; }
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) { real s = 0; for (int c = 0; c < 3; ++c) s += R[a][c] * P[c][b]; shape[t][3 * a + b] = s; }
    }
}

void SurfaceFV::dirac(const std::vector<Multivector>& F, std::vector<Multivector>& DF, std::vector<Multivector>* DSF) const {
    const std::size_t N = F.size();
    DF.assign(N, Multivector{}); if (DSF) DSF->assign(N, Multivector{});
    for (std::size_t t = 0; t < N; ++t) {
        std::array<CVec3, 8> g{};                                    // Gradient jeder Komponente
        for (const Edge& e : edges[t])
            for (int c = 0; c < 8; ++c) {
                const cplx w = F[e.nb].c[c] - F[t].c[c];
                g[c][0] += w * e.lsq.x; g[c][1] += w * e.lsq.y; g[c][2] += w * e.lsq.z;
            }
        const auto& S = shape[t];
        for (int c = 0; c < 8; ++c) {
            DF[t] = DF[t] + Multivector::vector(g[c]) * Multivector::blade(c);
            if (DSF) {
                CVec3 sg{};
                for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) sg[a] += S[3 * a + b] * g[c][b];
                (*DSF)[t] = (*DSF)[t] + Multivector::vector(sg) * Multivector::blade(c);
            }
        }
    }
}

void SurfaceFV::grad(const std::vector<cplx>& phi, std::vector<std::array<cplx, 3>>& g) const {
    g.assign(phi.size(), {cplx(0), cplx(0), cplx(0)});
    for (std::size_t t = 0; t < phi.size(); ++t)
        for (const Edge& e : edges[t]) {
            const cplx w = phi[e.nb] - phi[t];
            g[t][0] += w * e.lsq.x; g[t][1] += w * e.lsq.y; g[t][2] += w * e.lsq.z;
        }
}

void SurfaceFV::div(const std::vector<std::array<cplx, 3>>& v, std::vector<cplx>& dv) const {
    dv.assign(v.size(), cplx(0));
    for (std::size_t t = 0; t < v.size(); ++t) {
        for (const Edge& e : edges[t]) dv[t] += 0.5 * (dot_rc(e.conormal, v[t]) - dot_rc(e.nb_conormal, v[e.nb])) * e.len;
        dv[t] /= mesh.area[t];
    }
}

ThinLayerTransmissionOperator::ThinLayerTransmissionOperator(const TriangleMesh& m, const BoundaryOperator& E1, const BoundaryOperator& E2,
                                                             const Medium& in, const Medium& out, const std::vector<Coating>& coatings,
                                                             real f, real omega, ThinLayerModel model)
    : m_(m), E1_(E1), E2_(E2), in_(in), out_(out), fv_(m), omega_(omega), model_(model) {
    if (std::abs(in.chi) > 0 || std::abs(out.chi) > 0) throw std::invalid_argument("ThinLayer: nur achirale Medien");
    if (!(f >= 0.0 && f <= 1.0)) throw std::invalid_argument("ThinLayer: inner_fraction muss in [0, 1] liegen");
    for (const Coating& c : coatings) {
        if (std::abs(c.medium.chi) > 0) throw std::invalid_argument("ThinLayer: nur achirale Schichten");
        for (int side = 0; side < 2; ++side) {                  // Anteil innen (verdraengt in) und aussen (verdraengt out)
            const Medium& X = side == 0 ? in : out; const real d = c.thickness * (side == 0 ? f : 1.0 - f);
            aE_ += d * (1.0 / c.medium.eps - 1.0 / X.eps); bE_ += d * (c.medium.eps - X.eps);
            aH_ += d * (1.0 / c.medium.mu - 1.0 / X.mu);   bH_ += d * (c.medium.mu - X.mu);
        }
    }
    // Dirac-Form: Stapel von aussen nach innen, Grenzen von (1 - f) T bis -f T (T = Gesamtdicke)
    real T = 0; for (const Coating& c : coatings) T += c.thickness;
    stack_.push_back(out); bound_.push_back((1.0 - f) * T);
    for (auto it = coatings.rbegin(); it != coatings.rend(); ++it) { stack_.push_back(it->medium); bound_.push_back(bound_.back() - it->thickness); }
    stack_.push_back(in);
    auto same = [](const Medium& a, const Medium& b) { return a.eps == b.eps && a.mu == b.mu && a.chi == b.chi; };
    for (std::size_t l = 0; l + 1 < stack_.size();) {                        // gleiche Nachbarmedien: Grenze entfaellt
        if (same(stack_[l], stack_[l + 1])) { stack_.erase(stack_.begin() + l + 1); bound_.erase(bound_.begin() + l); }
        else ++l;
    }
    Jstep_.resize(stack_.size() - 1);
    for (std::size_t l = 0; l + 1 < stack_.size(); ++l) {
        Jstep_[l].resize(m.size());
        for (std::size_t t = 0; t < m.size(); ++t) Jstep_[l][t] = transmission_map(fv_.nsm[t], stack_[l + 1], stack_[l]);   // F_{l+1} aus F_l
    }
    Jsm_.resize(m.size());
    for (std::size_t t = 0; t < m.size(); ++t) Jsm_[t] = transmission_map(fv_.nsm[t], in, out);
    J_.resize(m.size()); P_.resize(m.size());
    for (std::size_t t = 0; t < m.size(); ++t) {
        J_[t] = transmission_map(m.normal[t], in, out);
        Mat8 A = J_[t]; for (int i = 0; i < 8; ++i) A[i * 8 + i] += 1.0;
        P_[t] = inverse8(A); for (auto& v : P_[t]) v *= 2.0;
    }
}

void ThinLayerTransmissionOperator::propagate(std::vector<Multivector>& F, const Medium& md, real nu0, real nu1) const {
    const real s = nu1 - nu0; if (s == 0.0) return;
    const std::size_t N = F.size(); const cplx ik = cplx(0, 1) * md.k(omega_);
    const bool second = model_ == ThinLayerModel::Dirac2;
    std::vector<Multivector> DF, DSF, G1(N), DG, G2;
    fv_.dirac(F, DF, second ? &DSF : nullptr);
    for (std::size_t t = 0; t < N; ++t) G1[t] = Multivector::vector(fv_.nsm[t]) * (F[t] * ik - DF[t]);     // B F
    if (second) fv_.dirac(G1, DG);
    for (std::size_t t = 0; t < N; ++t) {
        Multivector R = F[t] + G1[t] * s;
        if (second) {
            const Multivector n = Multivector::vector(fv_.nsm[t]);
            const Multivector BB = n * (G1[t] * ik - DG[t]);                                                  // B^2 F
            R = R + BB * (0.5 * s * s) + (n * DSF[t]) * (s * nu0 + 0.5 * s * s);                              // + B' F
        }
        F[t] = R;
    }
}

void ThinLayerTransmissionOperator::apply_dirac(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    const std::size_t N = m_.size();
    std::vector<Multivector> F(N);
    for (std::size_t t = 0; t < N; ++t) { const real s = 1.0 / std::sqrt(m_.area[t]); for (int c = 0; c < 8; ++c) F[t].c[c] = x[8 * t + c] * s; }
    const std::vector<Multivector> F0 = F;
    real nu = 0.0;
    for (std::size_t l = 0; l < stack_.size(); ++l) {
        const real target = l + 1 < stack_.size() ? bound_[l] : 0.0;   // naechste Grenze bzw. zurueck zu nu = 0 in Medium 1
        propagate(F, stack_[l], nu, target); nu = target;
        if (l + 1 < stack_.size())
            for (std::size_t t = 0; t < N; ++t) { Multivector G; cbem::apply(Jstep_[l][t], F[t].c.data(), G.c.data()); F[t] = G; }
    }
    y.resize(8 * N);
    for (std::size_t t = 0; t < N; ++t) {                                   // J(n) x + sqrt|tau| (Kette - J(nsm)) F0
        cplx a[8], b[8]; cbem::apply(J_[t], &x[8 * t], a); cbem::apply(Jsm_[t], F0[t].c.data(), b);
        const real s = std::sqrt(m_.area[t]);
        for (int c = 0; c < 8; ++c) y[8 * t + c] = a[c] + (F[t].c[c] - b[c]) * s;
    }
}

void ThinLayerTransmissionOperator::apply_Jeff(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    if (model_ != ThinLayerModel::Jump1) { apply_dirac(x, y); return; }
    const std::size_t N = m_.size();
    const cplx se2 = std::sqrt(out_.eps), sm2 = std::sqrt(out_.mu), se1 = std::sqrt(in_.eps), sm1 = std::sqrt(in_.mu);
    std::vector<cplx> phiE(N), phiH(N), divE, divH;
    std::vector<C3> Et(N), Ht(N), gE, gH;
    for (std::size_t t = 0; t < N; ++t) {
        const real s = 1.0 / std::sqrt(m_.area[t]); const Vec3& n = fv_.nsm[t];
        C3 e, h;
        for (int d = 0; d < 3; ++d) { e[d] = x[8 * t + VEC[d]] * s / se2; h[d] = x[8 * t + BIV[d]] * BS[d] * s / sm2; }
        const cplx en = dot_rc(n, e), hn = dot_rc(n, h);
        phiE[t] = out_.eps * en; phiH[t] = out_.mu * hn;                      // D_n, B_n (Aussenseite)
        for (int d = 0; d < 3; ++d) { Et[t][d] = e[d] - en * n[d]; Ht[t][d] = h[d] - hn * n[d]; }
    }
    fv_.grad(phiE, gE); fv_.grad(phiH, gH); fv_.div(Et, divE); fv_.div(Ht, divH);
    y.assign(8 * N, cplx(0));
    const cplx iw(0, omega_);
    for (std::size_t t = 0; t < N; ++t) {
        cbem::apply(J_[t], &x[8 * t], &y[8 * t]);
        const Vec3& n = fv_.nsm[t]; const real sa = std::sqrt(m_.area[t]);
        const C3 nxH = cross_rc(n, Ht[t]), nxE = cross_rc(n, Et[t]);
        const cplx jDn = -bE_ * divE[t], jBn = -bH_ * divH[t];
        for (int d = 0; d < 3; ++d) {
            const cplx jEt = aE_ * gE[t][d] - iw * bH_ * nxH[d];
            const cplx jHt = aH_ * gH[t][d] + iw * bE_ * nxE[d];
            const cplx dv = -se1 * (jEt + jDn / in_.eps * n[d]);              // Delta der Innenspur (Werte)
            const cplx dh = -sm1 * (jHt + jBn / in_.mu * n[d]);
            y[8 * t + VEC[d]] += sa * dv;
            y[8 * t + BIV[d]] += sa * dh * BS[d];
        }
    }
}

void ThinLayerTransmissionOperator::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    std::vector<cplx> Jx, E1Jx, E2x;
    apply_Jeff(x, Jx);
    E2_.apply(x, E2x); E1_.apply(Jx, E1Jx);
    y.resize(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) y[i] = 0.5 * (x[i] + E2x[i]) + 0.5 * (Jx[i] - E1Jx[i]);
}

void ThinLayerTransmissionOperator::precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    y.resize(x.size());
    for (std::size_t t = 0; t < m_.size(); ++t) cbem::apply(P_[t], &x[8 * t], &y[8 * t]);
}

ThinLayerScatteringProblem::ThinLayerScatteringProblem(const TriangleMesh& surface, const Medium& core, const std::vector<Coating>& coatings,
                                                       real omega, Medium outer, real inner_fraction, HMatrixParams hp, EntryParams ep,
                                                       ThinLayerModel model)
    : m_(surface), core_(core), outer_(outer), omega_(omega) {
    if (m_.normal.size() != m_.size()) m_.compute_geometry();
    K1_ = std::make_unique<KernelEntries>(m_, core.k(omega), ep); K2_ = std::make_unique<KernelEntries>(m_, outer.k(omega), ep);
    H1_ = std::make_unique<KernelHMatrix>(*K1_, hp); H2_ = std::make_unique<KernelHMatrix>(*K2_, hp);
    E1_ = std::make_unique<CauchyOperator>(m_, *H1_); E2_ = std::make_unique<CauchyOperator>(m_, *H2_);
    T_ = std::make_unique<ThinLayerTransmissionOperator>(m_, *E1_, *E2_, core, outer, coatings, inner_fraction, omega, model);
}

LayeredResult ThinLayerScatteringProblem::solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o) const {
    const cplx k = outer_.k(omega_);
    const std::vector<cplx> b = project_plane_wave(m_, k, outer_.eps, d, p);
    LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->apply(x, y); };
    LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T_->precondition(x, y); };
    LayeredResult r; GmresResult g = gmres(A, b, r.h, &M, o.tol, o.restart, o.max_iter);
    r.iterations = g.iterations; r.residual = g.rel_residual;
    std::vector<cplx> hs(b.size()); for (std::size_t i = 0; i < b.size(); ++i) hs[i] = r.h[i] - b[i];
    r.sigma_ext = extinction_cross_section(m_, hs, k, outer_.eps, d, p);
    r.forward = forward_amplitude(m_, hs, k, outer_.eps, d, p);
    return r;
}

}  // namespace cbem
