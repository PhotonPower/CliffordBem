#include "cbem/sources/near_field.hpp"
#include <chrono>
#include <cmath>
#include "cbem/geometry/quadrature.hpp"
#include "cbem/kernel/dirac_kernel.hpp"
#include "cbem/kernel/triangle_integrals.hpp"

namespace cbem {

namespace {
// Vorzeichen der Darstellung des Aussenfeldes durch das Cauchy-Integral: +1, konsistent mit far_field (Fernfeldgrenze)
// und mit der Randwertgrenze F_s -> h_s; gegen die Mie-Loesung in tests/test_near_field.cpp
constexpr real kSign = +1.0;
Multivector projector(int s) { return Multivector::blade(0, 0.5) + Multivector::blade(7, cplx(0, 0.5 * s)); }
}  // namespace

namespace {
// Eintraege K(x, t) = int_tau Phi_k(x - y) dS_y (Skalar, Vektor): fern 7-Punkt-Regel, nah analytisch fuer die singulaeren Anteile
struct EntryEval {
    const TriangleMesh& m; cplx k; MeshQuadrature q; std::vector<real> rad;
    EntryEval(const TriangleMesh& mesh, cplx kk) : m(mesh), k(kk), q(mesh, QuadRule::dunavant7()), rad(mesh.size()) {
        for (std::size_t t = 0; t < m.size(); ++t) { const auto tri = m.vertices(t); for (const Vec3& v : tri) rad[t] = std::max(rad[t], norm(v - m.centroid[t])); }
    }
    std::array<cplx, 4> operator()(const Vec3& x, std::size_t t) const {
        cplx S = 0; CVec3 V{};
        const Vec3* qp = q.points(t); const real* qw = q.weights(t);
        if (norm(x - m.centroid[t]) > 4 * rad[t]) {
            for (int a = 0; a < q.q; ++a) {
                const Vec3 z = x - qp[a]; const KernelValue kv = dirac_kernel_full(z, k);
                S += qw[a] * kv.s; const cplx c = qw[a] * kv.vcoef; V[0] += c * z.x; V[1] += c * z.y; V[2] += c * z.z;
            }
        } else {
            Vec3 Ig; real Ii; triangle_integrals(x, m.vertices(t), m.normal[t], Ig, Ii);
            S += -cplx(0, 1) * k * Ii / (4 * pi); V[0] += Ig.x / (4 * pi); V[1] += Ig.y / (4 * pi); V[2] += Ig.z / (4 * pi);
            for (int a = 0; a < q.q; ++a) {
                const Vec3 z = x - qp[a]; const KernelValue kv = dirac_kernel_remainder(z, k);
                S += qw[a] * kv.s; const cplx c = qw[a] * kv.vcoef; V[0] += c * z.x; V[1] += c * z.y; V[2] += c * z.z;
            }
        }
        return {S, V[0], V[1], V[2]};
    }
};
// n_tau u_tau je Dreieck
std::vector<Multivector> densities(const TriangleMesh& m, const std::vector<cplx>& hs) {
    std::vector<Multivector> G(m.size());
    for (std::size_t t = 0; t < m.size(); ++t) {
        Multivector u; for (int c = 0; c < 8; ++c) u.c[c] = hs[8 * t + c] / std::sqrt(m.area[t]);
        G[t] = Multivector::vector(m.normal[t]) * u;
    }
    return G;
}
}  // namespace

std::vector<Multivector> scattered_field(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, const std::vector<Vec3>& pts) {
    const EntryEval ev(m, k); const std::vector<Multivector> G = densities(m, hs);
    std::vector<Multivector> out(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        Multivector F;
        for (std::size_t t = 0; t < m.size(); ++t) {
            const auto K = ev(pts[i], t);
            F = F + (Multivector::blade(0, K[0]) + Multivector::vector(CVec3{K[1], K[2], K[3]})) * G[t];
        }
        out[i] = F * kSign;
    }
    return out;
}

NearFieldOperator::NearFieldOperator(const TriangleMesh& m, cplx k, const std::vector<Vec3>& pts, HMatrixParams prm)
    : m_(m), k_(k), prm_(prm), rows_(pts, std::vector<real>(pts.size(), 0.0), prm.leaf), cols_(m, prm.leaf) {
    auto t0 = std::chrono::steady_clock::now();
    const EntryEval ev(m, k);
    std::vector<std::pair<int, int>> adm, inadm;
    partition(0, 0, adm, inadm);
    dense_.resize(inadm.size()); lr_.resize(adm.size());
    for (std::size_t b = 0; b < inadm.size(); ++b) {
        Dense& D = dense_[b];
        D.R = rows_.indices(rows_.nodes[inadm[b].first]); D.C = cols_.indices(cols_.nodes[inadm[b].second]);
        D.K.resize(D.R.size() * D.C.size());
        for (std::size_t a = 0; a < D.R.size(); ++a) for (std::size_t c = 0; c < D.C.size(); ++c) D.K[a * D.C.size() + c] = ev(pts[D.R[a]], D.C[c]);
    }
    for (std::size_t b = 0; b < adm.size(); ++b) {
        LR& B = lr_[b];
        B.R = rows_.indices(rows_.nodes[adm[b].first]); B.C = cols_.indices(cols_.nodes[adm[b].second]);
        const std::size_t mm = B.R.size(), n = B.C.size();
        RowFn row = [&](std::size_t i, cplx* out) { for (std::size_t j = 0; j < n; ++j) { const auto K = ev(pts[B.R[i]], B.C[j]); for (int c = 0; c < 4; ++c) out[c * n + j] = K[c]; } };
        ColFn col = [&](std::size_t J, cplx* out) { const std::size_t c = J / n, j = J % n; for (std::size_t i = 0; i < mm; ++i) out[i] = ev(pts[B.R[i]], B.C[j])[c]; };
        B.f = aca_select(prm_.aca_plus, row, col, mm, 4 * n, prm_.eps); recompress(B.f, prm_.eps);
    }
    st_.n_dense = dense_.size(); st_.n_lowrank = lr_.size();
    std::size_t rsum = 0;
    for (auto& D : dense_) st_.entries_dense += 4 * D.K.size();
    for (auto& B : lr_) { st_.entries_lowrank += B.f.storage(); rsum += B.f.rank(); st_.max_rank = std::max(st_.max_rank, B.f.rank()); }
    st_.mean_rank = lr_.empty() ? 0 : double(rsum) / lr_.size();
    st_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void NearFieldOperator::partition(int t, int s, std::vector<std::pair<int, int>>& adm, std::vector<std::pair<int, int>>& inadm) const {
    const ClusterNode& a = rows_.nodes[t]; const ClusterNode& b = cols_.nodes[s];
    const real dist = box_distance(a, b);
    const bool ok = std::min(a.diam, b.diam) <= prm_.eta * dist && dist > 0 && dist > prm_.sep_factor * b.hmax &&
                    std::abs(k_) * std::max(a.diam, b.diam) <= prm_.max_kdiam;
    if (ok) { adm.emplace_back(t, s); return; }
    if (a.leaf() && b.leaf()) { inadm.emplace_back(t, s); return; }
    // in den groesseren (nicht Blatt-)Cluster absteigen
    if (b.leaf() || (!a.leaf() && a.diam >= b.diam)) { for (int x : a.child) partition(x, s, adm, inadm); }
    else { for (int y : b.child) partition(t, y, adm, inadm); }
}

std::vector<Multivector> NearFieldOperator::apply(const std::vector<cplx>& hs) const {
    const std::vector<Multivector> G = densities(m_, hs);
    const std::size_t N = m_.size();
    std::vector<cplx> Z(32 * N);                                        // Z_{t,0} = G_t, Z_{t,a} = e_a G_t
    for (std::size_t t = 0; t < N; ++t) {
        for (int q = 0; q < 8; ++q) Z[t * 32 + q] = G[t].c[q];
        const Vec3 e[3] = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
        for (int a = 0; a < 3; ++a) { const Multivector g = Multivector::vector(e[a]) * G[t]; for (int q = 0; q < 8; ++q) Z[t * 32 + (a + 1) * 8 + q] = g.c[q]; }
    }
    const std::size_t M = rows_.perm.size();
    std::vector<cplx> Y(8 * M, cplx(0));
    for (const auto& D : dense_) {
        const std::size_t n = D.C.size();
        for (std::size_t a = 0; a < D.R.size(); ++a) {
            cplx* y = &Y[D.R[a] * 8];
            for (std::size_t c = 0; c < n; ++c) { const Comp& K = D.K[a * n + c]; const cplx* z = &Z[D.C[c] * 32];
                for (int q = 0; q < 4; ++q) for (int r = 0; r < 8; ++r) y[r] += K[q] * z[q * 8 + r]; }
        }
    }
    std::vector<cplx> tmp;
    for (const auto& B : lr_) {
        const std::size_t m = B.R.size(), n = B.C.size(), r = B.f.rank();
        tmp.assign(r * 8, cplx(0));
        for (std::size_t k = 0; k < r; ++k) {
            const cplx* v = B.f.V.col(k);
            for (int c = 0; c < 4; ++c) for (std::size_t j = 0; j < n; ++j) { const cplx* z = &Z[B.C[j] * 32 + c * 8]; const cplx vv = v[c * n + j];
                for (int q = 0; q < 8; ++q) tmp[k * 8 + q] += vv * z[q]; }
        }
        for (std::size_t a = 0; a < m; ++a) { cplx* y = &Y[B.R[a] * 8]; for (std::size_t k = 0; k < r; ++k) { const cplx u = B.f.U(a, k); for (int q = 0; q < 8; ++q) y[q] += u * tmp[k * 8 + q]; } }
    }
    std::vector<Multivector> out(M);
    for (std::size_t i = 0; i < M; ++i) { for (int q = 0; q < 8; ++q) out[i].c[q] = Y[8 * i + q]; out[i] = out[i] * kSign; }
    return out;
}

std::vector<NearFieldPoint> exterior_near_field(const TriangleMesh& outer0, const std::vector<cplx>& h, const Medium& m, real omega,
                                                const Vec3& d0, const CVec3& p, const std::vector<Vec3>& pts, const NearFieldOptions& opt) {
    TriangleMesh tmp;                                                   // Geometrie (auch Elementgroessen) bei Bedarf berechnen
    const TriangleMesh& outer = (outer0.hmax.size() == outer0.size() && outer0.centroid.size() == outer0.size()) ? outer0 : (tmp = outer0, tmp.compute_geometry(), tmp);
    const Vec3 d = d0 / norm(d0);
    const PlaneWaveIncidence inc = plane_wave_incidence(m, omega, d, p);
    const std::vector<cplx> b = project_plane_wave(outer, inc.k, m.eps, d, p);
    std::vector<cplx> hs(h.size()); for (std::size_t i = 0; i < h.size(); ++i) hs[i] = h[i] - b[i];
    const bool useH = pts.size() >= opt.hmatrix_min_points;
    auto field = [&](const std::vector<cplx>& hh, cplx kk) {
        if (!useH) return scattered_field(outer, hh, kk, pts);
        HMatrixParams hp; hp.eps = opt.eps; return NearFieldOperator(outer, kk, pts, hp).apply(hh);
    };
    std::vector<Multivector> Fs;
    if (std::abs(m.chi) > 0) {                                          // je Helizitaet mit k_pm
        Fs.assign(pts.size(), Multivector{});
        for (int s : {+1, -1}) {
            const auto part = field(helicity_part(hs, s), m.k(omega, s));
            const Multivector P = projector(s);
            for (std::size_t i = 0; i < pts.size(); ++i) Fs[i] = Fs[i] + P * part[i];
        }
    } else Fs = field(hs, inc.k);
    const cplx se = std::sqrt(m.eps), sm = std::sqrt(m.mu);
    const CVec3 dxp{d.y * p[2] - d.z * p[1], d.z * p[0] - d.x * p[2], d.x * p[1] - d.y * p[0]};
    const real p2 = std::norm(p[0]) + std::norm(p[1]) + std::norm(p[2]);
    const real C0 = std::abs(std::real(se / sm)) * p2;                  // |Im(conj(E0).H0)| der zirkularen Welle, H0 = sqrt(eps/mu) d x E0
    std::vector<NearFieldPoint> out(pts.size());
    real hm = 0; for (real hh : outer.hmax) hm += hh; hm /= std::max<std::size_t>(1, outer.hmax.size());
    const std::vector<real> dist = distance_to_surface(outer, pts, 0.02 * hm);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        NearFieldPoint& r = out[i];
        r.inside = winding_number(outer, pts[i]) > 0.5;
        r.too_close = dist[i] < 0.02 * hm;
        const cplx ph = std::exp(cplx(0, 1) * inc.k * dot(d, pts[i]));
        const int VEC[3] = {1, 2, 4}, BIV[3] = {6, 5, 3}; const real BS[3] = {1, -1, 1};
        for (int a = 0; a < 3; ++a) {
            r.E[a] = Fs[i].c[VEC[a]] / se + p[a] * ph;
            r.H[a] = Fs[i].c[BIV[a]] * BS[a] / sm + (se / sm) * dxp[a] * ph;
        }
        real e2 = 0; cplx eh = 0;
        for (int a = 0; a < 3; ++a) { e2 += std::norm(r.E[a]); eh += std::conj(r.E[a]) * r.H[a]; }
        r.enhancement = e2 / p2;
        r.chirality = std::imag(eh) / C0;
    }
    return out;
}

}  // namespace cbem
