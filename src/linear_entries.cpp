#include "cbem/assembly/linear_entries.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>

#include "cbem/kernel/dirac_kernel.hpp"
#include "cbem/kernel/triangle_integrals.hpp"

namespace cbem {

namespace {

// K[a*3+b] += w la[a] lb[b] (s, v z)
inline void acc(LinearBlock& K, const real la[3], const real lb[3], real w, cplx s, cplx vc, const Vec3& z) {
    for (int a = 0; a < 3; ++a) {
        const real wa = w * la[a];
        for (int b = 0; b < 3; ++b) {
            const real ww = wa * lb[b];
            KernelComp& k = K[a * 3 + b];
            k[0] += ww * s; k[1] += ww * vc * z.x; k[2] += ww * vc * z.y; k[3] += ww * vc * z.z;
        }
    }
}

LinearBlock zero_block() {
    LinearBlock K;
    for (auto& k : K) k = KernelComp{0, 0, 0, 0};
    return K;
}

real point_triangle_distance(const Vec3& x, const std::array<Vec3, 3>& p) {
    Vec3 n = cross(p[1] - p[0], p[2] - p[0]); real a2 = norm(n); n = n / a2;
    real w = dot(x - p[0], n); Vec3 q = x - n * w;
    auto inside = [&](const Vec3& a, const Vec3& b) { return dot(cross(b - a, q - a), n) >= 0; };
    if (inside(p[0], p[1]) && inside(p[1], p[2]) && inside(p[2], p[0])) return std::abs(w);
    real d = 1e300;
    for (int e = 0; e < 3; ++e) {
        const Vec3& a = p[e]; const Vec3& b = p[(e + 1) % 3]; Vec3 ab = b - a;
        real t = std::max(0.0, std::min(1.0, dot(x - a, ab) / dot(ab, ab)));
        d = std::min(d, norm(x - (a + ab * t)));
    }
    return d;
}
real point_boundary_distance(const Vec3& x, const std::array<Vec3, 3>& p) {
    real d = 1e300;
    for (int e = 0; e < 3; ++e) {
        const Vec3& a = p[e]; const Vec3& b = p[(e + 1) % 3]; Vec3 ab = b - a;
        real t = std::max(0.0, std::min(1.0, dot(x - a, ab) / dot(ab, ab)));
        d = std::min(d, norm(x - (a + ab * t)));
    }
    return d;
}

}  // namespace

LinearKernelEntries::LinearKernelEntries(const TriangleMesh& mesh, cplx k, EntryParams prm)
    : m_(mesh), k_(k), prm_(prm), r7_(QuadRule::dunavant7()), rn_(QuadRule::subdivided(prm.near_subdivision)),
      q7_(mesh, r7_), qn_(mesh, rn_) {
    const std::size_t N = mesh.size();
    alpha_.resize(3 * N); beta_.resize(3 * N); S_.resize(N);
    for (std::size_t t = 0; t < N; ++t) {
        const auto v = mesh.vertices(t);
        const Vec3& n = mesh.normal[t];
        const real A2 = dot(cross(v[1] - v[0], v[2] - v[0]), n);
        for (int kk = 0; kk < 3; ++kk) {
            const Vec3& a = v[(kk + 1) % 3]; const Vec3& c = v[(kk + 2) % 3];
            beta_[3 * t + kk] = cross(n, c - a) / A2;                     // grad lambda_k
            alpha_[3 * t + kk] = -dot(beta_[3 * t + kk], a);              // lambda_k(a) = 0
        }
        S_[t] = psi_matrix(mesh.area[t]);
    }
    for (Adjacency a : {Adjacency::Vertex, Adjacency::Edge, Adjacency::Coincident})
        ss_[static_cast<int>(a)] = PairRule::sauter_schwab(a, prm_.ss_order);
    if (prm_.cache_near) build_near_cache();
}

std::vector<Vec3> LinearKernelEntries::index_points() const {
    std::vector<Vec3> p(3 * m_.size());
    for (std::size_t t = 0; t < m_.size(); ++t) {
        auto v = m_.vertices(t);
        for (int a = 0; a < 3; ++a) p[3 * t + a] = (m_.centroid[t] + v[a]) * 0.5;
    }
    return p;
}

std::vector<real> LinearKernelEntries::index_sizes() const {
    std::vector<real> h(3 * m_.size());
    for (std::size_t t = 0; t < m_.size(); ++t) for (int a = 0; a < 3; ++a) h[3 * t + a] = m_.hmax[t];
    return h;
}

bool LinearKernelEntries::is_near(std::size_t i, std::size_t j) const {
    real d = norm(m_.centroid[i] - m_.centroid[j]);
    return d < prm_.near_factor * std::max(m_.hmax[i], m_.hmax[j]);
}

Adjacency LinearKernelEntries::adjacency(std::size_t i, std::size_t j) const {
    if (i == j) return Adjacency::Coincident;
    int c = 0;
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) c += (m_.T[i][a] == m_.T[j][b]);
    return c == 0 ? Adjacency::None : (c == 1 ? Adjacency::Vertex : (c == 2 ? Adjacency::Edge : Adjacency::Coincident));
}

LinearBlock LinearKernelEntries::to_psi(std::size_t i, std::size_t j, const LinearBlock& L) const {
    const auto& Si = S_[i]; const auto& Sj = S_[j];
    LinearBlock T = zero_block(), R = zero_block();
    // T = Si L  (T_kb = sum_a Si_ka L_ab), R = T Sj^T (R_kl = sum_b T_kb Sj_lb)
    for (int kk = 0; kk < 3; ++kk) for (int b = 0; b < 3; ++b) for (int a = 0; a < 3; ++a)
        for (int c = 0; c < 4; ++c) T[kk * 3 + b][c] += Si[kk * 3 + a] * L[a * 3 + b][c];
    for (int kk = 0; kk < 3; ++kk) for (int l = 0; l < 3; ++l) for (int b = 0; b < 3; ++b)
        for (int c = 0; c < 4; ++c) R[kk * 3 + l][c] += T[kk * 3 + b][c] * Sj[l * 3 + b];
    return R;
}

void LinearKernelEntries::build_near_cache() {
    auto t0 = std::chrono::steady_clock::now();
    const std::size_t N = m_.size();
    const bool sym = prm_.sauter_schwab;                                  // wie KernelEntries: symmetrische Nahregeln
    std::vector<std::vector<std::size_t>> nb(N);
#ifdef CBEM_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (long i = 0; i < static_cast<long>(N); ++i)
        for (std::size_t j = sym ? i : 0; j < N; ++j) if (is_near(i, j)) nb[i].push_back(j);
    std::vector<std::vector<LinearBlock>> val(N);
#ifdef CBEM_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
    for (long i = 0; i < static_cast<long>(N); ++i) {
        val[i].resize(nb[i].size());
        for (std::size_t a = 0; a < nb[i].size(); ++a) val[i][a] = to_psi(i, nb[i][a], lambda_exact(i, nb[i][a]));
    }
    cache_.assign(N, {});
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t a = 0; a < nb[i].size(); ++a) {
            const std::size_t j = nb[i][a]; const LinearBlock& K = val[i][a];
            cache_[i].emplace_back(j, K); ++n_near_;
            if (sym && j != i) {                                          // K(j,i)_{ba} = (s, -v) von K(i,j)_{ab}
                LinearBlock Kt;
                for (int p = 0; p < 3; ++p) for (int q = 0; q < 3; ++q) {
                    const KernelComp& k = K[p * 3 + q];
                    Kt[q * 3 + p] = KernelComp{k[0], -k[1], -k[2], -k[3]};
                }
                cache_[j].emplace_back(i, Kt); ++n_near_;
            }
        }
    for (auto& row : cache_) std::sort(row.begin(), row.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    cached_ = true;
    t_near_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

LinearBlock LinearKernelEntries::block(std::size_t i, std::size_t j) const {
    if (cached_) {
        const auto& row = cache_[i];
        auto it = std::lower_bound(row.begin(), row.end(), j, [](const auto& e, std::size_t v) { return e.first < v; });
        if (it != row.end() && it->first == j) return it->second;
        return block_far(i, j);
    }
    return to_psi(i, j, lambda_exact(i, j));
}

LinearBlock LinearKernelEntries::block_far(std::size_t i, std::size_t j) const { return to_psi(i, j, lambda_far(i, j)); }

LinearBlock LinearKernelEntries::lambda_exact(std::size_t i, std::size_t j) const {
    if (!is_near(i, j)) return lambda_far(i, j);
    // Selbstterm: Sauter-Schwab (auch ohne prm.sauter_schwab, da Phi_0 hier nicht verschwindet); gestreckte Elemente
    // halbanalytisch wie KernelEntries (Sauter-Schwab stagniert dort bei etwa 3e-5, halbanalytisch konvergiert mit sa_order)
    if (i == j) return aspect(i) <= prm_.ss_max_aspect ? lambda_sauter_schwab(i, i, Adjacency::Coincident)
                                                       : lambda_semi_analytic(i, i, Adjacency::Coincident);
    if (prm_.sauter_schwab) {
        Adjacency a = adjacency(i, j);
        if (a != Adjacency::None) {
            if (aspect(i) <= prm_.ss_max_aspect && aspect(j) <= prm_.ss_max_aspect) return lambda_sauter_schwab(i, j, a);
            return lambda_semi_analytic(i, j, a);
        }
    }
    return lambda_near(i, j);
}

LinearBlock LinearKernelEntries::lambda_far(std::size_t i, std::size_t j) const {
    LinearBlock K = zero_block();
    const Vec3* xi = q7_.points(i); const real* wi = q7_.weights(i);
    const Vec3* yj = q7_.points(j); const real* wj = q7_.weights(j);
    for (int a = 0; a < q7_.q; ++a)
        for (int b = 0; b < q7_.q; ++b) {
            Vec3 z = xi[a] - yj[b]; KernelValue kv = dirac_kernel_full(z, k_);
            acc(K, r7_.bary[a].data(), r7_.bary[b].data(), wi[a] * wj[b], kv.s, kv.vcoef, z);
        }
    return K;
}

// Beitrag eines aeusseren Punktes x (Gewicht w, Formfunktionswerte lo des aeusseren Elements) mit exaktem Innenintegral
// ueber das Element inner. swap: aeusseres Element ist j (Ansatz), inneres i (Test).
void LinearKernelEntries::outer_point(const Vec3& x, real w, const real lo[3], std::size_t inner, bool swap, LinearBlock& K) const {
    auto tri = m_.vertices(inner); const Vec3& n = m_.normal[inner]; const cplx ik = cplx(0, 1) * k_;
    std::array<Vec3, 3> Ig; std::array<real, 3> Ii;
    triangle_integrals_linear(x, tri, n, Ig, Ii);
    const real sgn = swap ? -1.0 : 1.0;
    for (int p = 0; p < 3; ++p)                                          // p: Formfunktion des aeusseren Elements
        for (int q = 0; q < 3; ++q) {                                    // q: Formfunktion des inneren Elements
            KernelComp& k = swap ? K[q * 3 + p] : K[p * 3 + q];
            const real wp = w * lo[p];
            k[1] += sgn * wp * Ig[q].x / (4 * pi); k[2] += sgn * wp * Ig[q].y / (4 * pi); k[3] += sgn * wp * Ig[q].z / (4 * pi);
            k[0] += -ik * wp * Ii[q] / (4 * pi);
        }
    const Vec3* y = q7_.points(inner); const real* wy = q7_.weights(inner);
    for (int b = 0; b < q7_.q; ++b) {
        Vec3 z = swap ? (y[b] - x) : (x - y[b]);
        KernelValue kv = dirac_kernel_remainder(z, k_);
        if (swap) acc(K, r7_.bary[b].data(), lo, w * wy[b], kv.s, kv.vcoef, z);
        else acc(K, lo, r7_.bary[b].data(), w * wy[b], kv.s, kv.vcoef, z);
    }
}

void LinearKernelEntries::near_adaptive(const std::array<Vec3, 3>& o, std::size_t outer_t, std::size_t inner, int depth, bool swap,
                                        LinearBlock& K) const {
    const Vec3 c = (o[0] + o[1] + o[2]) / 3.0;
    real rho = 0; for (auto& v : o) rho = std::max(rho, norm(v - c));
    auto tri = m_.vertices(inner);
    real d = point_triangle_distance(c, tri) - rho;
    if (prm_.adapt_to_boundary) {
        const Vec3& n = m_.normal[inner];
        const real w0 = dot(o[0] - tri[0], n), w1 = dot(o[1] - tri[0], n), w2 = dot(o[2] - tri[0], n);
        const real tolw = 1e-12 * m_.hmax[inner];
        if ((w0 > tolw && w1 > tolw && w2 > tolw) || (w0 < -tolw && w1 < -tolw && w2 < -tolw))
            d = std::max(d, point_boundary_distance(c, tri) - rho);
    }
    if (depth < prm_.adapt_depth && !(rho < prm_.adapt_ratio * d)) {
        int e = 0; real L = -1;
        for (int q = 0; q < 3; ++q) { real l = norm(o[(q + 1) % 3] - o[q]); if (l > L) { L = l; e = q; } }
        const Vec3& a = o[e]; const Vec3& b = o[(e + 1) % 3]; const Vec3& cc = o[(e + 2) % 3];
        Vec3 mid = (a + b) * 0.5;
        near_adaptive({a, mid, cc}, outer_t, inner, depth + 1, swap, K);
        near_adaptive({mid, b, cc}, outer_t, inner, depth + 1, swap, K);
        return;
    }
    const real area = 0.5 * norm(cross(o[1] - o[0], o[2] - o[0]));
    for (std::size_t a = 0; a < r7_.w.size(); ++a) {
        Vec3 x = o[0] * r7_.bary[a][0] + o[1] * r7_.bary[a][1] + o[2] * r7_.bary[a][2];
        const real lo[3] = {lambda(outer_t, 0, x), lambda(outer_t, 1, x), lambda(outer_t, 2, x)};
        outer_point(x, r7_.w[a] * area, lo, inner, swap, K);
    }
}

LinearBlock LinearKernelEntries::lambda_near(std::size_t i, std::size_t j) const {
    LinearBlock K = zero_block();
    const bool swap = (i != j) && (m_.hmax[i] > 1.5 * m_.hmax[j]);
    const std::size_t outer = swap ? j : i, inner = swap ? i : j;
    if (prm_.adaptive_outer && i != j) {
        near_adaptive(m_.vertices(outer), outer, inner, 0, swap, K);
        return K;
    }
    // feste aeussere Regel (sub^2 * 7 Punkte), inneres Integral analytisch
    const Vec3* xo = qn_.points(outer); const real* wo = qn_.weights(outer);
    for (int a = 0; a < qn_.q; ++a) outer_point(xo[a], wo[a], rn_.bary[a].data(), inner, swap, K);
    return K;
}

LinearBlock LinearKernelEntries::lambda_semi_analytic(std::size_t i, std::size_t j, Adjacency adj) const {
    // wie KernelEntries::semi_analytic; der Selbstterm wird ueber das Innenintegral (Hauptwert in der Ebene) erfasst
    std::vector<real> g, gw; gauss_legendre01(prm_.sa_order, g, gw);
    LinearBlock K = zero_block();
    auto pt = [&](const Vec3& x, real w) {
        const real lo[3] = {lambda(i, 0, x), lambda(i, 1, x), lambda(i, 2, x)};
        outer_point(x, w, lo, j, false, K);
    };
    auto edge_patch = [&](const Vec3& p0, const Vec3& p1, const Vec3& c) {
        const real A2 = norm(cross(p1 - p0, c - p0));
        for (std::size_t a = 0; a < g.size(); ++a)
            for (std::size_t b = 0; b < g.size(); ++b) {
                real t = g[a], tau = g[b], sv = tau * tau * tau, ds = 3 * tau * tau;
                Vec3 x = (p0 + (p1 - p0) * t) * (1 - sv) + c * sv;
                pt(x, gw[a] * gw[b] * ds * (1 - sv) * A2);
            }
    };
    const auto& ti = m_.T[i]; const auto& tj = m_.T[j];
    if (adj == Adjacency::Coincident) {
        auto v = m_.vertices(i); const Vec3 G = m_.centroid[i];
        for (int e = 0; e < 3; ++e) edge_patch(v[e], v[(e + 1) % 3], G);
    } else if (adj == Adjacency::Edge) {
        int sh[2], ns = 0, other = -1;
        for (int a = 0; a < 3; ++a) { bool s = false; for (int b = 0; b < 3; ++b) s |= (ti[a] == tj[b]); if (s) sh[ns++] = ti[a]; else other = ti[a]; }
        edge_patch(m_.P[sh[0]], m_.P[sh[1]], m_.P[other]);
    } else {
        int V = -1, o[2], no = 0;
        for (int a = 0; a < 3; ++a) { bool s = false; for (int b = 0; b < 3; ++b) s |= (ti[a] == tj[b]); if (s) V = ti[a]; else o[no++] = ti[a]; }
        const Vec3 P0 = m_.P[V], B = m_.P[o[0]], C = m_.P[o[1]]; const real A2 = norm(cross(B - P0, C - P0));
        for (std::size_t a = 0; a < g.size(); ++a)
            for (std::size_t b = 0; b < g.size(); ++b) {
                real tau = g[a], u = tau * tau * tau, du = 3 * tau * tau, v = g[b];
                pt(P0 + ((B - P0) + (C - B) * v) * u, gw[a] * gw[b] * du * u * A2);
            }
    }
    return K;
}

LinearBlock LinearKernelEntries::lambda_sauter_schwab(std::size_t i, std::size_t j, Adjacency adj) const {
    std::array<int, 3> ti = m_.T[i], tj = m_.T[j];
    if (adj == Adjacency::Edge || adj == Adjacency::Vertex) {             // gemeinsame Ecken zuerst (wie KernelEntries)
        std::array<int, 3> oi{}, oj{}; int ns = 0;
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) if (ti[a] == tj[b]) { oi[ns] = ti[a]; oj[ns] = tj[b]; ++ns; }
        int ki = ns, kj = ns;
        for (int a = 0; a < 3; ++a) { bool sh = false; for (int q = 0; q < ns; ++q) sh |= (ti[a] == oi[q]); if (!sh) oi[ki++] = ti[a]; }
        for (int b = 0; b < 3; ++b) { bool sh = false; for (int q = 0; q < ns; ++q) sh |= (tj[b] == oj[q]); if (!sh) oj[kj++] = tj[b]; }
        ti = oi; tj = oj;
        if (adj == Adjacency::Edge) {
            auto ang = [&](int a, int b, int c) { Vec3 u = m_.P[b] - m_.P[a], v = m_.P[c] - m_.P[a]; return std::acos(dot(u, v) / (norm(u) * norm(v))); };
            real keep = ang(tj[0], tj[1], tj[2]) - ang(ti[0], ti[1], ti[2]);
            real flip = ang(tj[1], tj[0], tj[2]) - ang(ti[1], ti[0], ti[2]);
            if (flip > keep) { std::swap(ti[0], ti[1]); std::swap(tj[0], tj[1]); }
        }
    }
    const Vec3 A = m_.P[ti[0]], Bi = m_.P[ti[1]] - A, Ci = m_.P[ti[2]] - A;
    const Vec3 D = m_.P[tj[0]], Bj = m_.P[tj[1]] - D, Cj = m_.P[tj[2]] - D;
    const PairRule& R = ss_[static_cast<int>(adj)];
    const real jac = 4.0 * m_.area[i] * m_.area[j];
    LinearBlock K = zero_block();
    const bool coinc = (adj == Adjacency::Coincident);
    for (std::size_t q = 0; q < R.w.size(); ++q) {
        Vec3 x = A + Bi * R.x[q][0] + Ci * R.x[q][1], y = D + Bj * R.y[q][0] + Cj * R.y[q][1];
        Vec3 z = x - y; real r = norm(z);
        if (r < 1e-14) continue;
        KernelValue kv = dirac_kernel_full(z, k_);
        const real w = R.w[q] * jac;
        const real lx[3] = {lambda(i, 0, x), lambda(i, 1, x), lambda(i, 2, x)};
        const real ly[3] = {lambda(j, 0, y), lambda(j, 1, y), lambda(j, 2, y)};
        if (!coinc) { acc(K, lx, ly, w, kv.s, kv.vcoef, z); continue; }
        // Selbstterm: Rest normal, Phi_0 antisymmetrisiert: [lambda_a(x) lambda_b(y) - lambda_a(y) lambda_b(x)] / 2
        const cplx p0 = 1.0 / (4 * pi * r * r * r);
        acc(K, lx, ly, w, kv.s, kv.vcoef - p0, z);
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                const real f = 0.5 * w * (lx[a] * ly[b] - ly[a] * lx[b]);   // i == j: dieselben Funktionen
                KernelComp& k = K[a * 3 + b];
                k[1] += f * p0 * z.x; k[2] += f * p0 * z.y; k[3] += f * p0 * z.z;
            }
    }
    return K;
}

}  // namespace cbem
