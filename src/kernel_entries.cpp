#include "cbem/assembly/kernel_entries.hpp"
#include "cbem/kernel/dirac_kernel.hpp"
#include "cbem/kernel/triangle_integrals.hpp"
#include <algorithm>
#include <chrono>

namespace cbem {

KernelEntries::KernelEntries(const TriangleMesh& mesh, cplx k, EntryParams prm)
    : m_(mesh), k_(k), prm_(prm), q7_(mesh, QuadRule::dunavant7()), q2_(mesh, QuadRule::subdivided(2)),
      qn_(mesh, QuadRule::subdivided(prm.near_subdivision)) {
    rho_.resize(mesh.size());
    for (std::size_t t = 0; t < mesh.size(); ++t) { auto v = mesh.vertices(t); real r = 0; for (auto& p : v) r = std::max(r, norm(p - mesh.centroid[t])); rho_[t] = r; }
    if (prm_.sauter_schwab)
        for (Adjacency a : {Adjacency::Vertex, Adjacency::Edge, Adjacency::Coincident})
            ss_[static_cast<int>(a)] = PairRule::sauter_schwab(a, prm_.ss_order);
    if (prm_.cache_near) build_near_cache();
}

void KernelEntries::build_near_cache() {
    auto t0 = std::chrono::steady_clock::now();
    const std::size_t N = m_.size();
    // Nahpaare suchen: Schwerpunktabstand < near_factor * max(h_i, h_j). Mit Sauter-Schwab sind die
    // Nahfeldregeln symmetrisch (K_ji = (s, -v)), dann nur i <= j; sonst beide Haelften einzeln.
    const bool sym = prm_.sauter_schwab;
    std::vector<std::vector<std::size_t>> nb(N);
#ifdef CBEM_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (long i = 0; i < static_cast<long>(N); ++i)
        for (std::size_t j = sym ? i : 0; j < N; ++j) if (is_near(i, j)) nb[i].push_back(j);
    // Werte berechnen (obere Haelfte), untere Haelfte per Symmetrie
    std::vector<std::vector<KernelComp>> val(N);
#ifdef CBEM_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
    for (long i = 0; i < static_cast<long>(N); ++i) {
        val[i].resize(nb[i].size());
        for (std::size_t a = 0; a < nb[i].size(); ++a) val[i][a] = exact_uncached(i, nb[i][a]);
    }
    cache_.assign(N, {});
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t a = 0; a < nb[i].size(); ++a) {
            std::size_t j = nb[i][a]; const KernelComp& K = val[i][a];
            cache_[i].emplace_back(j, K); ++n_near_;
            if (sym && j != i) { cache_[j].emplace_back(i, KernelComp{K[0], -K[1], -K[2], -K[3]}); ++n_near_; }
        }
    for (auto& row : cache_) std::sort(row.begin(), row.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    cached_ = true;
    t_near_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

Adjacency KernelEntries::adjacency(std::size_t i, std::size_t j) const {
    if (i == j) return Adjacency::Coincident;
    int c = 0;
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) c += (m_.T[i][a] == m_.T[j][b]);
    return c == 0 ? Adjacency::None : (c == 1 ? Adjacency::Vertex : (c == 2 ? Adjacency::Edge : Adjacency::Coincident));
}

KernelComp KernelEntries::exact(std::size_t i, std::size_t j) const {
    if (cached_) {
        const auto& row = cache_[i];
        auto it = std::lower_bound(row.begin(), row.end(), j, [](const auto& e, std::size_t v) { return e.first < v; });
        if (it != row.end() && it->first == j) return it->second;
        return far(i, j);
    }
    return exact_uncached(i, j);
}

KernelComp KernelEntries::exact_uncached(std::size_t i, std::size_t j) const {
    if (!is_near(i, j)) return far(i, j);
    if (prm_.sauter_schwab) {
        Adjacency a = adjacency(i, j);
        if (a != Adjacency::None) {
            if (aspect(i) <= prm_.ss_max_aspect && aspect(j) <= prm_.ss_max_aspect) return sauter_schwab(i, j, a);
            return semi_analytic(i, j, a);
        }
    }
    return near(i, j);
}

namespace {
real point_triangle_distance(const Vec3& x, const std::array<Vec3, 3>& p) {
    Vec3 n = cross(p[1] - p[0], p[2] - p[0]); real a2 = norm(n); n = n / a2;
    real w = dot(x - p[0], n); Vec3 q = x - n * w;
    // baryzentrische Lage von q
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

// Beitrag eines aeusseren Punktes x (Gewicht w) mit exaktem Innenintegral ueber das Dreieck inner
void KernelEntries::outer_point(const Vec3& x, real w, std::size_t inner, bool skip_phi0, KernelComp& K) const {
    auto tri = m_.vertices(inner); const Vec3& n = m_.normal[inner]; const cplx ik = cplx(0, 1) * k_;
    Vec3 Ig; real Ii; triangle_integrals(x, tri, n, Ig, Ii);
    if (!skip_phi0) { K[1] += w * Ig.x / (4 * pi); K[2] += w * Ig.y / (4 * pi); K[3] += w * Ig.z / (4 * pi); }
    K[0] += -ik * w * Ii / (4 * pi);
    const Vec3* y = q7_.points(inner); const real* wy = q7_.weights(inner);
    for (int b = 0; b < q7_.q; ++b) {
        Vec3 z = x - y[b]; KernelValue kv = dirac_kernel_remainder(z, k_); real ww = w * wy[b];
        K[0] += ww * kv.s; K[1] += ww * kv.vcoef * z.x; K[2] += ww * kv.vcoef * z.y; K[3] += ww * kv.vcoef * z.z;
    }
}

// Halbanalytische Regel fuer benachbarte Paare: inneres Integral ueber tau_j exakt (Wilton), aeussere
// Regel auf tau_i in Koordinaten, die zur gemeinsamen Menge hin kubisch gradiert sind:
//   Kante:  x = (1 - s)[p0 + t (p1 - p0)] + s c,  s = tau^3   (Singularitaet bei s = 0)
//   Ecke:   x = V + u [(B - V) + v (C - B)],      u = tau^3
//   gleich: drei Teildreiecke (Schwerpunkt, Kante), je wie "Kante" zur Randkante gradiert.
// Das Innenintegral ist stetig; seine logarithmischen Singularitaeten liegen auf der gemeinsamen Menge.
KernelComp KernelEntries::semi_analytic(std::size_t i, std::size_t j, Adjacency adj) const {
    std::vector<real> g, gw; gauss_legendre01(prm_.sa_order, g, gw);
    KernelComp K{0, 0, 0, 0};
    auto edge_patch = [&](const Vec3& p0, const Vec3& p1, const Vec3& c, bool skip) {
        const real A2 = norm(cross(p1 - p0, c - p0));
        for (std::size_t a = 0; a < g.size(); ++a)
            for (std::size_t b = 0; b < g.size(); ++b) {
                real t = g[a], tau = g[b], sv = tau * tau * tau, ds = 3 * tau * tau;
                Vec3 x = (p0 + (p1 - p0) * t) * (1 - sv) + c * sv;
                outer_point(x, gw[a] * gw[b] * ds * (1 - sv) * A2, j, skip, K);
            }
    };
    const auto& ti = m_.T[i]; const auto& tj = m_.T[j];
    if (adj == Adjacency::Coincident) {
        auto v = m_.vertices(i); const Vec3 G = m_.centroid[i];
        for (int e = 0; e < 3; ++e) edge_patch(v[e], v[(e + 1) % 3], G, true);   // Phi_0-Selbstterm = 0
    } else if (adj == Adjacency::Edge) {
        int sh[2], ns = 0, other = -1;
        for (int a = 0; a < 3; ++a) { bool s = false; for (int b = 0; b < 3; ++b) s |= (ti[a] == tj[b]); if (s) sh[ns++] = ti[a]; else other = ti[a]; }
        edge_patch(m_.P[sh[0]], m_.P[sh[1]], m_.P[other], false);
    } else {
        int V = -1, o[2], no = 0;
        for (int a = 0; a < 3; ++a) { bool s = false; for (int b = 0; b < 3; ++b) s |= (ti[a] == tj[b]); if (s) V = ti[a]; else o[no++] = ti[a]; }
        const Vec3 P0 = m_.P[V], B = m_.P[o[0]], C = m_.P[o[1]]; const real A2 = norm(cross(B - P0, C - P0));
        for (std::size_t a = 0; a < g.size(); ++a)
            for (std::size_t b = 0; b < g.size(); ++b) {
                real tau = g[a], u = tau * tau * tau, du = 3 * tau * tau, v = g[b];
                Vec3 x = P0 + ((B - P0) + (C - B) * v) * u;
                outer_point(x, gw[a] * gw[b] * du * u * A2, j, false, K);
            }
    }
    return K;
}

void KernelEntries::near_adaptive(const std::array<Vec3, 3>& o, std::size_t inner, int depth, bool swap, bool, KernelComp& K) const {
    const Vec3 c = (o[0] + o[1] + o[2]) / 3.0;
    real rho = 0; for (auto& v : o) rho = std::max(rho, norm(v - c));
    auto tri = m_.vertices(inner);
    real d = point_triangle_distance(c, tri) - rho;
    if (prm_.adapt_to_boundary) {
        // Teilstueck ganz auf einer Seite der Ebene: Innenintegrale (Phi_0, 1/r) dort reell-analytisch mit
        // Singularitaeten nur auf dem Rand des inneren Dreiecks -> Abstand zum Rand massgeblich
        const Vec3& n = m_.normal[inner];
        const real w0 = dot(o[0] - tri[0], n), w1 = dot(o[1] - tri[0], n), w2 = dot(o[2] - tri[0], n);
        const real tolw = 1e-12 * m_.hmax[inner];
        if ((w0 > tolw && w1 > tolw && w2 > tolw) || (w0 < -tolw && w1 < -tolw && w2 < -tolw))
            d = std::max(d, point_boundary_distance(c, tri) - rho);
    }
    if (depth < prm_.adapt_depth && !(rho < prm_.adapt_ratio * d)) {
        // Halbierung der laengsten Kante
        int e = 0; real L = -1;
        for (int q = 0; q < 3; ++q) { real l = norm(o[(q + 1) % 3] - o[q]); if (l > L) { L = l; e = q; } }
        const Vec3& a = o[e]; const Vec3& b = o[(e + 1) % 3]; const Vec3& cc = o[(e + 2) % 3];
        Vec3 mid = (a + b) * 0.5;
        near_adaptive({a, mid, cc}, inner, depth + 1, swap, false, K);
        near_adaptive({mid, b, cc}, inner, depth + 1, swap, false, K);
        return;
    }
    // 7-Punkt-Regel auf dem Teilstueck; inneres Integral analytisch (Phi_0, 1/r), Rest mit 7 Punkten
    static const QuadRule R7 = QuadRule::dunavant7();
    const real area = 0.5 * norm(cross(o[1] - o[0], o[2] - o[0]));
    const Vec3& n = m_.normal[inner]; const cplx ik = cplx(0, 1) * k_;
    const real sgn = swap ? -1.0 : 1.0;
    const Vec3* yin = q7_.points(inner); const real* win = q7_.weights(inner); const int nin = q7_.q;
    for (std::size_t a = 0; a < R7.w.size(); ++a) {
        Vec3 x = o[0] * R7.bary[a][0] + o[1] * R7.bary[a][1] + o[2] * R7.bary[a][2]; real wx = R7.w[a] * area;
        Vec3 Ig; real Ii; triangle_integrals(x, tri, n, Ig, Ii);
        K[1] += sgn * wx * Ig.x / (4 * pi); K[2] += sgn * wx * Ig.y / (4 * pi); K[3] += sgn * wx * Ig.z / (4 * pi);
        K[0] += -ik * wx * Ii / (4 * pi);
        for (int b = 0; b < nin; ++b) {
            Vec3 z = swap ? (yin[b] - x) : (x - yin[b]);
            KernelValue kv = dirac_kernel_remainder(z, k_); real ww = wx * win[b];
            K[0] += ww * kv.s; K[1] += ww * kv.vcoef * z.x; K[2] += ww * kv.vcoef * z.y; K[3] += ww * kv.vcoef * z.z;
        }
    }
}

KernelComp KernelEntries::sauter_schwab(std::size_t i, std::size_t j, Adjacency adj) const {
    // Eckenreihenfolge: gemeinsame Ecken zuerst, in beiden Dreiecken in gleicher Reihenfolge
    std::array<int, 3> ti = m_.T[i], tj = m_.T[j];
    if (adj == Adjacency::Edge || adj == Adjacency::Vertex) {
        std::array<int, 3> oi{}, oj{}; int ns = 0;
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) if (ti[a] == tj[b]) { oi[ns] = ti[a]; oj[ns] = tj[b]; ++ns; }
        int ki = ns, kj = ns;
        for (int a = 0; a < 3; ++a) { bool sh = false; for (int q = 0; q < ns; ++q) sh |= (ti[a] == oi[q]); if (!sh) oi[ki++] = ti[a]; }
        for (int b = 0; b < 3; ++b) { bool sh = false; for (int q = 0; q < ns; ++q) sh |= (tj[b] == oj[q]); if (!sh) oj[kj++] = tj[b]; }
        ti = oi; tj = oj;
        if (adj == Adjacency::Edge) {
            // Orientierung der gemeinsamen Kante: die Regel konvergiert schnell, wenn das Ansatzdreieck an der
            // ersten gemeinsamen Ecke einen grossen, das Testdreieck einen kleinen Winkel hat (numerisch ermittelt)
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
    KernelComp K{0, 0, 0, 0};
    const bool coinc = (adj == Adjacency::Coincident);
    for (std::size_t q = 0; q < R.w.size(); ++q) {
        Vec3 x = A + Bi * R.x[q][0] + Ci * R.x[q][1], y = D + Bj * R.y[q][0] + Cj * R.y[q][1];
        Vec3 z = x - y; real r = norm(z);
        if (r < 1e-14) continue;
        KernelValue kv = dirac_kernel_full(z, k_);
        cplx vc = kv.vcoef;
        if (coinc) vc -= 1.0 / (4 * pi * r * r * r);          // Phi_0 im Selbstterm exakt 0 (Antisymmetrie)
        real w = R.w[q] * jac;
        K[0] += w * kv.s; K[1] += w * vc * z.x; K[2] += w * vc * z.y; K[3] += w * vc * z.z;
    }
    return K;
}

bool KernelEntries::is_near(std::size_t i, std::size_t j) const {
    real d = norm(m_.centroid[i] - m_.centroid[j]);
    return d < prm_.near_factor * std::max(m_.hmax[i], m_.hmax[j]);
}

KernelComp KernelEntries::far(std::size_t i, std::size_t j) const {
    KernelComp K{0, 0, 0, 0};
    const Vec3* xi = q7_.points(i); const real* wi = q7_.weights(i);
    const Vec3* yj = q7_.points(j); const real* wj = q7_.weights(j);
    for (int a = 0; a < q7_.q; ++a)
        for (int b = 0; b < q7_.q; ++b) {
            Vec3 z = xi[a] - yj[b]; KernelValue kv = dirac_kernel_full(z, k_); real ww = wi[a] * wj[b];
            K[0] += ww * kv.s; K[1] += ww * kv.vcoef * z.x; K[2] += ww * kv.vcoef * z.y; K[3] += ww * kv.vcoef * z.z;
        }
    return K;
}

KernelComp KernelEntries::near(std::size_t i, std::size_t j) const {
    KernelComp K{0, 0, 0, 0};
    const bool swap = (i != j) && (m_.hmax[i] > 1.5 * m_.hmax[j]);
    // aeussere Punkte auf dem kleineren Dreieck, inneres Integral analytisch ueber das andere
    std::size_t outer = swap ? j : i, inner = swap ? i : j;
    if (prm_.adaptive_outer && i != j) {
        KernelComp K{0, 0, 0, 0};
        near_adaptive(m_.vertices(outer), inner, 0, swap, false, K);
        return K;
    }
    const MeshQuadrature* qo = &qn_;
    const Vec3* xo = qo->points(outer); const real* wo = qo->weights(outer);
    auto tri = m_.vertices(inner); const Vec3& n = m_.normal[inner];
    const cplx ik = cplx(0, 1) * k_;
    Vec3 g{0, 0, 0}; real inv = 0;
    for (int a = 0; a < qo->q; ++a) {
        Vec3 Ig; real Ii; triangle_integrals(xo[a], tri, n, Ig, Ii);
        g += Ig * wo[a]; inv += Ii * wo[a];
    }
    // Phi_0-Anteil: ungerade in (x - y); Selbstterm exakt 0
    real sgn = swap ? -1.0 : 1.0;
    if (i != j) { K[1] += sgn * g.x / (4 * pi); K[2] += sgn * g.y / (4 * pi); K[3] += sgn * g.z / (4 * pi); }
    K[0] += -ik * inv / (4 * pi);
    // Rest mit Gauss: aeussere Punkte x (7 Punkte auf i bzw. Nahregel), z = x_i - y_j
    const Vec3* yi = swap ? q7_.points(i) : xo;  const real* wyi = swap ? q7_.weights(i) : wo;  int ni = swap ? q7_.q : qo->q;
    const Vec3* yj = swap ? xo : q7_.points(j); const real* wyj = swap ? wo : q7_.weights(j); int nj = swap ? qo->q : q7_.q;
    for (int a = 0; a < ni; ++a)
        for (int b = 0; b < nj; ++b) {
            Vec3 z = yi[a] - yj[b]; KernelValue kv = dirac_kernel_remainder(z, k_); real ww = wyi[a] * wyj[b];
            K[0] += ww * kv.s; K[1] += ww * kv.vcoef * z.x; K[2] += ww * kv.vcoef * z.y; K[3] += ww * kv.vcoef * z.z;
        }
    return K;
}

}  // namespace cbem
