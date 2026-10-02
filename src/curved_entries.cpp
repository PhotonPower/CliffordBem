#include "cbem/assembly/curved_entries.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>

#include "cbem/kernel/dirac_kernel.hpp"
#include "cbem/kernel/triangle_integrals.hpp"

namespace cbem {

namespace {

// Phi n = s n + v (z . n) + v (z ^ n) mit dem Kernkoeffizienten v (vcoef) und dem Skalar s
inline CurvedComp comps7(const Vec3& z, const Vec3& n, cplx s, cplx v) {
    return CurvedComp{v * dot(z, n), s * n.x, s * n.y, s * n.z,
                      v * (z.x * n.y - z.y * n.x), v * (z.x * n.z - z.z * n.x), v * (z.y * n.z - z.z * n.y)};
}

inline void acc(CurvedBlock& K, const std::array<real, 3>& la, const std::array<real, 3>& lb, real w, const CurvedComp& c) {
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            const real ww = w * la[a] * lb[b];
            CurvedComp& k = K[a * 3 + b];
            for (int q = 0; q < kCurvedComps; ++q) k[q] += ww * c[q];
        }
}

CurvedBlock zero_block() {
    CurvedBlock K;
    for (auto& k : K) k.fill(cplx(0));
    return K;
}

real point_triangle_distance(const Vec3& x, const std::array<Vec3, 3>& p) {
    Vec3 n = cross(p[1] - p[0], p[2] - p[0]); n = n / norm(n);
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

using Tri = std::array<std::array<real, 3>, 3>;   // Teildreieck in baryzentrischen Parametern
const Tri kRef = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};

}  // namespace

CurvedKernelEntries::CurvedKernelEntries(const QuadraticMesh& mesh, cplx k, EntryParams prm, CurvedNearParams np)
    : m_(mesh), k_(k), prm_(prm), np_(np), r7_(QuadRule::dunavant7()), ro_(np.outer_rule >= 2 ? QuadRule::conical(np.outer_rule) : r7_),
      rc_(np.correction_rule >= 2 ? QuadRule::conical(np.correction_rule) : r7_), q7_(mesh, r7_), S_(curved_psi_matrices(mesh)) {
    for (Adjacency a : {Adjacency::Vertex, Adjacency::Edge, Adjacency::Coincident}) {
        const auto& o = np_.ss_orders[static_cast<int>(a) - 1];
        ss_[static_cast<int>(a)] = o[0] > 0 ? PairRule::sauter_schwab(a, o) : PairRule::sauter_schwab(a, prm_.ss_order);
    }
    poly_.resize(m_.size());
    for (std::size_t t = 0; t < m_.size(); ++t) {
        // Formfunktionen in u = lambda_1, v = lambda_2 ausmultipliziert (Ecken V, Kantenmitten M01, M12, M20)
        const auto V = m_.flat.vertices(t); const auto& M = m_.mid[t];
        Poly& g = poly_[t];
        g.A = V[0];
        g.B = V[0] * -3.0 - V[1] + M[0] * 4.0;
        g.C = V[0] * -3.0 - V[2] + M[2] * 4.0;
        g.D = (V[0] + V[1]) * 2.0 - M[0] * 4.0;
        g.E = (V[0] - M[0] + M[1] - M[2]) * 4.0;
        g.F = (V[0] + V[2]) * 2.0 - M[2] * 4.0;
    }
    psiw_.resize(q7_.w.size());
    for (std::size_t t = 0; t < m_.size(); ++t)
        for (int p = 0; p < q7_.q; ++p)
            for (int a = 0; a < 3; ++a) {
                real s = 0;
                for (int k = 0; k < 3; ++k) s += S_[t][a * 3 + k] * q7_.lam[p][k];
                psiw_[t * q7_.q + p][a] = q7_.weights(t)[p] * s;
            }
    if (prm_.cache_near) build_near_cache();
}

std::vector<Vec3> CurvedKernelEntries::index_points() const {
    std::vector<Vec3> p(3 * m_.size());
    for (std::size_t t = 0; t < m_.size(); ++t) {
        const auto v = m_.flat.vertices(t);
        for (int a = 0; a < 3; ++a) p[3 * t + a] = (m_.flat.centroid[t] + v[a]) * 0.5;
    }
    return p;
}

bool CurvedKernelEntries::is_near(std::size_t i, std::size_t j) const {
    const real d = norm(m_.flat.centroid[i] - m_.flat.centroid[j]);
    return d < prm_.near_factor * std::max(m_.flat.hmax[i], m_.flat.hmax[j]);
}

Adjacency CurvedKernelEntries::adjacency(std::size_t i, std::size_t j) const {
    if (i == j) return Adjacency::Coincident;
    int c = 0;
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) c += (m_.flat.T[i][a] == m_.flat.T[j][b]);
    return c == 0 ? Adjacency::None : (c == 1 ? Adjacency::Vertex : (c == 2 ? Adjacency::Edge : Adjacency::Coincident));
}

real CurvedKernelEntries::distance_to_element(const Vec3& x, std::size_t j) const {
    // untere Schranke: Abstand zum Sehnendreieck minus groesste Abweichung der Flaeche davon (4/3 der Kantenmitten-Abweichung)
    return point_triangle_distance(x, m_.flat.vertices(j)) - 4.0 / 3.0 * m_.bulge[j];
}

CurvedBlock CurvedKernelEntries::to_psi(std::size_t i, std::size_t j, const CurvedBlock& L) const {
    const auto& Si = S_[i]; const auto& Sj = S_[j];
    CurvedBlock T = zero_block(), R = zero_block();
    for (int k = 0; k < 3; ++k) for (int b = 0; b < 3; ++b) for (int a = 0; a < 3; ++a)
        for (int c = 0; c < kCurvedComps; ++c) T[k * 3 + b][c] += Si[k * 3 + a] * L[a * 3 + b][c];
    for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) for (int b = 0; b < 3; ++b)
        for (int c = 0; c < kCurvedComps; ++c) R[k * 3 + l][c] += T[k * 3 + b][c] * Sj[l * 3 + b];
    return R;
}

void CurvedKernelEntries::build_near_cache() {
    auto t0 = std::chrono::steady_clock::now();
    const std::size_t N = m_.size();
    std::vector<std::vector<std::size_t>> nb(N);
    std::vector<std::vector<CurvedBlock>> val(N);
#ifdef CBEM_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
    for (long i = 0; i < static_cast<long>(N); ++i) {
        for (std::size_t j = 0; j < N; ++j) if (is_near(i, j)) nb[i].push_back(j);   // kein Symmetrieschluss: n(y) im Kern
        val[i].resize(nb[i].size());
        for (std::size_t a = 0; a < nb[i].size(); ++a) val[i][a] = to_psi(i, nb[i][a], lambda_exact(i, nb[i][a]));
    }
    cache_.assign(N, {});
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t a = 0; a < nb[i].size(); ++a) cache_[i].emplace_back(nb[i][a], val[i][a]);
        n_near_ += nb[i].size();
    }
    cached_ = true;
    t_near_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

CurvedBlock CurvedKernelEntries::block(std::size_t i, std::size_t j) const {
    if (cached_) {
        const auto& row = cache_[i];
        auto it = std::lower_bound(row.begin(), row.end(), j, [](const auto& e, std::size_t v) { return e.first < v; });
        if (it != row.end() && it->first == j) return it->second;
        return block_far(i, j);
    }
    return to_psi(i, j, lambda_exact(i, j));
}

CurvedBlock CurvedKernelEntries::block_far(std::size_t i, std::size_t j) const {
    // wie to_psi(lambda_far): Gewichte w psi_a statt w lambda_a, zweistufig kontrahiert (erst ueber y, dann ueber x)
    CurvedBlock K = zero_block();
    const int nq = q7_.q;
    const Vec3* x = q7_.points(i); const auto* gx = &psiw_[i * nq];
    const Vec3* y = q7_.points(j); const Vec3* ny = q7_.normals(j); const auto* gy = &psiw_[j * nq];
    for (int p = 0; p < nq; ++p) {
        std::array<CurvedComp, 3> T;
        for (auto& t : T) t.fill(cplx(0));
        for (int q = 0; q < nq; ++q) {
            const Vec3 z = x[p] - y[q]; const KernelValue kv = dirac_kernel_fast(z, k_);
            const CurvedComp c = comps7(z, ny[q], kv.s, kv.vcoef);
            for (int b = 0; b < 3; ++b) { const real g = gy[q][b]; for (int cc = 0; cc < kCurvedComps; ++cc) T[b][cc] += g * c[cc]; }
        }
        for (int a = 0; a < 3; ++a) {
            const real g = gx[p][a];
            for (int b = 0; b < 3; ++b) { CurvedComp& k = K[a * 3 + b]; for (int cc = 0; cc < kCurvedComps; ++cc) k[cc] += g * T[b][cc]; }
        }
    }
    return K;
}

CurvedBlock CurvedKernelEntries::lambda_exact(std::size_t i, std::size_t j) const {
    if (!is_near(i, j)) return lambda_far(i, j);
    const Adjacency a = adjacency(i, j);
    if (a != Adjacency::None) return lambda_sauter_schwab(i, j, a);
    return lambda_near(i, j);
}

CurvedBlock CurvedKernelEntries::lambda_far(std::size_t i, std::size_t j) const {
    CurvedBlock K = zero_block();
    const Vec3* x = q7_.points(i); const real* wx = q7_.weights(i);
    const Vec3* y = q7_.points(j); const Vec3* ny = q7_.normals(j); const real* wy = q7_.weights(j);
    for (int p = 0; p < q7_.q; ++p)
        for (int q = 0; q < q7_.q; ++q) {
            const Vec3 z = x[p] - y[q]; const KernelValue kv = dirac_kernel_fast(z, k_);
            acc(K, q7_.lam[p], q7_.lam[q], wx[p] * wy[q], comps7(z, ny[q], kv.s, kv.vcoef));
        }
    return K;
}

void CurvedKernelEntries::inner(const Vec3& x, std::size_t j, const Tri& tri, real aref, int depth, std::array<CurvedComp, 3>& out) const {
    const Vec3 c0 = gX(j, tri[0]), c1 = gX(j, tri[1]), c2 = gX(j, tri[2]);
    const Vec3 c = (c0 + c1 + c2) / 3.0;
    const real rho = std::max(norm(c0 - c), std::max(norm(c1 - c), norm(c2 - c)));
    const real d = norm(x - c) - rho;
    if (depth < np_.inner_depth && !(rho < np_.inner_ratio * d)) {
        const real l01 = norm(c1 - c0), l12 = norm(c2 - c1), l20 = norm(c0 - c2);
        int e = (l01 >= l12 && l01 >= l20) ? 0 : (l12 >= l20 ? 1 : 2);
        const auto& A = tri[e]; const auto& B = tri[(e + 1) % 3]; const auto& C = tri[(e + 2) % 3];
        const std::array<real, 3> M = {(A[0] + B[0]) / 2, (A[1] + B[1]) / 2, (A[2] + B[2]) / 2};
        inner(x, j, {A, M, C}, aref / 2, depth + 1, out);
        inner(x, j, {M, B, C}, aref / 2, depth + 1, out);
        return;
    }
    for (std::size_t p = 0; p < r7_.w.size(); ++p) {
        std::array<real, 3> l;
        for (int k = 0; k < 3; ++k) l[k] = r7_.bary[p][0] * tri[0][k] + r7_.bary[p][1] * tri[1][k] + r7_.bary[p][2] * tri[2][k];
        Vec3 n; const real J = gJ(j, l, &n);
        const Vec3 z = x - gX(j, l); const KernelValue kv = dirac_kernel_fast(z, k_);
        const CurvedComp cc = comps7(z, n, kv.s, kv.vcoef);
        const real w = r7_.w[p] * aref * J;
        for (int b = 0; b < 3; ++b) for (int q = 0; q < kCurvedComps; ++q) out[b][q] += w * l[b] * cc[q];
    }
}

void CurvedKernelEntries::outer(std::size_t i, std::size_t j, const Tri& tri, real aref, int depth, CurvedBlock& K) const {
    const Vec3 c0 = gX(i, tri[0]), c1 = gX(i, tri[1]), c2 = gX(i, tri[2]);
    const Vec3 c = (c0 + c1 + c2) / 3.0;
    const real rho = std::max(norm(c0 - c), std::max(norm(c1 - c), norm(c2 - c)));
    const real d = distance_to_element(c, j) - rho;
    const real ratio = np_.subtract ? np_.subtract_outer_ratio : np_.outer_ratio;
    if (depth < np_.outer_depth && !(rho < ratio * d)) {
        const real l01 = norm(c1 - c0), l12 = norm(c2 - c1), l20 = norm(c0 - c2);
        int e = (l01 >= l12 && l01 >= l20) ? 0 : (l12 >= l20 ? 1 : 2);
        const auto& A = tri[e]; const auto& B = tri[(e + 1) % 3]; const auto& C = tri[(e + 2) % 3];
        const std::array<real, 3> M = {(A[0] + B[0]) / 2, (A[1] + B[1]) / 2, (A[2] + B[2]) / 2};
        outer(i, j, {A, M, C}, aref / 2, depth + 1, K);
        outer(i, j, {M, B, C}, aref / 2, depth + 1, K);
        return;
    }
    for (std::size_t p = 0; p < ro_.w.size(); ++p) {
        std::array<real, 3> l;
        for (int k = 0; k < 3; ++k) l[k] = ro_.bary[p][0] * tri[0][k] + ro_.bary[p][1] * tri[1][k] + ro_.bary[p][2] * tri[2][k];
        const real w = ro_.w[p] * aref * gJ(i, l);
        std::array<CurvedComp, 3> in; for (auto& v : in) v.fill(cplx(0));
        if (np_.subtract) inner_subtracted(gX(i, l), j, in);
        else inner(gX(i, l), j, kRef, 0.5, 0, in);
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) for (int q = 0; q < kCurvedComps; ++q) K[a * 3 + b][q] += w * l[a] * in[b][q];
    }
}

CurvedKernelEntries::Tangent CurvedKernelEntries::tangent_at(const Vec3& x, std::size_t j) const {
    // Fusspunkt: Projektion auf das Sehnendreieck (baryzentrisch, in das Dreieck geklemmt), dann zwei Newton-Schritte auf der
    // gekruemmten Flaeche (Genauigkeit beeinflusst nur die Wirksamkeit der Subtraktion, nicht ihre Richtigkeit)
    const auto v = m_.flat.vertices(j); const Vec3& nf = m_.flat.normal[j];
    const real A2 = dot(cross(v[1] - v[0], v[2] - v[0]), nf);
    std::array<real, 3> l;
    for (int k = 0; k < 3; ++k) l[k] = dot(cross(v[(k + 2) % 3] - v[(k + 1) % 3], x - v[(k + 1) % 3]), nf) / A2;
    auto clamp = [](std::array<real, 3>& q) {
        real s = 0; for (auto& c : q) { c = std::max(c, 0.0); s += c; }
        for (auto& c : q) c /= s;
    };
    clamp(l);
    for (int it = 0; it < 2; ++it) {
        Vec3 Xu, Xv; gFrame(j, l, Xu, Xv);
        const Vec3 r = x - gX(j, l);
        const real a11 = dot(Xu, Xu), a12 = dot(Xu, Xv), a22 = dot(Xv, Xv), b1 = dot(r, Xu), b2 = dot(r, Xv);
        const real det = a11 * a22 - a12 * a12;
        const real du = (a22 * b1 - a12 * b2) / det, dv = (a11 * b2 - a12 * b1) / det;
        l = {l[0] - du - dv, l[1] + du, l[2] + dv};
        clamp(l);
    }
    Tangent T; T.lam = l;
    gFrame(j, l, T.Xu, T.Xv);
    const Vec3 c = cross(T.Xu, T.Xv); T.J = norm(c); T.n = c / T.J;
    T.X0 = gX(j, l) - T.Xu * l[1] - T.Xv * l[2];                     // X_aff(u, v) = X0 + Xu u + Xv v
    return T;
}

void CurvedKernelEntries::inner_subtracted(const Vec3& x, std::size_t j, std::array<CurvedComp, 3>& out) const {
    const Tangent T = tangent_at(x, j);
    // analytisch: singulaerer Kern ueber dem Tangentialdreieck mit n(u*), Gewichte lambda_b (dS = J(u*) du dv)
    const std::array<Vec3, 3> tri = {T.X0, T.X0 + T.Xu, T.X0 + T.Xv};
    std::array<Vec3, 3> Ig; std::array<real, 3> Ii;
    triangle_integrals_linear(x, tri, T.n, Ig, Ii);
    const cplx ik = cplx(0, 1) * k_;
    for (int b = 0; b < 3; ++b) {
        const CurvedComp c = comps7(Ig[b] / (4 * pi), T.n, -ik * Ii[b] / (4 * pi), 1.0);
        for (int q = 0; q < kCurvedComps; ++q) out[b][q] += c[q];
    }
    correction(x, j, T, kRef, 0.5, 0, out);
}

void CurvedKernelEntries::correction(const Vec3& x, std::size_t j, const Tangent& T, const Tri& tri, real aref, int depth,
                                     std::array<CurvedComp, 3>& out) const {
    const Vec3 c0 = gX(j, tri[0]), c1 = gX(j, tri[1]), c2 = gX(j, tri[2]);
    const Vec3 c = (c0 + c1 + c2) / 3.0;
    const real rho = std::max(norm(c0 - c), std::max(norm(c1 - c), norm(c2 - c)));
    const real d = norm(x - c) - rho;
    if (depth < np_.inner_depth && !(rho < np_.correction_ratio * d)) {
        const real l01 = norm(c1 - c0), l12 = norm(c2 - c1), l20 = norm(c0 - c2);
        int e = (l01 >= l12 && l01 >= l20) ? 0 : (l12 >= l20 ? 1 : 2);
        const auto& A = tri[e]; const auto& B = tri[(e + 1) % 3]; const auto& C = tri[(e + 2) % 3];
        const std::array<real, 3> M = {(A[0] + B[0]) / 2, (A[1] + B[1]) / 2, (A[2] + B[2]) / 2};
        correction(x, j, T, {A, M, C}, aref / 2, depth + 1, out);
        correction(x, j, T, {M, B, C}, aref / 2, depth + 1, out);
        return;
    }
    for (std::size_t p = 0; p < rc_.w.size(); ++p) {
        std::array<real, 3> l;
        for (int k = 0; k < 3; ++k) l[k] = rc_.bary[p][0] * tri[0][k] + rc_.bary[p][1] * tri[1][k] + rc_.bary[p][2] * tri[2][k];
        correction_point(x, j, T, l, rc_.w[p] * aref, out);
    }
}

// Rest im Parameterpunkt l mit Gewicht w (Referenzmass du dv): lambda_b (J f(x - X) n - J(u*) f_sing(x - X_aff) n(u*))
void CurvedKernelEntries::correction_point(const Vec3& x, std::size_t j, const Tangent& T, const std::array<real, 3>& l, real w,
                                           std::array<CurvedComp, 3>& out) const {
    const cplx ik = cplx(0, 1) * k_;
    Vec3 n; const real J = gJ(j, l, &n);
    const Vec3 z = x - gX(j, l); const KernelValue kv = dirac_kernel_fast(z, k_);
    const CurvedComp full = comps7(z, n, kv.s, kv.vcoef);
    const Vec3 za = x - (T.X0 + T.Xu * l[1] + T.Xv * l[2]);              // Tangentialdreieck, gleicher Parameter
    const real ra = norm(za);
    const CurvedComp sing = comps7(za, T.n, -ik / (4 * pi * ra), 1.0 / (4 * pi * ra * ra * ra));
    for (int b = 0; b < 3; ++b) {
        const real wb = w * l[b];
        for (int q = 0; q < kCurvedComps; ++q) out[b][q] += wb * (J * full[q] - T.J * sing[q]);
    }
}

CurvedBlock CurvedKernelEntries::lambda_near(std::size_t i, std::size_t j) const {
    CurvedBlock K = zero_block();
    outer(i, j, kRef, 0.5, 0, K);
    return K;
}

CurvedBlock CurvedKernelEntries::lambda_sauter_schwab(std::size_t i, std::size_t j, Adjacency adj) const {
    const auto& T = m_.flat.T;
    std::array<int, 3> ti = T[i], tj = T[j];
    if (adj == Adjacency::Edge || adj == Adjacency::Vertex) {             // gemeinsame Ecken zuerst (wie KernelEntries)
        std::array<int, 3> oi{}, oj{}; int ns = 0;
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) if (ti[a] == tj[b]) { oi[ns] = ti[a]; oj[ns] = tj[b]; ++ns; }
        int ki = ns, kj = ns;
        for (int a = 0; a < 3; ++a) { bool sh = false; for (int q = 0; q < ns; ++q) sh |= (ti[a] == oi[q]); if (!sh) oi[ki++] = ti[a]; }
        for (int b = 0; b < 3; ++b) { bool sh = false; for (int q = 0; q < ns; ++q) sh |= (tj[b] == oj[q]); if (!sh) oj[kj++] = tj[b]; }
        ti = oi; tj = oj;
        if (adj == Adjacency::Edge) {
            const auto& P = m_.flat.P;
            auto ang = [&](int a, int b, int c) { Vec3 u = P[b] - P[a], v = P[c] - P[a]; return std::acos(dot(u, v) / (norm(u) * norm(v))); };
            real keep = ang(tj[0], tj[1], tj[2]) - ang(ti[0], ti[1], ti[2]);
            real flip = ang(tj[1], tj[0], tj[2]) - ang(ti[1], ti[0], ti[2]);
            if (flip > keep) { std::swap(ti[0], ti[1]); std::swap(tj[0], tj[1]); }
        }
    }
    // Position der umgeordneten Ecke k in der Originalreihenfolge des Elements
    int perm_i[3], perm_j[3];
    for (int k = 0; k < 3; ++k) for (int a = 0; a < 3; ++a) { if (T[i][a] == ti[k]) perm_i[k] = a; if (T[j][a] == tj[k]) perm_j[k] = a; }
    const PairRule& R = ss_[static_cast<int>(adj)];
    CurvedBlock K = zero_block();
    const bool coinc = (adj == Adjacency::Coincident);
    for (std::size_t q = 0; q < R.w.size(); ++q) {
        const real lpx[3] = {1 - R.x[q][0] - R.x[q][1], R.x[q][0], R.x[q][1]};
        const real lpy[3] = {1 - R.y[q][0] - R.y[q][1], R.y[q][0], R.y[q][1]};
        std::array<real, 3> lx{}, ly{};
        for (int k = 0; k < 3; ++k) { lx[perm_i[k]] = lpx[k]; ly[perm_j[k]] = lpy[k]; }
        Vec3 nx, ny;
        const real Jx = gJ(i, lx, &nx), Jy = gJ(j, ly, &ny);
        const Vec3 x = gX(i, lx), y = gX(j, ly), z = x - y;
        const real r = norm(z);
        if (r < 1e-14) continue;
        const KernelValue kv = dirac_kernel_fast(z, k_);
        const real w = R.w[q] * Jx * Jy;
        if (!coinc) { acc(K, lx, ly, w, comps7(z, ny, kv.s, kv.vcoef)); continue; }
        // Selbstterm: Rest (ohne Phi_0) normal; Phi_0 = z/(4 pi r^3) n(y) = K_a + K_s
        const cplx p0 = 1.0 / (4 * pi * r * r * r);
        acc(K, lx, ly, w, comps7(z, ny, kv.s, kv.vcoef - p0));
        const Vec3 na = (nx + ny) * 0.5, ns = (ny - nx) * 0.5;
        const CurvedComp Ka = comps7(z, na, 0.0, p0), Ks = comps7(z, ns, 0.0, p0);
        acc(K, lx, ly, w, Ks);                                            // schwach singulaer
        for (int a = 0; a < 3; ++a)                                        // Hauptwert: antisymmetrisierte Gewichte
            for (int b = 0; b < 3; ++b) {
                const real f = 0.5 * w * (lx[a] * ly[b] - ly[a] * lx[b]);
                for (int c = 0; c < kCurvedComps; ++c) K[a * 3 + b][c] += f * Ka[c];
            }
    }
    return K;
}

}  // namespace cbem
