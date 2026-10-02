#include "cbem/geometry/quadratic_mesh.hpp"

#include <cmath>
#include <stdexcept>

namespace cbem {

Vec3 QuadraticMesh::X(std::size_t t, const std::array<real, 3>& l) const {
    const auto v = flat.vertices(t); const auto& M = mid[t];
    return v[0] * (l[0] * (2 * l[0] - 1)) + v[1] * (l[1] * (2 * l[1] - 1)) + v[2] * (l[2] * (2 * l[2] - 1))
         + M[0] * (4 * l[0] * l[1]) + M[1] * (4 * l[1] * l[2]) + M[2] * (4 * l[2] * l[0]);
}

void QuadraticMesh::frame(std::size_t t, const std::array<real, 3>& l, Vec3& Xu, Vec3& Xv) const {
    // u = lambda_1, v = lambda_2, lambda_0 = 1 - u - v
    const auto v = flat.vertices(t); const auto& M = mid[t];
    Xu = v[1] * (4 * l[1] - 1) - v[0] * (4 * l[0] - 1) + M[0] * (4 * (l[0] - l[1])) + M[1] * (4 * l[2]) - M[2] * (4 * l[2]);
    Xv = v[2] * (4 * l[2] - 1) - v[0] * (4 * l[0] - 1) - M[0] * (4 * l[1]) + M[1] * (4 * l[1]) + M[2] * (4 * (l[0] - l[2]));
}

real QuadraticMesh::jacobian(std::size_t t, const std::array<real, 3>& l, Vec3* normal) const {
    Vec3 Xu, Xv; frame(t, l, Xu, Xv);
    const Vec3 c = cross(Xu, Xv); const real J = norm(c);
    if (normal) *normal = c / J;
    return J;
}

void QuadraticMesh::compute_geometry() {
    if (mid.size() != flat.size()) throw std::invalid_argument("QuadraticMesh: je Element drei Kantenmitten");
    const QuadRule R = QuadRule::subdivided(2);
    area.assign(size(), 0.0); bulge.assign(size(), 0.0);
    for (std::size_t t = 0; t < size(); ++t) {
        for (std::size_t p = 0; p < R.w.size(); ++p) {
            const real J = jacobian(t, {R.bary[p][0], R.bary[p][1], R.bary[p][2]});
            if (!(J > 0)) throw std::runtime_error("QuadraticMesh: entartetes oder gefaltetes Element");
            area[t] += 0.5 * R.w[p] * J;
        }
        const auto v = flat.vertices(t);
        for (int e = 0; e < 3; ++e) bulge[t] = std::max(bulge[t], norm(mid[t][e] - (v[e] + v[(e + 1) % 3]) * 0.5));
    }
}

QuadraticMesh make_quadratic(const TriangleMesh& flat, const std::function<Vec3(const Vec3&)>& project) {
    QuadraticMesh q; q.flat = flat; q.mid.resize(flat.size());
    for (std::size_t t = 0; t < flat.size(); ++t) {
        const auto v = flat.vertices(t);
        for (int e = 0; e < 3; ++e) q.mid[t][e] = project((v[e] + v[(e + 1) % 3]) * 0.5);
    }
    q.compute_geometry();
    return q;
}

QuadraticMesh make_quadratic_sphere(const TriangleMesh& flat, const Vec3& c, real r) {
    return make_quadratic(flat, [&](const Vec3& p) { const Vec3 d = p - c; return c + d * (r / norm(d)); });
}

QuadraticMesh quadratic_icosphere(int n, real radius) { return make_quadratic_sphere(make_icosphere(n, radius), Vec3(0, 0, 0), radius); }

QuadraticMesh translated(const QuadraticMesh& m, const Vec3& shift, real scale) {
    QuadraticMesh q; q.flat = translated(m.flat, shift, scale); q.mid = m.mid;
    for (auto& M : q.mid) for (auto& p : M) p = p * scale + shift;
    q.compute_geometry();
    return q;
}

QuadraticMesh merge_quadratic(const std::vector<QuadraticMesh>& parts, std::vector<std::size_t>* body_begin) {
    std::vector<TriangleMesh> flats; for (auto& p : parts) flats.push_back(p.flat);
    MultiBodyMesh mb = make_multibody(flats);
    QuadraticMesh q; q.flat = mb.all;
    for (auto& p : parts) q.mid.insert(q.mid.end(), p.mid.begin(), p.mid.end());
    q.compute_geometry();
    if (body_begin) *body_begin = mb.body_begin;
    return q;
}

CurvedQuadrature::CurvedQuadrature(const QuadraticMesh& m, const QuadRule& r) : q(static_cast<int>(r.w.size())) {
    const std::size_t N = m.size();
    x.resize(N * q); n.resize(N * q); w.resize(N * q); lam.resize(q);
    for (int p = 0; p < q; ++p) lam[p] = {r.bary[p][0], r.bary[p][1], r.bary[p][2]};
    for (std::size_t t = 0; t < N; ++t)
        for (int p = 0; p < q; ++p) {
            x[t * q + p] = m.X(t, lam[p]);
            w[t * q + p] = 0.5 * r.w[p] * m.jacobian(t, lam[p], &n[t * q + p]);
        }
}

real signed_volume(const QuadraticMesh& m, int sub) {
    CurvedQuadrature Q(m, QuadRule::subdivided(sub));
    real V = 0;
    for (std::size_t i = 0; i < Q.x.size(); ++i) V += Q.w[i] * dot(Q.x[i], Q.n[i]);
    return V / 3.0;
}

std::vector<std::array<real, 9>> curved_psi_matrices(const QuadraticMesh& m) {
    CurvedQuadrature Q(m, QuadRule::subdivided(2));
    std::vector<std::array<real, 9>> S(m.size());
    for (std::size_t t = 0; t < m.size(); ++t) {
        real G[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int p = 0; p < Q.q; ++p)
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) G[a][b] += Q.weights(t)[p] * Q.lam[p][a] * Q.lam[p][b];
        // Cholesky G = L L^T
        real L[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j <= i; ++j) {
                real s = G[i][j]; for (int k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
                if (i == j) { if (!(s > 0)) throw std::runtime_error("curved_psi_matrices: Gram-Matrix nicht positiv definit"); L[i][i] = std::sqrt(s); }
                else L[i][j] = s / L[j][j];
            }
        // S = L^{-1} (untere Dreiecksmatrix)
        real Si[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i = 0; i < 3; ++i) {
            Si[i][i] = 1.0 / L[i][i];
            for (int j = 0; j < i; ++j) { real s = 0; for (int k = j; k < i; ++k) s -= L[i][k] * Si[k][j]; Si[i][j] = s / L[i][i]; }
        }
        for (int a = 0; a < 3; ++a) for (int k = 0; k < 3; ++k) S[t][a * 3 + k] = Si[a][k];
    }
    return S;
}

}  // namespace cbem
