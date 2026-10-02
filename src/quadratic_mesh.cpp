#include "cbem/geometry/quadratic_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

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

QuadraticMesh offset_surface(const QuadraticMesh& m0, real d) {
    QuadraticMesh m = m0;
    const std::size_t N = m0.size();
    std::vector<Vec3> nv(m0.flat.P.size(), Vec3{});
    std::map<std::pair<int, int>, Vec3> nm;                       // Kantenmitten, Schluessel (kleinere, groessere Ecke)
    auto key = [&](std::size_t t, int e) { const int a = m0.flat.T[t][e], b = m0.flat.T[t][(e + 1) % 3]; return std::make_pair(std::min(a, b), std::max(a, b)); };
    for (std::size_t t = 0; t < N; ++t)
        for (int e = 0; e < 3; ++e) {
            std::array<real, 3> l{0, 0, 0}; l[e] = 1;
            Vec3 n; m0.jacobian(t, l, &n); nv[m0.flat.T[t][e]] += n;
            l = {0, 0, 0}; l[e] = 0.5; l[(e + 1) % 3] = 0.5;
            m0.jacobian(t, l, &n); nm[key(t, e)] += n;
        }
    for (std::size_t v = 0; v < nv.size(); ++v) { const real a = norm(nv[v]); if (a > 0) m.flat.P[v] += nv[v] * (d / a); }
    for (std::size_t t = 0; t < N; ++t)
        for (int e = 0; e < 3; ++e) { const Vec3& n = nm[key(t, e)]; m.mid[t][e] = m0.mid[t][e] + n * (d / norm(n)); }
    m.flat.compute_geometry();
    m.compute_geometry();
    for (std::size_t t = 0; t < N; ++t) {
        const std::array<real, 3> c{1.0 / 3, 1.0 / 3, 1.0 / 3};
        Vec3 n0, n1; m0.jacobian(t, c, &n0); m.jacobian(t, c, &n1);
        if (!(dot(n0, n1) > 0.5) || !(m.area[t] > 1e-3 * m0.area[t]))
            throw std::runtime_error("offset_surface: Element " + std::to_string(t) + " klappt um oder entartet (|d| zu gross)");
    }
    if (!(signed_volume(m) * signed_volume(m0) > 0)) throw std::runtime_error("offset_surface: Flaeche stuelpt sich um (|d| zu gross)");
    // Faltung oder Durchdringung (wie im ebenen Fall): Knoten naeher als 0,8 |d| an einem anderen Teil der Originalflaeche;
    // Abstaende zum Sehnennetz, daher um die groesste Woelbung verringert
    real bmax = 0; for (real b : m0.bulge) bmax = std::max(bmax, b);
    const real lim = 0.8 * std::abs(d) - 4.0 / 3.0 * bmax;
    if (d != 0.0 && lim > 0) {
        std::vector<Vec3> q = m.flat.P; for (auto& M : m.mid) for (auto& p : M) q.push_back(p);
        const std::vector<real> dist = distance_to_surface(m0.flat, q, lim);
        for (std::size_t i = 0; i < q.size(); ++i)
            if (dist[i] < lim * (1 - 1e-9))
                throw std::runtime_error("offset_surface: Parallelflaeche faltet oder durchdringt sich (Knoten " + std::to_string(i) + ")");
    }
    return m;
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
