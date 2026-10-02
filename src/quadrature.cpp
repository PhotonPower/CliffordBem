#include "cbem/geometry/quadrature.hpp"
#include "cbem/geometry/sauter_schwab.hpp"

namespace cbem {

QuadRule QuadRule::dunavant7() {
    const real a1 = 0.059715871789770, b1 = 0.470142064105115;
    const real a2 = 0.797426985353087, b2 = 0.101286507323456;
    const real w1 = 0.132394152788506, w2 = 0.125939180544827;
    QuadRule r;
    r.bary = {{1. / 3, 1. / 3, 1. / 3}, {a1, b1, b1}, {b1, a1, b1}, {b1, b1, a1}, {a2, b2, b2}, {b2, a2, b2}, {b2, b2, a2}};
    r.w = {0.225, w1, w1, w1, w2, w2, w2};
    return r;
}

QuadRule QuadRule::conical(int n) {
    // lambda_1 = s, lambda_2 = t (1 - s), dA = (1 - s) ds dt; Gewichte auf Summe 1 (Referenzflaeche 1/2)
    std::vector<real> g, gw; gauss_legendre01(n, g, gw);
    QuadRule r;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
            const real l1 = g[a], l2 = g[b] * (1 - g[a]);
            r.bary.push_back({1 - l1 - l2, l1, l2}); r.w.push_back(2 * gw[a] * gw[b] * (1 - g[a]));
        }
    return r;
}

QuadRule QuadRule::subdivided(int sub) {
    if (sub <= 1) return dunavant7();
    QuadRule base = dunavant7(), r;
    auto add = [&](const std::array<std::array<real, 2>, 3>& tri) {
        for (std::size_t p = 0; p < base.w.size(); ++p) {
            real u = 0, v = 0;
            for (int k = 0; k < 3; ++k) { u += base.bary[p][k] * tri[k][0]; v += base.bary[p][k] * tri[k][1]; }
            r.bary.push_back({1 - u - v, u, v}); r.w.push_back(base.w[p] / (sub * sub));
        }
    };
    for (int i = 0; i < sub; ++i)
        for (int j = 0; j < sub - i; ++j) {
            real s = 1.0 / sub;
            add({{{i * s, j * s}, {(i + 1) * s, j * s}, {i * s, (j + 1) * s}}});
            if (j < sub - i - 1) add({{{(i + 1) * s, j * s}, {(i + 1) * s, (j + 1) * s}, {i * s, (j + 1) * s}}});
        }
    return r;
}

MeshQuadrature::MeshQuadrature(const TriangleMesh& m, const QuadRule& r) {
    q = static_cast<int>(r.w.size());
    x.resize(m.size() * q); w.resize(m.size() * q);
    for (std::size_t t = 0; t < m.size(); ++t) {
        auto v = m.vertices(t);
        for (int p = 0; p < q; ++p) {
            x[t * q + p] = v[0] * r.bary[p][0] + v[1] * r.bary[p][1] + v[2] * r.bary[p][2];
            w[t * q + p] = r.w[p] * m.area[t];
        }
    }
}

}  // namespace cbem
