#include "cbem/geometry/sauter_schwab.hpp"
#include <cmath>
#include <stdexcept>

namespace cbem {

void gauss_legendre01(int n, std::vector<real>& x, std::vector<real>& w) {
    x.resize(n); w.resize(n);
    for (int i = 0; i < n; ++i) {
        real z = std::cos(pi * (i + 0.75) / (n + 0.5)), dp = 0;
        for (int it = 0; it < 100; ++it) {
            real p0 = 1, p1 = z;
            for (int k = 2; k <= n; ++k) { real p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
            dp = n * (z * p1 - p0) / (z * z - 1);
            real dz = p1 / dp; z -= dz; if (std::abs(dz) < 1e-15) break;
        }
        x[i] = 0.5 * (1 - z); w[i] = 1.0 / ((1 - z * z) * dp * dp);   // auf [0,1]: w = 2/((1-z^2) p'^2) / 2
    }
}

PairRule PairRule::sauter_schwab(Adjacency adj, int order) {
    std::vector<real> g, gw; gauss_legendre01(order, g, gw);
    PairRule r;
    // Punkte zunaechst im Referenzdreieck {0 <= x2 <= x1 <= 1}; am Ende u = x1 - x2, v = x2
    auto add = [&](real a1, real a2, real b1, real b2, real w) { r.x.push_back({a1 - a2, a2}); r.y.push_back({b1 - b2, b2}); r.w.push_back(w); };
    for (int i0 = 0; i0 < order; ++i0) for (int i1 = 0; i1 < order; ++i1) for (int i2 = 0; i2 < order; ++i2) for (int i3 = 0; i3 < order; ++i3) {
        const real xi = g[i0], e1 = g[i1], e2 = g[i2], e3 = g[i3];
        const real W = gw[i0] * gw[i1] * gw[i2] * gw[i3];
        const real e12 = e1 * e2, e123 = e1 * e2 * e3;
        if (adj == Adjacency::Coincident) {
            const real w = W * xi * xi * xi * e1 * e1 * e2;
            add(xi, xi * (1 - e1 + e12), xi * (1 - e123), xi * (1 - e1), w);
            add(xi * (1 - e123), xi * (1 - e1), xi, xi * (1 - e1 + e12), w);
            add(xi, xi * (e1 - e12 + e123), xi * (1 - e12), xi * (e1 - e12), w);
            add(xi * (1 - e12), xi * (e1 - e12), xi, xi * (e1 - e12 + e123), w);
            add(xi * (1 - e123), xi * (e1 - e123), xi, xi * (e1 - e12), w);
            add(xi, xi * (e1 - e12), xi * (1 - e123), xi * (e1 - e123), w);
        } else if (adj == Adjacency::Edge) {
            const real w = W * xi * xi * xi * e1 * e1;
            add(xi, xi * e1 * e3, xi * (1 - e12), xi * e1 * (1 - e2), w);
            add(xi, xi * e1, xi * (1 - e123), xi * e12 * (1 - e3), w * e2);
            add(xi * (1 - e12), xi * e1 * (1 - e2), xi, xi * e123, w * e2);
            add(xi * (1 - e123), xi * e12 * (1 - e3), xi, xi * e1, w * e2);
            add(xi * (1 - e123), xi * e1 * (1 - e2 * e3), xi, xi * e12, w * e2);
        } else if (adj == Adjacency::Vertex) {
            const real w = W * xi * xi * xi * e2;
            add(xi, xi * e1, xi * e2, xi * e2 * e3, w);
            add(xi * e2, xi * e2 * e3, xi, xi * e1, w);
        } else throw std::invalid_argument("sauter_schwab: keine singulaere Nachbarschaft");
    }
    return r;
}

}  // namespace cbem
