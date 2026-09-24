#include "cbem/hmatrix/cluster_tree.hpp"
#include <algorithm>
#include <numeric>

namespace cbem {

ClusterTree::ClusterTree(const TriangleMesh& m, std::size_t leaf_size) {
    perm.resize(m.size()); std::iota(perm.begin(), perm.end(), 0);
    nodes.reserve(4 * m.size() / std::max<std::size_t>(1, leaf_size) + 4);
    build(m, 0, m.size(), leaf_size);
}

int ClusterTree::build(const TriangleMesh& m, std::size_t b, std::size_t e, std::size_t leaf_size) {
    int id = static_cast<int>(nodes.size()); nodes.emplace_back();
    ClusterNode c; c.begin = b; c.end = e;
    c.lo = m.centroid[perm[b]]; c.hi = c.lo; c.hmax = 0;
    for (std::size_t i = b; i < e; ++i) {
        const Vec3& p = m.centroid[perm[i]];
        for (int d = 0; d < 3; ++d) { c.lo[d] = std::min(c.lo[d], p[d]); c.hi[d] = std::max(c.hi[d], p[d]); }
        c.hmax = std::max(c.hmax, m.hmax[perm[i]]);
    }
    c.diam = norm(c.hi - c.lo);
    if (e - b > leaf_size) {
        Vec3 ext = c.hi - c.lo; int ax = 0; if (ext.y > ext[ax]) ax = 1; if (ext.z > ext[ax]) ax = 2;
        std::size_t mid = b + (e - b) / 2;
        std::nth_element(perm.begin() + b, perm.begin() + mid, perm.begin() + e,
                         [&](std::size_t u, std::size_t v) { return m.centroid[u][ax] < m.centroid[v][ax]; });
        if (mid > b && mid < e) {
            int l = build(m, b, mid, leaf_size); int r = build(m, mid, e, leaf_size);
            c.child[0] = l; c.child[1] = r;
        }
    }
    nodes[id] = c;
    return id;
}

real box_distance(const ClusterNode& a, const ClusterNode& b) {
    real s = 0;
    for (int d = 0; d < 3; ++d) {
        real g = std::max({0.0, a.lo[d] - b.hi[d], b.lo[d] - a.hi[d]}); s += g * g;
    }
    return std::sqrt(s);
}

}  // namespace cbem
