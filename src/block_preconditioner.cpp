#include "cbem/solvers/block_preconditioner.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include "cbem/operators/dense_blocks.hpp"
#include "cbem/hmatrix/cluster_tree.hpp"

namespace cbem {

FeatureSet FeatureSet::cube(real a, const Vec3& c) {
    FeatureSet f;
    for (int x : {-1, 1}) for (int y : {-1, 1}) for (int z : {-1, 1}) f.vertices.push_back(c + Vec3(a * x, a * y, a * z));
    for (int ax = 0; ax < 3; ++ax)
        for (int s1 : {-1, 1}) for (int s2 : {-1, 1}) {
            Vec3 p, q; int u = (ax + 1) % 3, v = (ax + 2) % 3;
            p[ax] = -a; q[ax] = a; p[u] = q[u] = a * s1; p[v] = q[v] = a * s2;
            f.edges.emplace_back(p + c, q + c);
        }
    return f;
}

static real seg_dist(const Vec3& x, const Vec3& p, const Vec3& q) {
    Vec3 d = q - p; real t = dot(x - p, d) / dot(d, d); t = std::max(0.0, std::min(1.0, t));
    return norm(x - (p + d * t));
}

std::vector<std::vector<std::size_t>> group_by_features(const TriangleMesh& m, const FeatureSet& f, real R) {
    const std::size_t nv = f.vertices.size();
    std::vector<std::vector<std::size_t>> G(nv + f.edges.size());
    for (std::size_t t = 0; t < m.size(); ++t) {
        const Vec3& c = m.centroid[t];
        std::size_t bv = 0; real dv = std::numeric_limits<real>::max();
        for (std::size_t v = 0; v < nv; ++v) { real d = norm(c - f.vertices[v]); if (d < dv) { dv = d; bv = v; } }
        if (dv < R) { G[bv].push_back(t); continue; }
        std::size_t be = 0; real de = std::numeric_limits<real>::max();
        for (std::size_t e = 0; e < f.edges.size(); ++e) { real d = seg_dist(c, f.edges[e].first, f.edges[e].second); if (d < de) { de = d; be = e; } }
        if (de < R) G[nv + be].push_back(t);
    }
    std::vector<std::vector<std::size_t>> out;
    for (auto& g : G) if (!g.empty()) out.push_back(std::move(g));
    return out;
}

std::vector<std::vector<std::size_t>> group_by_clusters(const TriangleMesh& m, std::size_t max_size) {
    ClusterTree tree(m, max_size);
    std::vector<std::vector<std::size_t>> out;
    for (const auto& c : tree.nodes) if (c.leaf()) out.push_back(tree.indices(c));
    return out;
}

BlockPreconditioner::BlockPreconditioner(const TriangleMesh& m, const KernelEntries& E1, const KernelEntries& E2,
                                         const TransmissionOperator& T, const std::vector<std::vector<std::size_t>>& groups)
    : T_(T), groups_(groups) {
    build(m, [&E1](const std::vector<std::size_t>& B) { return cauchy_block(E1, B, B); }, E2);
}

BlockPreconditioner::BlockPreconditioner(const TriangleMesh& m, const BlockFn& inner, const KernelEntries& E2,
                                         const TransmissionOperator& T, const std::vector<std::vector<std::size_t>>& groups)
    : T_(T), groups_(groups) { build(m, inner, E2); }

void BlockPreconditioner::build(const TriangleMesh& m, const BlockFn& inner, const KernelEntries& E2) {
    auto t0 = std::chrono::steady_clock::now();
    const auto& T = T_;
    const auto& J = T.J(); std::size_t tot = 0;
    lu_.resize(groups_.size()); piv_.resize(groups_.size());
#ifdef CBEM_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (long g = 0; g < static_cast<long>(groups_.size()); ++g) {
        const auto& B = groups_[g]; const std::size_t nb = 8 * B.size();
        Matrix A1 = inner(B), A2 = cauchy_block(E2, B, B);
        Matrix T_BB(nb, nb);
        // T = 1/2 (I + E2) + 1/2 (J - E1 J),  J blockdiagonal
        for (std::size_t bc = 0; bc < B.size(); ++bc) {
            const Mat8& Jc = J[B[bc]];
            for (std::size_t r = 0; r < nb; ++r)
                for (int q = 0; q < 8; ++q) {
                    cplx s = 0; for (int p = 0; p < 8; ++p) s += A1(r, 8 * bc + p) * Jc[p * 8 + q];
                    T_BB(r, 8 * bc + q) = 0.5 * A2(r, 8 * bc + q) - 0.5 * s;
                }
            for (int r = 0; r < 8; ++r) for (int q = 0; q < 8; ++q)
                T_BB(8 * bc + r, 8 * bc + q) += 0.5 * ((r == q ? 1.0 : 0.0) + Jc[r * 8 + q]);
        }
        lu_factor(T_BB, piv_[g]); lu_[g] = std::move(T_BB);
    }
    for (auto& B : groups_) { tot += B.size(); max_block_ = std::max(max_block_, 8 * B.size()); }
    frac_ = double(tot) / m.size();
    sec_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void BlockPreconditioner::apply(const std::vector<cplx>& x, std::vector<cplx>& y) const {
    T_.precondition(x, y);                                   // punktweise fuer alle
    std::vector<cplx> b;
    for (std::size_t g = 0; g < groups_.size(); ++g) {
        const auto& B = groups_[g]; b.resize(8 * B.size());
        for (std::size_t a = 0; a < B.size(); ++a) for (int q = 0; q < 8; ++q) b[8 * a + q] = x[8 * B[a] + q];
        lu_solve(lu_[g], piv_[g], b.data());
        for (std::size_t a = 0; a < B.size(); ++a) for (int q = 0; q < 8; ++q) y[8 * B[a] + q] = b[8 * a + q];
    }
}

}  // namespace cbem
