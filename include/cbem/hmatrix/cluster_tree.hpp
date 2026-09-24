#pragma once
// Geometrischer Clusterbaum ueber Dreiecksschwerpunkten (Halbierung entlang der laengsten Achse).
#include <vector>
#include "cbem/geometry/mesh.hpp"

namespace cbem {

struct ClusterNode {
    std::size_t begin = 0, end = 0;      // Bereich in ClusterTree::perm
    Vec3 lo, hi;                          // Bounding Box der Schwerpunkte
    real diam = 0, hmax = 0;              // Boxdiagonale, groesste Elementkante
    int child[2] = {-1, -1};
    bool leaf() const { return child[0] < 0; }
    std::size_t size() const { return end - begin; }
};

class ClusterTree {
public:
    ClusterTree(const TriangleMesh& m, std::size_t leaf_size);
    std::vector<std::size_t> perm;        // Dreiecksindizes, clusterweise zusammenhaengend
    std::vector<ClusterNode> nodes;       // nodes[0] = Wurzel
    std::vector<std::size_t> indices(const ClusterNode& c) const { return {perm.begin() + c.begin, perm.begin() + c.end}; }
private:
    int build(const TriangleMesh& m, std::size_t b, std::size_t e, std::size_t leaf_size);
};

real box_distance(const ClusterNode& a, const ClusterNode& b);

}  // namespace cbem
