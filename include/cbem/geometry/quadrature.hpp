#pragma once
// Dreiecksquadratur: Dunavant-Regel Grad 5 (7 Punkte), optional auf sub^2 Teildreiecke verfeinert.
#include <vector>
#include "cbem/geometry/mesh.hpp"

namespace cbem {

struct QuadRule {                      // Regel auf dem Referenzdreieck (baryzentrisch)
    std::vector<std::array<real, 3>> bary;
    std::vector<real> w;               // Summe der Gewichte = 1
    static QuadRule dunavant7();
    static QuadRule subdivided(int sub);   // sub^2 Teildreiecke mit je 7 Punkten
};

// Quadraturpunkte aller Dreiecke eines Netzes (Gewichte inkl. Flaeche), flach gespeichert.
struct MeshQuadrature {
    int q = 0;                          // Punkte je Dreieck
    std::vector<Vec3> x;                // x[t*q + p]
    std::vector<real> w;
    MeshQuadrature() = default;
    MeshQuadrature(const TriangleMesh& m, const QuadRule& r);
    const Vec3* points(std::size_t t) const { return &x[t * q]; }
    const real* weights(std::size_t t) const { return &w[t * q]; }
};

}  // namespace cbem
