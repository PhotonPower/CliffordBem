#pragma once
// Dreiecksnetze (ebene Dreiecke, nach aussen orientiert) und Netzgeneratoren.
#include <array>
#include <vector>
#include "cbem/core/types.hpp"

namespace cbem {

struct TriangleMesh {
    std::vector<Vec3> P;                    // Knoten
    std::vector<std::array<int, 3>> T;      // Dreiecke (gegen den Uhrzeigersinn um die Aussennormale)
    // abgeleitete Groessen (compute_geometry)
    std::vector<Vec3> normal, centroid;
    std::vector<real> area, hmax;
    std::size_t size() const { return T.size(); }
    std::array<Vec3, 3> vertices(std::size_t t) const { return {P[T[t][0]], P[T[t][1]], P[T[t][2]]}; }
    void compute_geometry();
    void orient_outward(const Vec3& center = {0, 0, 0});   // fuer sternfoermige Gebiete um center
};

TriangleMesh make_icosphere(int n, real radius = 1.0);        // 20 n^2 Dreiecke
TriangleMesh translated(const TriangleMesh& m, const Vec3& shift, real scale = 1.0);  // x -> scale x + shift

// Mehrere Koerper: Vereinigung der Netze; Dreiecke eines Koerpers liegen zusammenhaengend,
// body_begin[b] .. body_begin[b+1] sind die Dreiecke von Koerper b.
struct MultiBodyMesh {
    TriangleMesh all;
    std::vector<TriangleMesh> parts;
    std::vector<std::size_t> body_begin;
    std::size_t bodies() const { return parts.size(); }
};
MultiBodyMesh make_multibody(const std::vector<TriangleMesh>& parts);
TriangleMesh make_cube_graded(int L);                         // Wuerfel [-1,1]^3, zu den Kanten gradiert, 48 (L+1)^2 Dreiecke
TriangleMesh make_cube_uniform(int n);                        // Wuerfel [-1,1]^3, gleichmaessig, 12 n^2 Dreiecke

// Parallelflaeche im Abstand d (d > 0 nach aussen, d < 0 nach innen) mit gleicher Topologie. Jeder Knoten wird
// um delta verschoben mit n_f . delta = d fuer alle angrenzenden Flaechennormalen n_f (kleinste Quadrate,
// Winkelgewichte, Pseudoinverse): auf glatten Stuecken delta = d n (gemittelte Normale), an Kanten und Ecken
// "auf Gehrung", sodass jede ebene Seite genau um d verschoben wird (Wuerfel -> Wuerfel mit Kante + 2d).
// Wirft std::runtime_error, wenn ein Dreieck umklappt oder entartet (|d| zu gross gegen Kruemmungsradius/Netz).
TriangleMesh offset_surface(const TriangleMesh& m, real d);

}  // namespace cbem
