#pragma once
// Quadratische (gekruemmte) Dreieckselemente (v0.47, Stufe 2b): 6 Knoten je Element wie Gmsh ElementOrder 2 -- Ecken
// V_0, V_1, V_2 (die des ebenen Netzes) und Kantenmitten M_01, M_12, M_20 auf der Flaeche. Parametrisierung ueber dem
// Referenzdreieck (u, v), lambda = (1 - u - v, u, v):
//   X(lambda) = sum_k lambda_k (2 lambda_k - 1) V_k + 4 lambda_0 lambda_1 M_01 + 4 lambda_1 lambda_2 M_12 + 4 lambda_2 lambda_0 M_20,
// Flaechenelement dS = |X_u x X_v| du dv, Normale n = X_u x X_v / |X_u x X_v| (orientiert wie das ebene Netz).
// Gemeinsame Kanten haengen nur von ihren Ecken und ihrer Mitte ab: das Netz ist konform, wenn benachbarte Elemente dieselbe
// Mitte verwenden (make_quadratic projiziert die Sehnenmitte, also fuer beide Nachbarn denselben Punkt).
// Liegen die Mitten auf den Sehnen, ist das Element eben (Gegenprobe gegen die ebenen Elemente).
#include <array>
#include <functional>
#include <vector>

#include "cbem/geometry/mesh.hpp"
#include "cbem/geometry/quadrature.hpp"

namespace cbem {

struct QuadraticMesh {
    TriangleMesh flat;                          // Ecken, Konnektivitaet, Sehnendreiecke (Clusterbaum, Nachbarschaft, h_max)
    std::vector<std::array<Vec3, 3>> mid;       // je Element M_01, M_12, M_20
    std::vector<real> area;                     // Flaeche der gekruemmten Elemente
    std::vector<real> bulge;                    // groesster Abstand der Kantenmitten von den Sehnenmitten
    std::size_t size() const { return flat.size(); }
    Vec3 X(std::size_t t, const std::array<real, 3>& lam) const;
    // Tangenten X_u, X_v und Jacobi-Determinante |X_u x X_v|, Normale
    void frame(std::size_t t, const std::array<real, 3>& lam, Vec3& Xu, Vec3& Xv) const;
    real jacobian(std::size_t t, const std::array<real, 3>& lam, Vec3* normal = nullptr) const;
    void compute_geometry();                    // area, bulge
};

// aus einem ebenen Netz: Kantenmitten = project(Sehnenmitte)
QuadraticMesh make_quadratic(const TriangleMesh& flat, const std::function<Vec3(const Vec3&)>& project);
// Kantenmitten radial auf die Kugel um center mit Radius r (Ikosaederkugel: wie Gmsh fuer eine Kugel)
QuadraticMesh make_quadratic_sphere(const TriangleMesh& flat, const Vec3& center, real radius);
QuadraticMesh quadratic_icosphere(int n, real radius = 1.0);
QuadraticMesh translated(const QuadraticMesh& m, const Vec3& shift, real scale = 1.0);
// mehrere Koerper hintereinander (wie make_multibody)
QuadraticMesh merge_quadratic(const std::vector<QuadraticMesh>& parts, std::vector<std::size_t>* body_begin = nullptr);
// Volumen mit Vorzeichen (Gauss, int x.n dS / 3)
real signed_volume(const QuadraticMesh& m, int sub = 3);
// Parallelflaeche im Abstand d (v0.62; d > 0 nach aussen) fuer glatte Flaechen: Ecken und Kantenmitten entlang der
// gemittelten Normalen der angrenzenden gekruemmten Elemente an diesem Knoten verschoben (Kugel: Radiusfehler etwa d h^2).
// Wirft wie offset_surface(TriangleMesh), wenn Elemente umklappen, die Flaeche sich umstuelpt oder faltet.
QuadraticMesh offset_surface(const QuadraticMesh& m, real d);

// Quadraturpunkte aller Elemente in Parameterkoordinaten (Gewichte inkl. |X_u x X_v| und Referenzflaeche 1/2)
struct CurvedQuadrature {
    int q = 0;
    std::vector<Vec3> x, n;                     // Punkte, Normalen
    std::vector<real> w;
    std::vector<std::array<real, 3>> lam;       // baryzentrische Parameter je Punkt (gleich fuer alle Elemente)
    CurvedQuadrature() = default;
    CurvedQuadrature(const QuadraticMesh& m, const QuadRule& r);
    const Vec3* points(std::size_t t) const { return &x[t * q]; }
    const Vec3* normals(std::size_t t) const { return &n[t * q]; }
    const real* weights(std::size_t t) const { return &w[t * q]; }
};

// je Element orthonormierte Basis psi_a = sum_k S[a * 3 + k] lambda_k bezueglich int psi_a psi_b dS (gekruemmt):
// S = L^{-1} mit der Cholesky-Zerlegung G = L L^T der Gram-Matrix der lambda (|X_u x X_v| ist kein Polynom: unterteilte
// Regel mit 4 x 7 Punkten; ebene Elemente: exakt)
std::vector<std::array<real, 9>> curved_psi_matrices(const QuadraticMesh& m);

}  // namespace cbem
