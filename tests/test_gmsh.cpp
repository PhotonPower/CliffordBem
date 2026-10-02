// Gmsh-Ein-/Ausgabe: Rundreise 2.2, Lesen einer 4.1-Datei, Orientierung, gleiche Streuloesung.
// Zweite Ordnung (v0.48): eine echte Gmsh-Datei (CAD-Kugel, Mesh.ElementOrder 2, mit Punkten, Linien zweiter Ordnung und
// parametrischen Koordinaten; data/meshes), Kantenmitten auf der Kugel, Orientierung, Volumen, Rundreise mit umgedrehter
// Flaeche, Ablehnung erster Ordnung, Streuung auf gekruemmten Elementen gegen Mie.
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/curved_problem.hpp"
#include "cbem/problems/scattering_problem.hpp"
#include <cstdio>
#include <fstream>
#include "check.hpp"
using namespace cbem;
static void test_second_order() {
    const std::string f = std::string(CBEM_MESH_DIR) + "/sphere_r1_order2_gmsh41.msh";
    auto qb = read_gmsh_quadratic(f);
    auto fb = read_gmsh(f);                                              // dieselbe Datei, nur die Ecken
    CHECK(qb.size() == 1 && fb.size() == 1 && qb[0].mesh.size() == 320 && fb[0].mesh.size() == 320, "Gmsh 2. Ordnung: Elementzahl");
    const QuadraticMesh& q = qb[0].mesh;
    double dm = 0, out = 1e9;
    for (std::size_t t = 0; t < q.size(); ++t) {
        for (const auto& M : q.mid[t]) dm = std::max(dm, std::abs(norm(M) - 1));
        out = std::min(out, dot(q.flat.normal[t], q.flat.centroid[t]));
    }
    const double vq = signed_volume(q) / (4 * pi / 3) - 1, vf = signed_volume(q.flat) / (4 * pi / 3) - 1;
    std::printf("  Gmsh 2. Ordnung: %zu Elemente, Kantenmitten |M| - 1 = %.1e, Volumen quadratisch %+.2e, eben %+.2e\n", q.size(), dm, vq, vf);
    CHECK(dm < 1e-12, "Kantenmitten liegen nicht auf der CAD-Kugel");
    CHECK(out > 0, "nicht nach aussen orientiert");
    CHECK(std::abs(vq) < 2e-4 && std::abs(vf) > 100 * std::abs(vq), "Volumen der quadratischen Elemente");
    // Rundreise mit zwei Koerpern, der zweite absichtlich innen orientiert (Mitten muessen mitgetauscht werden)
    QuadraticMesh a = quadratic_icosphere(3), b = translated(quadratic_icosphere(2), Vec3(4, 0, 0), 0.5), bf = b;
    for (std::size_t t = 0; t < bf.size(); ++t) { std::swap(bf.flat.T[t][1], bf.flat.T[t][2]); auto& M = bf.mid[t]; M = {M[2], M[1], M[0]}; }
    bf.flat.compute_geometry(); bf.compute_geometry();
    write_gmsh22_quadratic("/tmp/cbem_test_o2.msh", {a, bf});
    auto rb = read_gmsh_quadratic("/tmp/cbem_test_o2.msh");
    double dmid = 0, dvol = std::abs(signed_volume(rb[1].mesh) - signed_volume(b));
    for (std::size_t t = 0; t < a.size(); ++t) for (int e = 0; e < 3; ++e) dmid = std::max(dmid, norm(rb[0].mesh.mid[t][e] - a.mid[t][e]));
    std::printf("  Rundreise 2. Ordnung: Mitten %.1e, Volumen des umgedrehten Koerpers %.1e\n", dmid, dvol);
    CHECK(rb.size() == 2 && dmid < 1e-15 && dvol < 1e-13 && signed_volume(rb[1].mesh) > 0, "Rundreise zweiter Ordnung");
    bool rejected = false;
    try { read_gmsh_quadratic("/tmp/cbem_test.msh"); } catch (const std::runtime_error&) { rejected = true; }
    CHECK(rejected, "Netz erster Ordnung als quadratisch angenommen");
    // Streuung: Glas auf der Gmsh-Kugel (heute -5,83 % mit den Ecken derselben Datei)
    SolveOptions so; so.tol = 1e-9;
    CurvedScatteringProblem C({q}, {Medium{2.25}}, 1.0);
    const real e = C.solve_plane_wave(Vec3(0, 0, 1), CVec3{1, 0, 0}, so).sigma_ext / pi / 0.2150978 - 1;
    std::printf("  Glas auf der Gmsh-Kugel (320 gekruemmte Elemente): %+.4f %% gegen Mie\n", 100 * e);
    CHECK(std::abs(e) < 3e-4, "Streuung auf dem Gmsh-Netz zweiter Ordnung: %.5f", e);
}

int main() {
    TriangleMesh a = make_icosphere(3), b = translated(make_cube_uniform(3), Vec3(4, 0, 0), 0.8);
    // absichtlich eine Flaeche falsch orientieren: der Leser muss sie umdrehen
    TriangleMesh bf = b; for (auto& t : bf.T) std::swap(t[1], t[2]); bf.compute_geometry();
    write_gmsh22("/tmp/cbem_test.msh", {a, bf});
    auto bodies = read_gmsh("/tmp/cbem_test.msh");
    CHECK(bodies.size() == 2, "zwei Koerper erwartet (%zu)", bodies.size());
    CHECK(bodies[0].mesh.size() == a.size() && bodies[1].mesh.size() == b.size(), "Dreieckszahl");
    CHECK(signed_volume(bodies[0].mesh) > 0 && signed_volume(bodies[1].mesh) > 0, "Orientierung nicht nach aussen");
    std::printf("  Volumen: Kugel %.6f, Wuerfel %.6f (exakt %.6f)\n", signed_volume(bodies[0].mesh), signed_volume(bodies[1].mesh), 1.6 * 1.6 * 1.6);
    // Format 4.1 (von Hand, ein Tetraeder als geschlossene Flaeche, Entitaet 7 mit physikalischem Tag 3)
    std::ofstream o("/tmp/cbem_test41.msh");
    o << "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Entities\n0 0 1 0\n7 0 0 0 1 1 1 1 3 0\n$EndEntities\n"
         "$Nodes\n1 4 1 4\n2 7 0 4\n1\n2\n3\n4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n"
         "$Elements\n1 4 1 4\n2 7 2 4\n1 1 3 2\n2 1 2 4\n3 1 4 3\n4 2 3 4\n$EndElements\n";
    o.close();
    auto tet = read_gmsh("/tmp/cbem_test41.msh");
    CHECK(tet.size() == 1 && tet[0].tag == 3 && tet[0].mesh.size() == 4, "Format 4.1");
    CHECK(std::abs(signed_volume(tet[0].mesh) - 1.0 / 6) < 1e-14, "Tetraedervolumen %.3e", signed_volume(tet[0].mesh));
    // gleiche Streuloesung fuer das eingelesene Kugelnetz
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10; const Medium g{2.25, 1.0, 0.0};
    real s0 = ScatteringProblem({a}, {g}, 1.0, {}, hp).solve_plane_wave(Vec3(0, 0, 1), CVec3{1.0, 0.0, 0.0}, so).sigma_ext;
    real s1 = ScatteringProblem({bodies[0].mesh}, {g}, 1.0, {}, hp).solve_plane_wave(Vec3(0, 0, 1), CVec3{1.0, 0.0, 0.0}, so).sigma_ext;
    std::printf("  sigma_ext: erzeugt %.10f, eingelesen %.10f\n", s0, s1);
    CHECK(std::abs(s0 - s1) < 1e-9 * s0, "eingelesenes Netz liefert andere Loesung");
    test_second_order();
    REPORT();
}
