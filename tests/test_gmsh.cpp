// Gmsh-Ein-/Ausgabe: Rundreise 2.2, Lesen einer 4.1-Datei, Orientierung, gleiche Streuloesung.
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/scattering_problem.hpp"
#include <cstdio>
#include <fstream>
#include "check.hpp"
using namespace cbem;
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
    REPORT();
}
