#pragma once
// Gmsh-Netze (ASCII, Format 2.2 und 4.1): Oberflaechendreiecke (Elementtyp 2) werden nach physikalischer
// Gruppe (2.2) bzw. Entitaet/physikalischem Tag (4.1) in Koerper getrennt. Jeder Koerper wird ueber das
// Vorzeichen seines Volumens nach aussen orientiert (geschlossene, konsistent orientierte Flaechen).
// Zweite Ordnung (v0.48): 6-Knoten-Dreiecke (Elementtyp 9; Knoten 3, 4, 5 auf den Kanten 0-1, 1-2, 2-0) liest
// read_gmsh_quadratic als QuadraticMesh; read_gmsh liest davon nur die Ecken. Andere Elementtypen (Punkte, Linien,
// Volumenelemente, auch zweiter Ordnung) werden uebersprungen.
#include <string>
#include <vector>
#include "cbem/geometry/mesh.hpp"
#include "cbem/geometry/quadratic_mesh.hpp"

namespace cbem {

struct GmshBody { int tag; TriangleMesh mesh; };

std::vector<GmshBody> read_gmsh(const std::string& path, real scale = 1.0);
void write_gmsh22(const std::string& path, const std::vector<TriangleMesh>& bodies);   // physikalische Gruppen 1, 2, ...
real signed_volume(const TriangleMesh& m);

struct GmshQuadraticBody { int tag; QuadraticMesh mesh; };
// Netze zweiter Ordnung (Elementtyp 9); Fehler, wenn ein Koerper Dreiecke erster Ordnung enthaelt
std::vector<GmshQuadraticBody> read_gmsh_quadratic(const std::string& path, real scale = 1.0);
// Gmsh 2.2 zweiter Ordnung (Typ 9), physikalische Gruppen 1, 2, ...; Kantenmitten je Kante einmal (konform)
void write_gmsh22_quadratic(const std::string& path, const std::vector<QuadraticMesh>& bodies);

}  // namespace cbem
