#pragma once
// Gmsh-Netze (ASCII, Format 2.2 und 4.1): Oberflaechendreiecke (Elementtyp 2) werden nach physikalischer
// Gruppe (2.2) bzw. Entitaet/physikalischem Tag (4.1) in Koerper getrennt. Jeder Koerper wird ueber das
// Vorzeichen seines Volumens nach aussen orientiert (geschlossene, konsistent orientierte Flaechen).
#include <string>
#include <vector>
#include "cbem/geometry/mesh.hpp"

namespace cbem {

struct GmshBody { int tag; TriangleMesh mesh; };

std::vector<GmshBody> read_gmsh(const std::string& path, real scale = 1.0);
void write_gmsh22(const std::string& path, const std::vector<TriangleMesh>& bodies);   // physikalische Gruppen 1, 2, ...
real signed_volume(const TriangleMesh& m);

}  // namespace cbem
