#pragma once
// Richtungsquadraturen auf der Einheitskugel fuer Orientierungsmittelung (Lebedev, exakt fuer Polynome
// vom Grad 3, 5, 7 mit 6, 14, 26 Punkten). Mittelung ueber Einfallsrichtungen bei festem Koerper ist
// aequivalent zur Mittelung ueber Koerperorientierungen bei fester Einfallsrichtung.
#include <vector>
#include "cbem/core/types.hpp"

namespace cbem {

struct Direction { Vec3 d; real w; };
std::vector<Direction> lebedev(int npoints);   // 1 (nur +z, Gewicht 1), 6, 14, 26

}  // namespace cbem
