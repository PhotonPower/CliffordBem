#pragma once
// Sauter-Schwab-Quadratur fuer singulaere Galerkin-Doppelintegrale auf Dreieckspaaren
// (S. A. Sauter, C. Schwab: Boundary Element Methods, Springer 2011, Kap. 5.2):
// identische Dreiecke, gemeinsame Kante, gemeinsame Ecke. Die Regeln bilden [0,1]^4 so auf
// tau x tau ab, dass die Singularitaet durch die Jacobi-Determinante kompensiert wird.
// Punkte in Koordinaten des Einheitsdreiecks {u, v >= 0, u + v <= 1}; Summe der Gewichte = 1/4.
// Konventionen: gemeinsame Kante = Referenzkante (0,0)-(1,0) in beiden Dreiecken (gleiche
// Richtung), gemeinsame Ecke = Referenzecke (0,0).
#include <array>
#include <vector>
#include "cbem/core/types.hpp"

namespace cbem {

enum class Adjacency { None, Vertex, Edge, Coincident };

struct PairRule {
    std::vector<std::array<real, 2>> x, y;   // Test- und Ansatzpunkte (u, v)
    std::vector<real> w;
    static PairRule sauter_schwab(Adjacency a, int order);
    // anisotrop (v0.53): Gauss-Punkte je Richtung (xi, eta1, eta2, eta3) der Sauter-Schwab-Abbildung
    static PairRule sauter_schwab(Adjacency a, const std::array<int, 4>& orders);
};

// Gauss-Legendre auf [0,1]
void gauss_legendre01(int n, std::vector<real>& x, std::vector<real>& w);

}  // namespace cbem
