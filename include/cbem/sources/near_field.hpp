#pragma once
// Nahfeld im Aussenraum (v0.24).
//
// Das Streufeld an einem Punkt x ausserhalb aller Koerper ist das Cauchy-Integral der Streuspur h_s = h - h_inc auf der
// Flaeche zum Aussenraum:  F_s(x) = sign * sum_tau int_tau Phi_k(x - y) dS_y  n_tau u_tau  (u_tau Wert der stueckweise
// konstanten Dichte). Quadratur: fern die 7-Punkt-Regel; nahe (Abstand < 3 h) die Zerlegung Phi_k = Phi_0 + (-ik/(4 pi r))
// + Rest mit analytischen Dreiecksintegralen fuer die singulaeren Anteile (triangle_integrals) und der 7-Punkt-Regel fuer den
// glatten Rest -- genau bis dicht an die Oberflaeche. Chirales Aussenmedium: je Helizitaet mit k_pm und P_pm h_s.
// E = vec(F) / sqrt(eps), H aus dem Bivektoranteil (I sqrt(mu) H). Punkte innerhalb eines Koerpers (Windungszahl der
// Aussenflaeche > 1/2) werden markiert; dort liefert die Aussendarstellung nicht das physikalische Feld. Ebenso Punkte naeher als
// 2 % der mittleren Elementgroesse am Netz (Kantenspitzen der stueckweise konstanten Dichte; ab etwa 3 % der Elementgroesse
// trifft das Nahfeld die Mie-Loesung auf wenige Prozent).
#include <vector>
#include "cbem/sources/chiral_incidence.hpp"

namespace cbem {

struct NearFieldPoint {
    CVec3 E{}, H{};              // Gesamtfeld (einfallend + gestreut)
    bool inside = false;         // innerhalb eines Koerpers (Feld nicht definiert)
    bool too_close = false;      // naeher als 2 % der Elementgroesse am Netz: stueckweise konstante Dichten springen an den
                                 // Kanten, das Cauchy-Integral hat dort logarithmische Spitzen (Wert unzuverlaessig)
    real enhancement = 0;        // |E|^2 / |E_0|^2
    real chirality = 0;          // Im(conj(E).H) / |Im(conj(E_0).H_0)| (optische Chiralitaet relativ zur zirkularen ebenen Welle)
};

// Streufeld F_s an Punkten (Multivektor), Streuspur hs auf der Flaeche m, Wellenzahl k (ein Helizitaetskanal)
std::vector<Multivector> scattered_field(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, const std::vector<Vec3>& pts);

// Gesamtes Nahfeld im Aussenraum: outer = Flaeche(n) zum Aussenraum, h = Gesamtspur darauf (wie im Loeser), Aussenmedium m.
// Die chirale Normierung bezieht sich auf die zirkulare ebene Welle derselben Richtung im selben Medium (auch bei linearer p).
std::vector<NearFieldPoint> exterior_near_field(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega,
                                                const Vec3& d, const CVec3& p, const std::vector<Vec3>& pts);

}  // namespace cbem
