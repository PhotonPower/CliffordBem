#pragma once
// Elektrischer Dipol als Anregung, Zerfallsraten eines Emitters vor Nanostrukturen (v0.26).
//
// Konventionen wie im Kern: e^{-i omega t}, omega = k0 (eps0 = mu0 = 1), F = sqrt(eps) E + I sqrt(mu) H, H so normiert, dass fuer
// die ebene Welle H = sqrt(eps/mu) d x E gilt (d. h. H ist das physikalische H in Einheiten mit eps0 = mu0 = 1). Dipol p am
// Ort r0 im (achiralen) Aussenmedium, Strom J = -i omega p delta(x - r0), G = e^{ikr}/(4 pi r), n = (x - r0)/r:
//   H = omega k (n x p) G (1 - 1/(ikr)),
//   E = (1/eps) [ k^2 (n x p) x n + (3 n (n.p) - p)(1/r^2 - ik/r) ] G.
// Zerfallsraten relativ zum freien Dipol im selben Medium:
//   gesamt:     gamma/gamma0 = 1 + 6 pi eps Im(conj(p).E_s(r0)) / (k^3 |p|^2)   (Streufeld am Ort des Dipols),
//   strahlend:  Leistung des Fernfeldes (Dipol + Streufeld) / Leistung des freien Dipols (Gauss-Legendre x gleichmaessig in phi),
//   nicht strahlend (Absorption im Koerper): gesamt - strahlend.
// Die induzierte Ladung ist auf einen Fleck der Groesse etwa d (Abstand zur Oberflaeche) konzentriert; das Netz muss dort
// Elemente deutlich kleiner als d haben. Die Projektion des Dipolfeldes auf die Dreiecke unterteilt nahe Dreiecke adaptiv
// (Feld ~ 1/r^3), das behebt den Fehler der rechten Seite, nicht den der Diskretisierung.
#include <vector>
#include "cbem/sources/near_field.hpp"

namespace cbem {

// E und H des Dipols p am Ort r0 im Medium m (achiral) am Punkt x
void dipole_field(const Medium& m, real omega, const Vec3& r0, const CVec3& p, const Vec3& x, CVec3& E, CVec3& H);

// Projektion der Spur F = sqrt(eps) E + I sqrt(mu) H des Dipolfeldes auf die stueckweise konstanten Dichten (wie
// project_plane_wave), nahe Dreiecke adaptiv unterteilt (sub ~ 4 h / Abstand, hoechstens max_sub)
std::vector<cplx> project_dipole(const TriangleMesh& mesh, const Medium& m, real omega, const Vec3& r0, const CVec3& p, int max_sub = 24);

struct DipoleRates {
    real total = 0, radiative = 0, nonradiative = 0;   // relativ zum freien Dipol
    real distance = 0;                                  // Abstand des Dipols vom Netz
    bool too_close = false;                             // Abstand < 2 % der Elementgroesse (Streufeld am Dipolort unzuverlaessig)
};

// Raten aus der geloesten Gesamtspur h auf der Flaeche zum Aussenraum und der projizierten Dipolspur b (h_s = h - b)
DipoleRates dipole_rates(const TriangleMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b, const Medium& m,
                         real omega, const Vec3& r0, const CVec3& p, int ntheta = 40);

// Fluoreszenzverstaerkung eines fest, aber zufaellig orientierten Emitters (v0.27): Anregung exc[a] = |E_a|^2/|E0|^2 bei der
// Anregungswellenlaenge, Raten rates[a] bei der Emission, jeweils fuer die Richtungen a = x, y, z; q0 intrinsische Quantenausbeute:
//   F/F0 = sum_a exc[a] q_a / q0,   q_a = gamma_rad,a / (gamma_tot,a + (1 - q0)/q0).
// Anregung und Quantenausbeute werden gemeinsam gemittelt (nicht das Produkt der Mittelwerte).
real fluorescence_enhancement(const real exc[3], const DipoleRates rates[3], real q0);

}  // namespace cbem
