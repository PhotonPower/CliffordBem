#pragma once
// Optische Kraefte ueber den Maxwellschen Spannungstensor (v0.33).
//
// Zeitgemittelt, eps0 = mu0 = 1, achirales Aussenmedium (eps, mu):
//   <T> = 1/2 Re[ eps E (x) E* + mu H (x) H* - 1/2 (eps |E|^2 + mu |H|^2) I ],   F = int_S <T>.n dS  (n nach aussen)
// ueber eine geschlossene Flaeche S im Aussenraum, die genau den betrachteten Koerper umschliesst. Zwei Wege:
//   - Randspuren (force_from_traces): S = die Flaeche des Koerpers selbst, Felder aus den stueckweise konstanten Spuren
//     (Aussenseite); billig, aber dort sind die Felder am schaerfsten.
//   - umschliessende Flaeche (force_on_sphere, force_on_offset): Kugel um den Koerper bzw. Parallelflaeche im Abstand delta,
//     Felder aus dem Nahfeld (glatt); fuer einzelne Koerper eines Dimers die Parallelflaeche mit delta kleiner als der halbe
//     Spalt.
// Normierung: Die ebene Welle mit Amplitude |E0| traegt die Impulsstromdichte 1/2 eps |E0|^2; der Strahlungsdruck auf einen
// Koerper ist F_z = sigma_pr 1/2 eps |E0|^2, sigma_pr = sigma_ext - <cos theta> sigma_sca.
#include <vector>
#include "cbem/sources/near_field.hpp"

namespace cbem {

// <T>.n fuer Felder E, H und Normale n im Medium m
Vec3 stress_dot_normal(const CVec3& E, const CVec3& H, const Vec3& n, const Medium& m);

// Kraft je Koerper aus den Aussenspuren h auf der Vereinigung der Flaechen (body_begin: erster Dreiecksindex je Koerper,
// zuletzt die Gesamtzahl)
std::vector<Vec3> force_from_traces(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m,
                                    const std::vector<std::size_t>& body_begin);

// Kraft auf alles innerhalb der Kugel (Mittelpunkt c, Radius R), Gauss-Legendre in cos theta x gleichmaessig in phi, Felder aus
// dem Nahfeld der ebenen Welle (d, p)
Vec3 force_on_sphere(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d, const CVec3& p,
                     const Vec3& c, real R, int ntheta = 32, const NearFieldOptions& o = {});

// Kraft auf den von der Flaeche body umschlossenen Koerper: Parallelflaeche im Abstand delta, Mittelpunktsregel auf ihren
// Dreiecken (Felder an den Schwerpunkten), subdiv-fach unterteilt
Vec3 force_on_offset(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d, const CVec3& p,
                     const TriangleMesh& body, real delta, const NearFieldOptions& o = {});

}  // namespace cbem
