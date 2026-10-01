#pragma once
// Optische Kraefte ueber den Maxwellschen Spannungstensor (v0.33).
//
// Zeitgemittelt, eps0 = mu0 = 1, Aussenmedium (eps, mu; seit v0.36 auch chiral, D = eps E + i chi H, B = mu H - i chi E,
// Minkowski-Tensor 1/2 Re[E (x) D* + H (x) B* - 1/2 (E.D* + H.B*) I]); achiral:
//   <T> = 1/2 Re[ eps E (x) E* + mu H (x) H* - 1/2 (eps |E|^2 + mu |H|^2) I ],   F = int_S <T>.n dS  (n nach aussen)
// ueber eine geschlossene Flaeche S im Aussenraum, die genau den betrachteten Koerper umschliesst. Zwei Wege:
//   - Randspuren (force_from_traces): S = die Flaeche des Koerpers selbst, Felder aus den stueckweise konstanten Spuren
//     (Aussenseite); billig, aber dort sind die Felder am schaerfsten.
//   - umschliessende Flaeche (force_on_sphere, force_on_offset): Kugel um den Koerper bzw. Parallelflaeche im Abstand delta,
//     Felder aus dem Nahfeld (glatt); fuer einzelne Koerper eines Dimers die Parallelflaeche mit delta kleiner als der halbe
//     Spalt.
// Normierung: Die ebene Welle mit Amplitude |E0| traegt die Impulsstromdichte 1/2 eps |E0|^2; der Strahlungsdruck auf einen
// Koerper ist F_z = sigma_pr 1/2 eps |E0|^2, sigma_pr = sigma_ext - <cos theta> sigma_sca.
#include <functional>
#include <memory>
#include <vector>
#include "cbem/sources/near_field.hpp"

namespace cbem {

// Kraefte bei Dipolanregung (v0.37). Die Kraft auf einen Koerper liefern force_on_offset / force_on_sphere mit einem
// Feldauswerter fuer das Dipolfeld. Dazu: Kraft auf den Emitter (nur das Streufeld, das Eigenfeld ist am Ort singulaer)
//   F_em = 1/2 Re[sum p_j grad E_s,j* + sum m_j grad H_s,j*] - (omega k^3 / 12 pi) Re(p x m*),
// und der abgestrahlte Impuls Pi = int 1/2 eps |E_inf|^2 x dOmega (freier Dipol plus Streufeld). Impulserhaltung:
// F_Koerper + F_em + Pi = 0. b = project_dipole(outer, m, omega, r0, p, 24, md).
class DipoleField;
Vec3 emitter_force(const TriangleMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b, const Medium& m, real omega,
                   const DipoleField& src, real delta);
Vec3 radiated_momentum(const TriangleMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b, const Medium& m, real omega,
                       const DipoleField& src, int ntheta = 32);

// Feldauswerter (v0.37): liefert zu Punkten das Gesamtfeld; damit arbeiten alle Kraftfunktionen mit beliebiger Anregung
using NearFieldEval = std::function<std::vector<NearFieldPoint>(const std::vector<Vec3>&)>;
// allgemeine Anregung: b = Projektion des einfallenden Feldes (inc->project(outer, m) bzw. project_dipole); outer und h werden
// per Referenz gehalten und muessen leben, solange der Auswerter benutzt wird
NearFieldEval make_near_field_eval(const TriangleMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b, const Medium& m,
                                   real omega, std::shared_ptr<const IncidentField> inc, const NearFieldOptions& o = {});
NearFieldEval make_plane_wave_eval(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d,
                                   const CVec3& p, const NearFieldOptions& o = {});
Vec3 force_on_sphere(const NearFieldEval& nf, const Medium& m, const Vec3& c, real R, int ntheta = 32);
Vec3 force_on_offset(const NearFieldEval& nf, const Medium& m, const TriangleMesh& body, real delta);

// <T>.n fuer Felder E, H und Normale n im Medium m
Vec3 stress_dot_normal(const CVec3& E, const CVec3& H, const Vec3& n, const Medium& m);

// Kraft je Koerper aus den Aussenspuren h auf der Vereinigung der Flaechen (body_begin: erster Dreiecksindex je Koerper,
// zuletzt die Gesamtzahl)
std::vector<Vec3> force_from_traces(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m,
                                    const std::vector<std::size_t>& body_begin);

// Kraft aus der Impulsbilanz im Fernfeld (v0.36), unabhaengig vom Spannungstensor: Der Impulsstrom einer Helizitaetswelle ist
// (k_s/omega) mal ihr Energiestrom 1/2 sqrt(eps/mu) |E|^2; F = c [k_inc |p|^2 sigma_ext d - sum_s' k_s' int |E_inf,s'|^2 x dOmega],
// c = 1/2 sqrt(eps/mu)/omega (achiral ein Kanal). Kraft auf alle Koerper zusammen; sigma_ext aus der Loesung.
Vec3 force_from_far_field(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d, const CVec3& p,
                          real sigma_ext, int ntheta = 32);

// Kraft auf alles innerhalb der Kugel (Mittelpunkt c, Radius R), Gauss-Legendre in cos theta x gleichmaessig in phi, Felder aus
// dem Nahfeld der ebenen Welle (d, p)
Vec3 force_on_sphere(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d, const CVec3& p,
                     const Vec3& c, real R, int ntheta = 32, const NearFieldOptions& o = {});

// Kraft auf den von der Flaeche body umschlossenen Koerper: Parallelflaeche im Abstand delta, Mittelpunktsregel auf ihren
// Dreiecken (Felder an den Schwerpunkten), subdiv-fach unterteilt
Vec3 force_on_offset(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d, const CVec3& p,
                     const TriangleMesh& body, real delta, const NearFieldOptions& o = {});

// ---------------------------------------------------------------------------------------------------------------------------
// Kraft auf ein kleines Teilchen in Dipolnaeherung (v0.33, Stufe 2). Duale Groessen e = sqrt(eps) E, h = sqrt(mu) H:
//   p / sqrt(eps) = A_e e + i A_c h,   m / sqrt(mu) = A_m h - i A_c e
// (A_e, A_m: elektrische und magnetische Polarisierbarkeit, fuer eine Kugel A_e = 6 pi i a_1 / k^3, A_m = 6 pi i b_1 / k^3;
// A_c: chiraler Anteil, (A_++ - A_--)/2 in der Helizitaetsbasis). Zeitgemittelte Kraft
//   F = 1/2 Re[ sum_j p_j grad E_j* + sum_j m_j grad H_j* ] - (omega k^3 / (12 pi)) Re(p x m*),
// der letzte Term ist der Rueckstoss durch die Interferenz von p und m (aus der Impulsbilanz im Fernfeld). Fuer die ebene Welle
// ergibt das exakt den Mie-Strahlungsdruck mit n = 1 einschliesslich Re(a_1 b_1*). Das Teilchen wirkt nicht auf das Feld zurueck.
struct DipolePolarizability { cplx Ae = 0, Am = 0, Ac = 0; };
// aus den Mie-Koeffizienten a_1, b_1 der (achiralen) Kugel im Medium mit Wellenzahl k
inline DipolePolarizability polarizability_from_mie(cplx a1, cplx b1, cplx k) {
    const cplx c = 6 * pi * cplx(0, 1) / (k * k * k); return {c * a1, c * b1, 0.0};
}
struct FieldGradient { CVec3 E{}, H{}; CVec3 dE[3], dH[3]; };   // dE[i][j] = d_i E_j
// allgemeine Anregung ueber einen Feldauswerter (v0.37)
std::vector<FieldGradient> fields_with_gradients(const NearFieldEval& nf, const std::vector<Vec3>& x, real delta);
// E, H und Gradienten an den Punkten x (zentrale Differenzen mit Schritt delta) aus der Aussenspur, ebene Welle (d, p)
std::vector<FieldGradient> fields_with_gradients(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega,
                                                 const Vec3& d, const CVec3& p, const std::vector<Vec3>& x, real delta,
                                                 const NearFieldOptions& o = {});
// Kraft auf das Dipolteilchen im Feld g
Vec3 dipole_particle_force(const FieldGradient& g, const DipolePolarizability& a, const Medium& m, real omega);

}  // namespace cbem
