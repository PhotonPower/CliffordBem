#pragma once
// Nahfeld und Kraefte auf gekruemmten Elementen (v0.58).
//
// Streufeld im Aussenraum aus der Streuspur h_s = h - b mit 24 Koeffizienten je Element (unstetig lineare Dichten in der je
// Element orthonormierten Basis psi_a):
//   F_s(x) = sum_tau sum_a [ int_tau psi_a(y) Phi_k(x - y) n(y) dS_y ] u_{tau,a},   u_{tau,a} = h_s[8 (3 tau + a) + 0..7],
// mit der Normalen im Integral (wie CurvedKernelEntries) und demselben Vorzeichen wie der Pfad konstanter Dichten
// (near_field.hpp). Quadratur: fern die 7-Punkt-Regel je Element, nah die Singularitaetssubtraktion der nahen Paare
// (CurvedKernelEntries::point_integrals). Direkte Summation O(Punkte x Elemente); ab 4000 Punkten H-Matrix (v0.61, siehe
// scattered_field_curved_hmatrix).
//
// Markierungen beziehen sich auf das Sehnennetz (outer.flat): innen (Windungszahl > 1/2) und zu nah (naeher als 2 % der
// Elementgroesse oder als die groesste Woelbung der Elemente am Sehnennetz; dort kann ein Punkt zwischen Sehne und gekruemmter
// Flaeche liegen). Chirales Aussenmedium (v0.59): Streufeld je Helizitaet mit k_pm und P_pm h_s (wie near_field.hpp).
//
// Kraefte: Die Feldauswerter make_near_field_eval_curved / make_plane_wave_eval_curved liefern das Gesamtfeld und arbeiten mit
// force_on_sphere, force_on_offset und fields_with_gradients aus optical_force.hpp.
#include <memory>
#include <vector>

#include "cbem/geometry/quadratic_mesh.hpp"
#include "cbem/sources/optical_force.hpp"

namespace cbem {

// Spur eines beliebigen einfallenden Feldes F = sqrt(eps) E + I sqrt(mu) H in der psi-Basis (24 je Element; unterteilte
// 7-Punkt-Regel mit sub^2 Teildreiecken auf der gekruemmten Flaeche; sub = 2 wie project_plane_wave_curved im Loeser)
std::vector<cplx> project_incident_curved(const QuadraticMesh& m, const IncidentField& inc, const Medium& med, int sub = 2);
// dieselbe Projektion in zwei Schritten (v0.60, fuer Felder aus Python): Quadraturpunkte (elementweise, N q) und Spur aus den
// Feldwerten E, H an diesen Punkten
std::vector<Vec3> projection_points_curved(const QuadraticMesh& m, int sub = 2);
std::vector<cplx> project_samples_curved(const QuadraticMesh& m, const std::vector<CVec3>& E, const std::vector<CVec3>& H, const Medium& med,
                                         int sub = 2);

// Optionen der Nahfeldauswertung gekruemmter Elemente (v0.61): H-Matrix ab 4000 Punkten (darunter, etwa fuer die 2048 Punkte
// von force_on_sphere, ist die direkte Summation schneller), ACA-Toleranz 1e-6 (Fehler etwa 3e-7 bezogen auf max |F_s|, weit
// unter dem Diskretisierungsfehler); eben 2000 Punkte und 1e-4
inline NearFieldOptions curved_near_field_options() { NearFieldOptions o; o.hmatrix_min_points = 4000; o.eps = 1e-6; return o; }

// Streufeld F_s an Punkten (Multivektor), Streuspur hs auf der gekruemmten Flaeche m, Wellenzahl k; direkte Summation
std::vector<Multivector> scattered_field_curved(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, const std::vector<Vec3>& pts);

// dasselbe Streufeld mit einer H-Matrix (v0.61): Zeilen = Punkte, Spalten = Elemente (Clusterbaeume mit Blattgroesse 64).
// Zulaessige Bloecke (min(diam) <= 2 dist, dist > 3 h_max, dist > 4 r_max, |k| diam <= 20) liegen ganz im Bereich der Fernregel
// von point_integrals; dort ist das Feld eine Summe ueber die 7 Quadraturpunkte je Element, F = sum_q (s + v z) W_q mit
// W_q = n_q sum_a w_q psi_a(y_q) u_a, und ACA+ approximiert den Dirac-Kern (s, v z) zwischen Punkten und Quadraturpunkten
// (4 Komponenten). Unzulaessige Bloecke direkt mit point_integrals. Jeder Block wird nach dem Aufbau sofort angewandt und
// verworfen (ohne Nachkompression; Speicher je Thread ein Block). 5120 Elemente x 40 000 Punkte, 12 Threads: 3,7 s statt
// 27 s direkt (eps = 1e-6, Fehler 3e-7); eps = 1e-4: 2,2 s, Fehler 4e-5. stats: Blockzahlen, Raenge, Zeit (nichts gespeichert).
std::vector<Multivector> scattered_field_curved_hmatrix(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, const std::vector<Vec3>& pts,
                                                        real eps = 1e-6, HStats* stats = nullptr);

// Gesamtes Nahfeld im Aussenraum: h = Gesamtspur (wie im Loeser), b = Projektion des einfallenden Feldes; H-Matrix ab
// opt.hmatrix_min_points Punkten mit ACA-Toleranz opt.eps
std::vector<NearFieldPoint> exterior_near_field_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b,
                                                       const Medium& m, real omega, const IncidentField& inc, const std::vector<Vec3>& pts,
                                                       const NearFieldOptions& opt = curved_near_field_options());
// ebene Welle (d, p): b wie im Loeser (project_plane_wave_curved)
std::vector<NearFieldPoint> exterior_near_field_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega,
                                                       const Vec3& d, const CVec3& p, const std::vector<Vec3>& pts,
                                                       const NearFieldOptions& opt = curved_near_field_options());

// Feldauswerter fuer die Kraftfunktionen; outer und h werden per Referenz gehalten
NearFieldEval make_near_field_eval_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b, const Medium& m,
                                          real omega, std::shared_ptr<const IncidentField> inc,
                                          const NearFieldOptions& opt = curved_near_field_options());
NearFieldEval make_plane_wave_eval_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d,
                                          const CVec3& p, const NearFieldOptions& opt = curved_near_field_options());

}  // namespace cbem
