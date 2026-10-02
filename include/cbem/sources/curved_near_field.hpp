#pragma once
// Nahfeld und Kraefte auf gekruemmten Elementen (v0.58).
//
// Streufeld im Aussenraum aus der Streuspur h_s = h - b mit 24 Koeffizienten je Element (unstetig lineare Dichten in der je
// Element orthonormierten Basis psi_a):
//   F_s(x) = sum_tau sum_a [ int_tau psi_a(y) Phi_k(x - y) n(y) dS_y ] u_{tau,a},   u_{tau,a} = h_s[8 (3 tau + a) + 0..7],
// mit der Normalen im Integral (wie CurvedKernelEntries) und demselben Vorzeichen wie der Pfad konstanter Dichten
// (near_field.hpp). Quadratur: fern die 7-Punkt-Regel je Element, nah die Singularitaetssubtraktion der nahen Paare
// (CurvedKernelEntries::point_integrals). Direkte Summation O(Punkte x Elemente).
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

// Streufeld F_s an Punkten (Multivektor), Streuspur hs auf der gekruemmten Flaeche m, Wellenzahl k
std::vector<Multivector> scattered_field_curved(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, const std::vector<Vec3>& pts);

// Gesamtes Nahfeld im Aussenraum: h = Gesamtspur (wie im Loeser), b = Projektion des einfallenden Feldes
std::vector<NearFieldPoint> exterior_near_field_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b,
                                                       const Medium& m, real omega, const IncidentField& inc, const std::vector<Vec3>& pts);
// ebene Welle (d, p): b wie im Loeser (project_plane_wave_curved)
std::vector<NearFieldPoint> exterior_near_field_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega,
                                                       const Vec3& d, const CVec3& p, const std::vector<Vec3>& pts);

// Feldauswerter fuer die Kraftfunktionen; outer und h werden per Referenz gehalten
NearFieldEval make_near_field_eval_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b, const Medium& m,
                                          real omega, std::shared_ptr<const IncidentField> inc);
NearFieldEval make_plane_wave_eval_curved(const QuadraticMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega, const Vec3& d,
                                          const CVec3& p);

}  // namespace cbem
