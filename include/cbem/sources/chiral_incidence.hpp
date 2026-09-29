#pragma once
// Ebene Welle und Extinktion in einem (moeglicherweise chiralen) Aussenmedium (v0.21).
//
// Im Pasteur-Medium (D = eps E + i chi H, B = mu H - i chi E) sind nur Helizitaetswellen Eigenmoden: F = P_+ F + P_- F mit
// den zentralen Projektoren P_pm = (1 +- iI)/2 und (nabla - i k_pm) F_pm = 0, k_pm = omega (sqrt(eps mu) -+ chi). Die
// Wellenimpedanz eta = sqrt(mu/eps) ist fuer beide gleich, daher bleibt F_inc = sqrt(eps) (p + I d x p) e^{i k d.x} fuer eine
// zirkulare Polarisation p gueltig; sie liegt ganz im Bild eines Projektors P_sigma und laeuft mit k_sigma. Linear polarisiertes
// Licht ist keine Eigenmode (es dreht seine Polarisationsebene auf dem Weg) und wird abgelehnt.
// Das Streufeld strahlt in beiden Kanaelen (k_+, k_-). Mit der einfallenden Welle interferiert nur der Kanal derselben
// Helizitaet; der andere hat eine andere Wellenzahl, seine Kreuzterme mitteln sich im Fernfeld heraus. Optisches Theorem je
// Kanal: sigma_ext = (4 pi / k_sigma) Im conj(p).E_inf^sigma(d) / |p|^2 mit dem Fernfeld von P_sigma h_s bei k_sigma.
// Achirales Aussenmedium: unveraendert (k, kein Projektor).
#include <vector>
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/sources/fields.hpp"

namespace cbem {

struct PlaneWaveIncidence {
    cplx k;        // Wellenzahl der einfallenden Welle
    int proj = 0;  // Projektor der einfallenden Welle: +1 (P_+), -1 (P_-), 0 (achiral, kein Projektor)
};

// Wellenzahl und Helizitaetskanal einer ebenen Welle (Richtung d, Polarisation p) im Medium m; wirft std::invalid_argument,
// wenn m chiral und p keine Helizitaetswelle ist
PlaneWaveIncidence plane_wave_incidence(const Medium& m, real omega, const Vec3& d, const CVec3& p);

// P_sigma h je Dreieck (Koeffizienten, 8 je Dreieck); sigma = 0: unveraendert
std::vector<cplx> helicity_part(const std::vector<cplx>& h, int sigma);

// Extinktion und Vorwaertsamplitude einer Streuspur h_s im Aussenmedium m (optisches Theorem je Kanal)
real extinction_in_medium(const TriangleMesh& mesh, const std::vector<cplx>& hs, const Medium& m, const PlaneWaveIncidence& inc,
                          const Vec3& d, const CVec3& p);
cplx forward_amplitude_in_medium(const TriangleMesh& mesh, const std::vector<cplx>& hs, const Medium& m, const PlaneWaveIncidence& inc,
                                 const Vec3& d, const CVec3& p);

}  // namespace cbem
