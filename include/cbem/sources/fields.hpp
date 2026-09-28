#pragma once
// Anregungen und Auswertung:
//  - L^2-Projektion einer ebenen Welle F = sqrt(eps)(p + I d x p) e^{ik d.x} auf die Basis,
//  - Fernfeld F_inf(xh) = -(ik/4pi) int e^{-ik xh.y} (1 + xh) n h^s dS  der aeusseren Streuspur h^s,
//  - Extinktionsquerschnitt ueber das optische Theorem, Vorwaertsamplitude (Betrag und Phase).
#include <vector>
#include "cbem/clifford/multivector.hpp"
#include "cbem/geometry/mesh.hpp"

namespace cbem {

std::vector<cplx> project_plane_wave(const TriangleMesh& m, cplx k, cplx eps, const Vec3& d, const Vec3& p, int sub = 2);
// komplexe Polarisation (z. B. zirkular: p = e1 + i s e2)
std::vector<cplx> project_plane_wave(const TriangleMesh& m, cplx k, cplx eps, const Vec3& d, const CVec3& p, int sub = 2);
// zirkulare Polarisation zur Ausbreitungsrichtung d: p = u + i s v mit (u, v, d) rechtshaendig
CVec3 circular_polarization(const Vec3& d, int s);
Multivector far_field(const TriangleMesh& m, const std::vector<cplx>& h_scat, cplx k, const Vec3& xhat, int sub = 2);
// sigma_ext = 4 pi / k Im( p . E_inf(d) ),  E_inf = Vektoranteil / sqrt(eps)
real extinction_cross_section(const TriangleMesh& m, const std::vector<cplx>& h_scat, cplx k, cplx eps, const Vec3& d, const Vec3& p);
// allgemein: sigma_ext = 4 pi / k Im( conj(p) . E_inf(d) ) / |p|^2
real extinction_cross_section(const TriangleMesh& m, const std::vector<cplx>& h_scat, cplx k, cplx eps, const Vec3& d, const CVec3& p);
// Vorwaertsamplitude S(0) = -i k conj(p).E_inf(d) / |p|^2 in der Normierung von Bohren-Huffman
// (E_s ~ e^{ikr}/(-ikr) S E_0, sigma_ext = 4 pi / k^2 Re S(0); Kugel: S(0) = sum (2n+1)/2 (a_n + b_n)).
// arg S(0) ist die Phase des vorwaerts gestreuten Lichts relativ zur einfallenden Welle.
cplx forward_amplitude(const TriangleMesh& m, const std::vector<cplx>& h_scat, cplx k, cplx eps, const Vec3& d, const CVec3& p);

}  // namespace cbem
