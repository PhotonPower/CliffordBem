#pragma once
// Anregungen und Auswertung:
//  - L^2-Projektion einer ebenen Welle F = sqrt(eps)(p + I d x p) e^{ik d.x} auf die Basis,
//  - Fernfeld F_inf(xh) = -(ik/4pi) int e^{-ik xh.y} (1 + xh) n h^s dS  der aeusseren Streuspur h^s,
//  - Extinktionsquerschnitt ueber das optische Theorem.
#include <vector>
#include "cbem/clifford/multivector.hpp"
#include "cbem/geometry/mesh.hpp"

namespace cbem {

std::vector<cplx> project_plane_wave(const TriangleMesh& m, cplx k, cplx eps, const Vec3& d, const Vec3& p, int sub = 2);
Multivector far_field(const TriangleMesh& m, const std::vector<cplx>& h_scat, cplx k, const Vec3& xhat, int sub = 2);
// sigma_ext = 4 pi / k Im( p . E_inf(d) ),  E_inf = Vektoranteil / sqrt(eps)
real extinction_cross_section(const TriangleMesh& m, const std::vector<cplx>& h_scat, cplx k, cplx eps, const Vec3& d, const Vec3& p);

}  // namespace cbem
