#pragma once
// Analytische Dreiecksintegrale (Wilton-Formeln):
//   I_grad(x) = int_T (x - y)/|x - y|^3 dS_y,   I_inv(x) = int_T 1/|x - y| dS_y.
#include <array>
#include "cbem/core/types.hpp"

namespace cbem {

void triangle_integrals(const Vec3& x, const std::array<Vec3, 3>& p, const Vec3& n, Vec3& Igrad, real& Iinv);

// Dasselbe fuer die linearen Formfunktionen lambda_b (baryzentrische Koordinaten zur Ecke p[b]; v0.45):
//   Igrad[b] = int_T lambda_b(y) (x - y)/|x - y|^3 dS_y,   Iinv[b] = int_T lambda_b(y)/|x - y| dS_y.
// Mit dem Fusspunkt x* von x in der Ebene, rho = y - x*, R^2 = |rho|^2 + w^2 und lambda_b(y) = lambda_b(x*) + g_b.rho:
//   int lambda_b / R          = lambda_b(x*) Iinv + g_b . V1,                    V1 = int rho / R
//   int lambda_b (x - y)/R^3  = lambda_b(x*) Igrad + w n (g_b . V3) - M3 g_b,    V3 = int rho / R^3, M3 = int rho rho^T / R^3,
// ueber den Satz von Gauss in der Ebene als Kantenintegrale (m_e: aeussere Kantennormale in der Ebene, l: Bogenlaenge):
//   V1 = sum_e m_e int_e R dl,   V3 = -sum_e m_e int_e dl/R,   M3 = Iinv (1 - n n^T) - sum_e m_e (x) int_e rho/R dl.
// Summe ueber b gleich triangle_integrals (sum_b lambda_b = 1).
void triangle_integrals_linear(const Vec3& x, const std::array<Vec3, 3>& p, const Vec3& n, std::array<Vec3, 3>& Igrad,
                               std::array<real, 3>& Iinv);

}  // namespace cbem
