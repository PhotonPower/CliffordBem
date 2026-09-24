#pragma once
// Analytische Dreiecksintegrale (Wilton-Formeln):
//   I_grad(x) = int_T (x - y)/|x - y|^3 dS_y,   I_inv(x) = int_T 1/|x - y| dS_y.
#include <array>
#include "cbem/core/types.hpp"

namespace cbem {

void triangle_integrals(const Vec3& x, const std::array<Vec3, 3>& p, const Vec3& n, Vec3& Igrad, real& Iinv);

}  // namespace cbem
