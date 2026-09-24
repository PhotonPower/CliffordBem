#pragma once
// Neugestartetes GMRES mit Rechtsvorkonditionierung (Givens-Rotationen, modifiziertes Gram-Schmidt).
#include <functional>
#include <vector>
#include "cbem/core/types.hpp"

namespace cbem {

using LinOp = std::function<void(const std::vector<cplx>&, std::vector<cplx>&)>;

struct GmresResult { int iterations = 0; real rel_residual = 0; bool converged = false; };

// Loest A x = b; M (optional) Rechtsvorkonditionierer: A M u = b, x = M u.
GmresResult gmres(const LinOp& A, const std::vector<cplx>& b, std::vector<cplx>& x, const LinOp* M = nullptr,
                  real tol = 1e-6, int restart = 200, int max_iter = 2000);

}  // namespace cbem
