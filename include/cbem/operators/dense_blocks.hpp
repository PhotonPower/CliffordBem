#pragma once
// Dichte Teilbloecke von Randoperatoren mit exakten Eintraegen (fuer Vorkonditionierer und Referenzen).
#include <vector>
#include "cbem/assembly/kernel_entries.hpp"
#include "cbem/clifford/multivector.hpp"
#include "cbem/linalg/dense.hpp"

namespace cbem {

// E_{RC} (8|R| x 8|C|) des Cauchy-Operators in der orthonormalen Basis
Matrix cauchy_block(const KernelEntries& E, const std::vector<std::size_t>& R, const std::vector<std::size_t>& C);

}  // namespace cbem
