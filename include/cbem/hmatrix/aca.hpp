#pragma once
// Adaptive Kreuzapproximation mit partieller Pivotisierung und QR/SVD-Nachkompression.
#include <functional>
#include "cbem/linalg/dense.hpp"

namespace cbem {

struct LowRank {                       // A ~ U V^T  (U: m x r, V: n x r)
    Matrix U, V;
    std::size_t rank() const { return U.cols; }
    std::size_t storage() const { return U.a.size() + V.a.size(); }
};

using RowFn = std::function<void(std::size_t i, cplx* out)>;   // Zeile i (Laenge n)
using ColFn = std::function<void(std::size_t j, cplx* out)>;   // Spalte j (Laenge m)

LowRank aca_partial(const RowFn& row, const ColFn& col, std::size_t m, std::size_t n, real eps, std::size_t rmax = 300);
void recompress(LowRank& lr, real eps);

}  // namespace cbem
