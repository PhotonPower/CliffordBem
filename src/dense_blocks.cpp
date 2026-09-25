#include "cbem/operators/dense_blocks.hpp"
#include <cmath>

namespace cbem {

Matrix cauchy_block(const KernelEntries& E, const std::vector<std::size_t>& R, const std::vector<std::size_t>& C) {
    const TriangleMesh& m = E.mesh();
    Matrix B(8 * R.size(), 8 * C.size());
    const Multivector e[4] = {Multivector::blade(0), Multivector::blade(1), Multivector::blade(2), Multivector::blade(4)};
    std::vector<std::array<Mat8, 4>> L(C.size());
    for (std::size_t b = 0; b < C.size(); ++b) {
        Multivector n = Multivector::vector(m.normal[C[b]]);
        for (int c = 0; c < 4; ++c) L[b][c] = (e[c] * n).left_matrix();
    }
    for (std::size_t a = 0; a < R.size(); ++a)
        for (std::size_t b = 0; b < C.size(); ++b) {
            KernelComp K = E.exact(R[a], C[b]);
            const real s = -2.0 / std::sqrt(m.area[R[a]] * m.area[C[b]]);
            for (int r = 0; r < 8; ++r)
                for (int q = 0; q < 8; ++q) {
                    cplx v = 0; for (int c = 0; c < 4; ++c) v += K[c] * L[b][c][r * 8 + q];
                    B(8 * a + r, 8 * b + q) = s * v;
                }
        }
    return B;
}

}  // namespace cbem
