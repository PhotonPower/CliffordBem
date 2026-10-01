#pragma once
// Cauchy-Randoperator E_k fuer unstetig lineare Dichten (v0.45): Basis psi je Element orthonormal (LinearKernelEntries),
// Koeffizienten x[8 I + q], I = 3 t + a. Galerkin-Matrix
//   (E x)_I = -2 sum_J sum_c K_c(I, J) (e_c n_{t(J)}) x_J
// wie CauchyOperator, ohne Flaechenfaktoren (in psi enthalten). Ebene Elemente: Normale je Element konstant.
#include "cbem/operators/cauchy_operator.hpp"

namespace cbem {

class LinearCauchyOperator : public BoundaryOperator {
public:
    LinearCauchyOperator(const TriangleMesh& mesh, const KernelHMatrix& H);   // H aus LinearKernelEntries
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const override;   // x, y: 24 N
    std::size_t size() const override { return 8 * NI_; }
private:
    std::size_t NI_;                // 3 N Basisfunktionen
    const KernelHMatrix& H_;
    std::vector<Mat8> Ln_;          // NI*4 Matrizen L(e_c n_t)
};

}  // namespace cbem
