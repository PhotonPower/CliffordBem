#pragma once
// Cauchy-Randoperator E_k in der orthonormalen stueckweise konstanten Basis:
//   (E x)_i = -2/sqrt|tau_i| sum_c sum_j K_c(i,j) z_{j,c},   z_{j,c} = L(e_c n_j) x_j / sqrt|tau_j|
// (e_0 = 1, e_1..3 = Basisvektoren). Die Kernkomponenten liefert eine KernelHMatrix.
#include <vector>
#include "cbem/clifford/multivector.hpp"
#include "cbem/hmatrix/hmatrix.hpp"

namespace cbem {

// Schnittstelle fuer Randoperatoren auf der orthonormalen stueckweise konstanten Basis (8 Komponenten je Dreieck)
class BoundaryOperator {
public:
    virtual ~BoundaryOperator() = default;
    virtual void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const = 0;
    virtual std::size_t size() const = 0;
};

class CauchyOperator : public BoundaryOperator {
public:
    CauchyOperator(const TriangleMesh& mesh, const KernelHMatrix& H);
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const override;   // x, y: 8N
    std::size_t size() const override { return 8 * N_; }
    // Z aus x (fuer Referenzrechnungen mit denselben Konventionen)
    void make_Z(const std::vector<cplx>& x, std::vector<cplx>& Z) const;
private:
    std::size_t N_;
    const KernelHMatrix& H_;
    std::vector<Mat8> Ln_;          // N*4 Matrizen L(e_c n_j)/sqrt|tau_j|
    std::vector<real> isqA_;
};

}  // namespace cbem
