#pragma once
// Innerer Cauchy-Operator eines chiralen (Pasteur-)Mediums (AP 1, Schritt 1.6):
//   F = F_+ + F_-,  F_pm = P_pm F,  P_pm = (1 +- iI)/2,  (nabla - i k_pm) F_pm = 0,
//   k_pm = omega (sqrt(eps) sqrt(mu) -+ chi).
// Da P_pm zentral sind, ist E_1 h = P_+ E_{k+} h + P_- E_{k-} h wieder eine Involution.
#include "cbem/operators/cauchy_operator.hpp"

namespace cbem {

Mat8 helicity_projector(int sign);   // P_pm als 8x8-Matrix (Links- = Rechtsmultiplikation)

class ChiralCauchyOperator : public BoundaryOperator {
public:
    // beliebige Randoperatoren mit 8 Komponenten je Basisfunktion (konstante oder lineare Dichten; v0.45)
    ChiralCauchyOperator(const BoundaryOperator& E_plus, const BoundaryOperator& E_minus) : Ep_(E_plus), Em_(E_minus),
        Pp_(helicity_projector(+1)), Pm_(helicity_projector(-1)) {}
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const override;
    std::size_t size() const override { return Ep_.size(); }
private:
    const BoundaryOperator& Ep_;
    const BoundaryOperator& Em_;
    Mat8 Pp_, Pm_;
};

}  // namespace cbem
