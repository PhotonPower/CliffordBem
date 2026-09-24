#pragma once
// Resonanzfreie Transmissionsgleichung (AP 1): T_1 h = (E_2^+ + E_1^- J) h
//   = 1/2 (h + E_2 h) + 1/2 (J h - E_1 J h) = h_inc,
// h = Randspur des aeusseren Gesamtfelds (orthonormale stueckweise konstante Basis, 8 Komponenten je Dreieck).
// J: dreiecksweise Transmissionsabbildung mit Standardwahl der Hilfsparameter
//   a_H = sqrt(mu1)/sqrt(mu2) (Skalar), b_H = sqrt(eps1)/sqrt(eps2) (Pseudoskalar).
#include <vector>
#include "cbem/clifford/multivector.hpp"
#include "cbem/operators/cauchy_operator.hpp"

namespace cbem {

struct Medium { cplx eps{1.0}, mu{1.0}; };

// J als 8x8-Matrix zur Normalen n (Medium 1 innen, Medium 2 aussen)
Mat8 transmission_map(const Vec3& n, const Medium& inner, const Medium& outer);

class TransmissionOperator {
public:
    TransmissionOperator(const TriangleMesh& m, const CauchyOperator& E_inner, const CauchyOperator& E_outer,
                         const Medium& inner, const Medium& outer);
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    // punktweise Vorkonditionierung P = 2 (1 + J)^{-1}
    void precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    std::size_t size() const { return 8 * N_; }
    const std::vector<Mat8>& J() const { return J_; }
private:
    std::size_t N_;
    const CauchyOperator& E1_;
    const CauchyOperator& E2_;
    std::vector<Mat8> J_, P_;
};

}  // namespace cbem
