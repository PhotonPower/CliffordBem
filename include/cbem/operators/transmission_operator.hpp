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

// Pasteur-Medium: D = eps E + i chi H, B = mu H - i chi E (chi = 0: achiral)
struct Medium {
    cplx eps{1.0}, mu{1.0}, chi{0.0};
    cplx k(real omega, int helicity = 0) const {          // helicity = 0: achiral; +-1: k_pm = omega (sqrt(eps mu) -+ chi)
        return omega * (std::sqrt(eps) * std::sqrt(mu) - real(helicity) * chi);
    }
};

// J als 8x8-Matrix zur Normalen n (Medium 1 innen, Medium 2 aussen). Tangentialanteile mit sqrt(eps1/eps2)
// bzw. sqrt(mu1/mu2); Normalanteile (e_n, h_n)_1 = D_1 C_1^{-1} C_2 D_2^{-1} (e_n, h_n)_2 mit
// C_j = [[eps_j, i chi_j], [-i chi_j, mu_j]], D_j = diag(sqrt(eps_j), sqrt(mu_j)) (AP 1, Lemma Jchiral).
Mat8 transmission_map(const Vec3& n, const Medium& inner, const Medium& outer);

class TransmissionOperator {
public:
    TransmissionOperator(const TriangleMesh& m, const BoundaryOperator& E_inner, const BoundaryOperator& E_outer,
                         const Medium& inner, const Medium& outer);
    // mehrere Koerper: Innenmedium je Dreieck
    TransmissionOperator(const TriangleMesh& m, const BoundaryOperator& E_inner, const BoundaryOperator& E_outer,
                         const std::vector<Medium>& inner_per_triangle, const Medium& outer);
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    // punktweise Vorkonditionierung P = 2 (1 + J)^{-1}
    void precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    std::size_t size() const { return 8 * N_; }
    const std::vector<Mat8>& J() const { return J_; }
private:
    void setup(const TriangleMesh& m, const std::vector<Medium>& in, const Medium& out);
    std::size_t N_;
    const BoundaryOperator& E1_;
    const BoundaryOperator& E2_;
    std::vector<Mat8> J_, P_;
};

}  // namespace cbem
