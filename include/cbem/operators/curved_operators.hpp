#pragma once
// Operatoren fuer gekruemmte (quadratische) Elemente mit unstetig linearen Dichten (v0.47, Stufe 2b).
//  - CurvedCauchyOperator: (E x)_I = -2 sum_J sum_c K_c(I, J) e_c x_J mit den sieben Komponenten von CurvedKernelEntries
//    (die Normale steht im Integral).
//  - CurvedTransmissionOperator: T_1 = 1/2 (1 + E_2) + 1/2 (1 - E_1) J_G mit der Galerkin-Projektion J_G der
//    Transmissionsabbildung je Element, (J_G)_(ab) = int psi_a psi_b J(n(x)) dS (24 x 24; die Normale variiert im Element).
//    Die exakte Komposition E_1 J waere ein Kern mit 64 Komponenten; die Projektion aendert die Extinktion der Kugel um etwa
//    20 % des Diskretisierungsfehlers und konvergiert gleich schnell (Prototyp, docs/results_curved.md).
//    Vorkonditionierung 2 (1 + J_G)^{-1} je Element.
#include <array>
#include <vector>

#include "cbem/geometry/quadratic_mesh.hpp"
#include "cbem/hmatrix/curved_hmatrix.hpp"
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/operators/cauchy_operator.hpp"

namespace cbem {

class CurvedCauchyOperator : public BoundaryOperator {
public:
    CurvedCauchyOperator(const QuadraticMesh& mesh, const CurvedHMatrix& H);
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const override;
    std::size_t size() const override { return 8 * NI_; }
private:
    std::size_t NI_;
    const CurvedHMatrix& H_;
    Mat8 Lb_[kCurvedComps];
};

// J_G = Galerkin-Projektion der Transmissionsabbildung (24 x 24 je Element, zeilenweise) und P = 2 (1 + J_G)^{-1} fuer
// Innen- und Aussenmedium je Element (v0.62; auch fuer geschichtete Koerper, deren Flaechen verschiedene Mediumpaare haben)
void curved_transmission_blocks(const QuadraticMesh& m, const std::vector<std::array<real, 9>>& S, const std::vector<Medium>& inner,
                                const std::vector<Medium>& outer, std::vector<std::vector<cplx>>& JG, std::vector<std::vector<cplx>>& P);
// y = B x mit 24 x 24-Bloecken je Element
void curved_block_apply(const std::vector<std::vector<cplx>>& B, const std::vector<cplx>& x, std::vector<cplx>& y);

class CurvedTransmissionOperator {
public:
    CurvedTransmissionOperator(const QuadraticMesh& mesh, const std::vector<std::array<real, 9>>& S, const BoundaryOperator& E_inner,
                               const BoundaryOperator& E_outer, const std::vector<Medium>& inner_per_element, const Medium& outer);
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    void precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    void apply_J(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    std::size_t size() const { return 24 * N_; }
    const std::vector<std::vector<cplx>>& JG() const { return JG_; }   // je Element 24 x 24, zeilenweise
private:
    std::size_t N_;
    const BoundaryOperator& E1_;
    const BoundaryOperator& E2_;
    std::vector<std::vector<cplx>> JG_, P_;
};

}  // namespace cbem
