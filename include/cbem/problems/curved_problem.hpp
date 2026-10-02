#pragma once
// Streuproblem auf quadratischen (gekruemmten) Elementen mit unstetig linearen Dichten (v0.47, Stufe 2b der gekruemmten
// Elemente): 24 Unbekannte je Element, h[8 (3 t + a) + q] mit der je Element orthonormierten Basis psi (curved_psi_matrices).
// Umfang wie LinearScatteringProblem: ein oder mehrere Koerper (auch chiral), achirales Aussenmedium, ebene Wellen und
// beliebige rechte Seiten, Extinktion und Vorwaertsamplitude. Prueffall: Kugel gegen die Vorhersage aus Stufe 1b
// (docs/results_curved.md).
#include <memory>
#include <vector>

#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/curved_operators.hpp"
#include "cbem/operators/multibody_operator.hpp"
#include "cbem/problems/scattering_problem.hpp"

namespace cbem {

std::vector<cplx> project_plane_wave_curved(const QuadraticMesh& m, cplx k, cplx eps, const Vec3& d, const CVec3& p, int sub = 2);
Multivector far_field_curved(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, const Vec3& xhat, int sub = 2);
real extinction_cross_section_curved(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p);
cplx forward_amplitude_curved(const QuadraticMesh& m, const std::vector<cplx>& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p);

class CurvedScatteringProblem {
public:
    CurvedScatteringProblem(const std::vector<QuadraticMesh>& bodies, const std::vector<Medium>& media, real omega, Medium outer = {},
                            HMatrixParams hp = curved_hmatrix_params(), EntryParams ep = {}, CurvedNearParams np = {});
    PlaneWaveResult solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o = {}) const;
    PlaneWaveResult solve_rhs(const std::vector<cplx>& b, const SolveOptions& o = {}) const;
    const QuadraticMesh& mesh() const { return all_; }
    const std::vector<std::size_t>& body_begin() const { return begin_; }
    const CurvedTransmissionOperator& T() const { return *T_; }
    std::size_t unknowns() const { return T_->size(); }
    double hmatrix_bytes() const;
    double near_seconds() const;
private:
    std::vector<QuadraticMesh> parts_;
    QuadraticMesh all_;
    std::vector<std::size_t> begin_;
    real omega_;
    Medium outer_;
    std::vector<std::unique_ptr<CurvedKernelEntries>> ents_;
    std::vector<std::unique_ptr<CurvedHMatrix>> hms_;
    std::vector<std::unique_ptr<CurvedCauchyOperator>> cops_;
    std::vector<std::unique_ptr<ChiralCauchyOperator>> chops_;
    std::unique_ptr<BlockDiagonalOperator> inner_;
    const BoundaryOperator* inner_ptr_ = nullptr;
    const BoundaryOperator* outer_op_ = nullptr;
    std::unique_ptr<CurvedTransmissionOperator> T_;
};

}  // namespace cbem
