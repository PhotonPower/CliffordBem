#pragma once
// Streuproblem mit unstetig linearen Dichten auf ebenen Dreiecken (v0.45, Stufe 2a der gekruemmten Elemente).
// Gleiche Formulierung wie ScatteringProblem (T_1 = 1/2 (1 + E_2) + 1/2 (1 - E_1) J), Basis je Element orthonormal
// (psi = S lambda, LinearKernelEntries): 24 Unbekannte je Dreieck, Koeffizienten h[8 (3 t + a) + q].
// Umfang: ein oder mehrere Koerper (auch chiral), seit v0.59 auch chirales Aussenmedium (Helizitaetswellen mit k_pm,
// optisches Theorem je Kanal), ebene Wellen und beliebige rechte Seiten, Extinktion und Vorwaertsamplitude. Nahfeld, Kraefte
// und Block-/HODLR-Vorkonditionierung folgen.
// Ergebnis der Messung (docs/results_curved.md): auf ebenen Elementen allein kein Gewinn an der Kugel (der Geometriefehler
// dominiert); Grundlage fuer gekruemmte Elemente (Stufe 2b).
#include <memory>
#include <vector>

#include "cbem/assembly/linear_entries.hpp"
#include "cbem/operators/linear_cauchy_operator.hpp"
#include "cbem/operators/multibody_operator.hpp"
#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/problems/scattering_problem.hpp"

namespace cbem {

// L2-Projektion der ebenen Welle sqrt(eps)(p + I d x p) e^{ik d.x} auf die linearen Dichten (24 N)
std::vector<cplx> project_plane_wave_linear(const TriangleMesh& m, cplx k, cplx eps, const Vec3& d, const CVec3& p, int sub = 2);
// Fernfeld, Extinktion, Vorwaertsamplitude der Streuspur hs (24 N); wie far_field usw. fuer konstante Dichten
Multivector far_field_linear(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, const Vec3& xhat, int sub = 2);
real extinction_cross_section_linear(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p);
cplx forward_amplitude_linear(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p);
// Dichte (Multivektor) der Spur h (24 N) im Punkt x des Dreiecks t
Multivector linear_trace_value(const TriangleMesh& m, const std::vector<cplx>& h, std::size_t t, const Vec3& x);
// Mittelwert je Dreieck als konstante Spur (8 N, Normierung wie die konstanten Dichten): Vergleich mit ScatteringProblem
std::vector<cplx> linear_to_constant(const TriangleMesh& m, const std::vector<cplx>& h);

class LinearScatteringProblem {
public:
    LinearScatteringProblem(const std::vector<TriangleMesh>& bodies, const std::vector<Medium>& media, real omega, Medium outer = {},
                            HMatrixParams hp = {}, EntryParams ep = {});
    PlaneWaveResult solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o = {}) const;
    PlaneWaveResult solve_rhs(const std::vector<cplx>& b, const SolveOptions& o = {}) const;   // b: 24 N
    const TriangleMesh& mesh() const { return mb_.all; }
    const MultiBodyMesh& multibody() const { return mb_; }
    const TransmissionOperator& T() const { return *T_; }
    std::size_t unknowns() const { return T_->size(); }
    double hmatrix_bytes() const;
    std::size_t near_pairs() const;
private:
    MultiBodyMesh mb_;
    real omega_;
    Medium outer_;
    std::vector<std::unique_ptr<LinearKernelEntries>> ents_;
    std::vector<std::unique_ptr<KernelHMatrix>> hms_;
    std::vector<std::unique_ptr<LinearCauchyOperator>> cops_;
    std::vector<std::unique_ptr<ChiralCauchyOperator>> chops_;
    std::unique_ptr<BlockDiagonalOperator> inner_;
    const BoundaryOperator* inner_ptr_ = nullptr;
    const BoundaryOperator* outer_op_ = nullptr;
    std::unique_ptr<TransmissionOperator> T_;
};

}  // namespace cbem
