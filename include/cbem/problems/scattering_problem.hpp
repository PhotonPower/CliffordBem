#pragma once
// Streuproblem fuer einen oder mehrere Koerper in Vakuum (Aussenmedium beliebig achiral):
// baut alle Randoperatoren (aussen auf der Vereinigung, innen blockdiagonal je Koerper, chiral wo chi != 0),
// den Transmissionsoperator T_1 und loest fuer ebene Wellen. Zentrale Einstiegsklasse fuer Anwendungen.
#include <memory>
#include <vector>
#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/multibody_operator.hpp"
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/solvers/block_preconditioner.hpp"
#include "cbem/solvers/hodlr.hpp"
#include "cbem/solvers/gmres.hpp"

namespace cbem {

struct SolveOptions { real tol = 1e-6; int restart = 300; int max_iter = 3000; };
struct PlaneWaveResult { real sigma_ext = 0; cplx forward = 0; int iterations = 0; real residual = 0; std::vector<cplx> h; };   // forward: S(0), siehe fields.hpp

class ScatteringProblem {
public:
    // bodies: Teilnetze (getrennte, nach aussen orientierte geschlossene Flaechen), media: Innenmedium je Koerper
    ScatteringProblem(const std::vector<TriangleMesh>& bodies, const std::vector<Medium>& media, real omega,
                      Medium outer = {}, HMatrixParams hp = {}, EntryParams ep = {}, bool union_interior = false);
    PlaneWaveResult solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o = {}) const;
    // Innerer Operator E_1 auf B x B (blockdiagonal ueber Koerper, chiral: P+ E_{k+} + P- E_{k-}), dicht
    Matrix inner_block(const std::vector<std::size_t>& B) const;
    // Blockvorkonditionierung einschalten (Gruppen z. B. aus group_by_clusters / group_by_features)
    void use_block_preconditioner(const std::vector<std::vector<std::size_t>>& groups);
    const BlockPreconditioner* block_preconditioner() const { return prec_.get(); }
    // 8x8-Block T_ij des Gesamtsystems (exakte Eintraege) und hierarchische Faktorisierung als Vorkonditionierer
    Mat8 system_entry(std::size_t i, std::size_t j) const;
    void use_hodlr_preconditioner(HodlrParams p = {});
    const HodlrSolver* hodlr() const { return hodlr_.get(); }
    const TriangleMesh& mesh() const { return mb_.all; }
    const MultiBodyMesh& multibody() const { return mb_; }
    const TransmissionOperator& T() const { return *T_; }
    double hmatrix_bytes() const;
private:
    MultiBodyMesh mb_;
    real omega_; Medium outer_;
    std::vector<std::unique_ptr<KernelEntries>> ents_;
    std::vector<std::unique_ptr<KernelHMatrix>> hms_;
    std::vector<std::unique_ptr<CauchyOperator>> cops_;
    std::vector<std::unique_ptr<ChiralCauchyOperator>> chops_;
    std::unique_ptr<BoundaryOperator> inner_;
    const BoundaryOperator* inner_ptr_ = nullptr;
    const BoundaryOperator* outer_op_ = nullptr;           // chirales Aussenmedium: P+ E_{k+} + P- E_{k-}
    std::unique_ptr<TransmissionOperator> T_;
    std::unique_ptr<BlockPreconditioner> prec_;
    std::unique_ptr<HodlrSolver> hodlr_;
    std::vector<std::array<Mat8, 4>> Lcn_;          // L(e_c n_j) je Dreieck
    std::vector<std::size_t> body_of_;              // Koerper je Dreieck (Innenoperator)
    const KernelEntries* outer_entries_ = nullptr;
    // je Koerper: Teile des Innenoperators (Eintraege, Helizitaet 0/+1/-1) und Dreiecksbereich
    struct InnerPart { const KernelEntries* E; int helicity; };
    std::vector<std::vector<InnerPart>> inner_parts_;
    std::vector<InnerPart> outer_parts_;                   // Anteile des Aussenoperators (Einzeleintraege)
    std::vector<std::size_t> inner_begin_;
};

}  // namespace cbem
