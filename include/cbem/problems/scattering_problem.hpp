#pragma once
// Streuproblem fuer einen oder mehrere Koerper in Vakuum (Aussenmedium beliebig achiral):
// baut alle Randoperatoren (aussen auf der Vereinigung, innen blockdiagonal je Koerper, chiral wo chi != 0),
// den Transmissionsoperator T_1 und loest fuer ebene Wellen. Zentrale Einstiegsklasse fuer Anwendungen.
#include <memory>
#include <vector>
#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/multibody_operator.hpp"
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/solvers/gmres.hpp"

namespace cbem {

struct SolveOptions { real tol = 1e-6; int restart = 300; int max_iter = 3000; };
struct PlaneWaveResult { real sigma_ext = 0; int iterations = 0; real residual = 0; std::vector<cplx> h; };

class ScatteringProblem {
public:
    // bodies: Teilnetze (getrennte, nach aussen orientierte geschlossene Flaechen), media: Innenmedium je Koerper
    ScatteringProblem(const std::vector<TriangleMesh>& bodies, const std::vector<Medium>& media, real omega,
                      Medium outer = {}, HMatrixParams hp = {}, EntryParams ep = {}, bool union_interior = false);
    PlaneWaveResult solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o = {}) const;
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
    const CauchyOperator* outer_op_ = nullptr;
    std::unique_ptr<TransmissionOperator> T_;
};

}  // namespace cbem
