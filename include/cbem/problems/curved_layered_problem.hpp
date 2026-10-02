#pragma once
// Geschichtete Koerper auf gekruemmten (quadratischen) Elementen (v0.62): Kern-Schale, Mehrfachschichten, mehrere Koerper.
//
// Formulierung wie LayeredScatteringProblem (layered_problem.hpp): Grenzflaechengraph aus geschlossenen, nach aussen
// orientierten Flaechen, je Flaeche die Aussenspur h_s (24 Unbekannte je Element in der psi-Basis), Innenspur J_G h_s mit der
// Galerkin-Projektion der Transmissionsabbildung (Mediumpaar der Flaeche). Je Gebiet R ein CurvedCauchyOperator auf seinem
// Rand (Flaechen mit ihren eigenen Aussennormalen, Vorzeichen sigma in der Dichte); Zeile s:
//   1/2 (h_s - [E_out(sigma u)]_s) + 1/2 (J_G h_s - [E_in(sigma u)]_s) = b_s   (b nur an Flaechen zum Aussenraum).
// Fuer einen homogenen Koerper ist das T_1 von CurvedScatteringProblem. Vorkonditionierung 2 (1 + J_G)^{-1} je Element.
//
// Duenne Schichten (Flaechen im Abstand d << h): Nahquadratur mit Randabstand und Abstaenden zur gekruemmten Flaeche
// (CurvedNearParams::adapt_to_boundary, Voreinstellung curved_layered_near_params()). Goldkern mit Glasschale gegen
// Aden-Kerker, 2 x 1280 Elemente: sigma_ext auf 2e-6 (d = 0,2) bis 2,5e-5 (d = 0,01 ... 0,05), die Wirkung der Schicht direkt
// gegen die Rechnung ohne Schicht auf 2e-4 (konstante Dichten: 0,3-0,9 % bzw. 1-2 % gegen eine neutrale Vergleichsrechnung;
// docs/results_curved.md). Der Aufbau waechst mit abnehmendem d (Band der Breite d um die Kanten uebereinanderliegender
// Elemente): 17 s (d = 0,2) bis 440 s (d = 0,01) mit 12 Threads.
#include <memory>
#include <vector>

#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/curved_operators.hpp"
#include "cbem/problems/layered_problem.hpp"

namespace cbem {

struct CurvedLayeredGeometry {
    std::vector<QuadraticMesh> surfaces;        // geschlossene, nach aussen orientierte Flaechen
    std::vector<int> inside, outside;           // Gebiet innen / aussen je Flaeche
    std::vector<Medium> region_medium;          // Medium je Gebiet; Gebiet 0 = Aussenraum
    explicit CurvedLayeredGeometry(Medium exterior = {}) : region_medium{exterior} {}
    int add_region(const Medium& m) { region_medium.push_back(m); return static_cast<int>(region_medium.size()) - 1; }
    int add_surface(const QuadraticMesh& m, int in, int out) {
        surfaces.push_back(m); inside.push_back(in); outside.push_back(out); return static_cast<int>(surfaces.size()) - 1;
    }
    std::size_t elements() const { std::size_t n = 0; for (auto& s : surfaces) n += s.size(); return n; }
};

// wie fuer LayeredGeometry: homogener Koerper, verschachtelte Flaechen (von aussen nach innen), beschichteter Koerper
// (Parallelflaechen ueber offset_surface(QuadraticMesh), fuer glatte Flaechen)
int add_body(CurvedLayeredGeometry& g, const QuadraticMesh& surface, const Medium& inner, int parent = 0);
std::vector<int> add_layered_body(CurvedLayeredGeometry& g, const std::vector<QuadraticMesh>& surfaces_outer_first,
                                  const std::vector<Medium>& media, int parent = 0);
std::vector<int> add_coated_body(CurvedLayeredGeometry& g, const QuadraticMesh& surface, const Medium& core,
                                 const std::vector<Coating>& coatings, bool outward = true, int parent = 0);

// Nahquadratur fuer geschichtete Koerper: Randabstand und Abstaende zur gekruemmten Flaeche
inline CurvedNearParams curved_layered_near_params() { CurvedNearParams np; np.adapt_to_boundary = true; return np; }

class CurvedLayeredTransmissionOperator {
public:
    struct Region {
        const BoundaryOperator* E;               // Cauchy-Operator auf dem Gebietsrand (lokale Nummerierung)
        std::vector<std::size_t> el;             // globaler Elementindex je lokalem Element
        std::vector<char> is_inner;              // 1: Gebiet liegt innen an dieser Flaeche (Spur J_G h, sigma = +1)
    };
    // JG, P: je Element 24 x 24 (curved_transmission_blocks)
    CurvedLayeredTransmissionOperator(std::vector<Region> regions, std::vector<std::vector<cplx>> JG, std::vector<std::vector<cplx>> P);
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    void precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const;   // 2 (1 + J_G)^{-1} je Element
    void apply_J(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    std::size_t size() const { return 24 * JG_.size(); }
private:
    std::vector<Region> R_;
    std::vector<std::vector<cplx>> JG_, P_;
};

class CurvedLayeredScatteringProblem {
public:
    CurvedLayeredScatteringProblem(const CurvedLayeredGeometry& g, real omega, HMatrixParams hp = curved_hmatrix_params(), EntryParams ep = {},
                                   CurvedNearParams np = curved_layered_near_params());
    LayeredResult solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o = {}) const;
    const QuadraticMesh& mesh() const { return all_; }                 // alle Flaechen (globale Nummerierung)
    std::size_t surface_begin(std::size_t s) const { return begin_[s]; }
    const CurvedLayeredTransmissionOperator& T() const { return *T_; }
    std::size_t unknowns() const { return T_->size(); }
    double hmatrix_bytes() const;
    double near_seconds() const;
private:
    CurvedLayeredGeometry g_;
    real omega_;
    QuadraticMesh all_;
    std::vector<std::size_t> begin_;
    std::vector<std::unique_ptr<QuadraticMesh>> region_mesh_;
    std::vector<std::unique_ptr<CurvedKernelEntries>> ents_;
    std::vector<std::unique_ptr<CurvedHMatrix>> hms_;
    std::vector<std::unique_ptr<CurvedCauchyOperator>> cops_;
    std::vector<std::unique_ptr<ChiralCauchyOperator>> chops_;
    std::unique_ptr<CurvedLayeredTransmissionOperator> T_;
    std::vector<std::size_t> ext_el_;           // globale Elemente am Aussenraum (in der Reihenfolge von region_mesh_[0])
};

}  // namespace cbem
