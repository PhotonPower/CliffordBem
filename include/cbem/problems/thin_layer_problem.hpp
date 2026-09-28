#pragma once
// Duenne Schichten in erster Ordnung: eine Flaeche, T_1 mit effektiver Transmissionsabbildung J_eff = J + L_d.
//
// Eine Schicht (Medium c, Dicke d) auf der Referenzflaeche Gamma aendert die Spruenge der Felder an Gamma. Mit dem
// verdraengten Medium X (aussen: X = 2, Schicht waechst nach aussen; innen: X = 1, Schicht verdraengt den Kern)
// ergibt die Integration der Maxwell-Gleichungen ueber die Schicht (e^{-i omega t}, Sprung [f] = f_2 - f_1 an Gamma):
//   [E_t] = d ( grad_G (D_n (1/eps_c - 1/eps_X)) - i omega (mu_c - mu_X) n x H_t )
//   [H_t] = d ( grad_G (B_n (1/mu_c  - 1/mu_X )) + i omega (eps_c - eps_X) n x E_t )
//   [D_n] = -d (eps_c - eps_X) div_G E_t,     [B_n] = -d (mu_c - mu_X) div_G H_t      (+ O(d^2)).
// Kruemmungsterme heben sich in erster Ordnung heraus (E_t, D_n sind in erster Naeherung stetig). Die Innenspur ist
// J h + Delta(h); T_eff = E_2^+ + E_1^- J_eff. Das entspricht einer Greenschen Funktion der Schichtfolge in erster
// Ordnung: E_1 (L_d h) ist das Cauchy-Integral des Kerns Phi_k1 nach partieller Integration der Ableitungen auf die
// Dichte. Die Flaechenableitungen stueckweise konstanter Dichten sind lokale Operatoren ueber die drei Kantennachbarn:
// grad_G als Kleinste-Quadrate-Gradient in der Tangentialebene (exakt fuer lineare Funktionen; der Green-Gauss-
// Gradient mit Kantenmitteln ist auf Ikosaedernetzen inkonsistent), div_G in Flussform mit Kantenmitteln (konservativ,
// keine Nettoladung, Ordnung 1); die Cauchy-Operatoren bleiben unveraendert.
// Lage der Referenzflaeche: der Anteil f (inner_fraction) der Schicht liegt innerhalb (verdraengt Medium 1), 1 - f
// ausserhalb (verdraengt Medium 2); beide Anteile addieren sich in erster Ordnung. f = 0: Referenz = Kernoberflaeche
// (Schicht waechst nach aussen), f = 1: Referenz = Aussenflaeche, f = 1/2: Schichtmitte (kleinster Fehler O(d^2),
// docs/results_coated.md).
// Gueltig fuer d << Kruemmungsradius und d |k_t| << 1 (auf dem Netz: d < h); Fehler O(d^2). Achirale Medien.
#include <memory>
#include <vector>
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/problems/layered_problem.hpp"

namespace cbem {

// Finite-Volumen-Flaechenoperatoren auf einem geschlossenen, konformen Dreiecksnetz (Werte je Dreieck)
struct SurfaceFV {
    struct Edge { std::size_t nb; real len; Vec3 conormal; Vec3 nb_conormal; Vec3 lsq; };   // Nachbar, Laenge, Konormalen
    std::vector<std::array<Edge, 3>> edges;                                    // (nach aussen), LSQ-Gewicht fuer grad
    explicit SurfaceFV(const TriangleMesh& m);
    const TriangleMesh& mesh;
    // grad_G phi (tangential im Dreieck), div_G v (v tangential je Dreieck)
    void grad(const std::vector<cplx>& phi, std::vector<std::array<cplx, 3>>& g) const;
    void div(const std::vector<std::array<cplx, 3>>& v, std::vector<cplx>& dv) const;
};

class ThinLayerTransmissionOperator {
public:
    ThinLayerTransmissionOperator(const TriangleMesh& m, const BoundaryOperator& E_inner, const BoundaryOperator& E_outer,
                                  const Medium& inner, const Medium& outer, const std::vector<Coating>& coatings,
                                  real inner_fraction, real omega);
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    void precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const;   // 2 (1 + J)^{-1} (ohne Schichtanteil)
    void apply_Jeff(const std::vector<cplx>& x, std::vector<cplx>& y) const;     // J_eff x
    std::size_t size() const { return 8 * m_.size(); }
private:
    const TriangleMesh& m_;
    const BoundaryOperator& E1_;
    const BoundaryOperator& E2_;
    Medium in_, out_;
    SurfaceFV fv_;
    std::vector<Mat8> J_, P_;
    cplx aE_ = 0, bE_ = 0, aH_ = 0, bH_ = 0;   // d-gewichtete Kontraste (siehe oben)
    real omega_;
};

// Streuproblem: ein homogener Koerper mit duennen Schichten erster Ordnung (eine Flaeche)
class ThinLayerScatteringProblem {
public:
    ThinLayerScatteringProblem(const TriangleMesh& surface, const Medium& core, const std::vector<Coating>& coatings,
                               real omega, Medium outer = {}, real inner_fraction = 0.0, HMatrixParams hp = {}, EntryParams ep = {});
    LayeredResult solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o = {}) const;
    const TriangleMesh& mesh() const { return m_; }
    double hmatrix_bytes() const { return H1_->stats().bytes() + H2_->stats().bytes(); }
private:
    TriangleMesh m_;
    Medium core_, outer_;
    real omega_;
    std::unique_ptr<KernelEntries> K1_, K2_;
    std::unique_ptr<KernelHMatrix> H1_, H2_;
    std::unique_ptr<CauchyOperator> E1_, E2_;
    std::unique_ptr<ThinLayerTransmissionOperator> T_;
};

}  // namespace cbem
