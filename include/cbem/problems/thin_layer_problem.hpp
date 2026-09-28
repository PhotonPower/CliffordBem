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
//
// Dirac-Form (Modelle Dirac1, Dirac2): Auf Parallelflaechen ist nabla = n d_nu + D^(nu) mit D^(nu) = D - nu D_S + O(nu^2),
// D = sum_i e_i d_i (tangentialer Dirac-Operator), D_S = sum_i kappa_i e_i d_i (mit dem Formoperator S = grad_G n
// gewichtet). Aus (nabla - ik) F = 0 folgt d_nu F = B(nu) F, B = n (ik - D^(nu)), also fuer einen Schritt s ab nu0
//   U = 1 + s B + (s nu0 + s^2/2) B' + s^2/2 B^2 + O(s^3),   B' = n D_S.
// J_eff ist das Produkt aus Schritten und Transmissionsabbildungen von der Aussenspur (Medium 2 bei nu = 0) durch alle
// Schichten bis in Medium 1 und zurueck zu nu = 0. Dirac1 ist die erste Ordnung (dieselbe Naeherung wie die Sprungform,
// anders diskretisiert), Dirac2 die zweite Ordnung einschliesslich Kruemmung (Fehler O(d^3)). Zweite Flaechenableitungen
// sind hintereinandergeschaltete Kleinste-Quadrate-Gradienten (im L^2-Mittel konsistent).
// Chirale Medien (Kern und Schichten, v0.17): Mit den zentralen Helizitaetsprojektoren P_pm = (1 +- iI)/2 und
// (nabla - i k_pm) F_pm = 0 wird ik durch den zentralen Multivektor iK, K = k_+ P_+ + k_- P_-, ersetzt; alle Umformungen
// bleiben gueltig (K vertauscht mit n und D), B^2 F = -K^2 F - DDF - 2H n (iKF - DF). Das Aussenmedium bleibt achiral.
// Normalen: Die Facettennormale eines ebenen Dreiecks weicht um O(h) von der glatten Normalen am Schwerpunkt ab; abgeleitet
// ueber den Abstand h wird daraus ein Fehler O(1) (Formoperator auf der Kugel: 15 % von 1/R, nicht konvergent). Alle
// Schichtkorrekturen verwenden deshalb die glatte Normale nsm (Mittel der winkelgewichteten Knotennormalen). Damit ohne
// Schicht exakt T_1 bleibt, ist J_eff = J(n) + [Kette(nsm) - J(nsm)] (Dirac-Form) bzw. J(n) + L_d(nsm) (Sprungform).
#include <memory>
#include <vector>
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/multibody_operator.hpp"
#include "cbem/problems/layered_problem.hpp"

namespace cbem {

// Finite-Volumen-Flaechenoperatoren auf einem geschlossenen, konformen Dreiecksnetz (Werte je Dreieck)
struct SurfaceFV {
    struct Edge { std::size_t nb; real len; Vec3 conormal; Vec3 nb_conormal; Vec3 lsq; };   // Nachbar, Laenge, Konormalen
    std::vector<std::array<Edge, 3>> edges;                                    // (nach aussen), LSQ-Gewicht fuer grad
    std::vector<Vec3> nsm;                                                     // glatte Normale am Schwerpunkt (Mittel der Knotennormalen)
    std::vector<std::array<real, 9>> shape;                                    // Formoperator S = grad_G nsm (symmetrisch, tangential)
    // Quadratische Anpassung ueber alle Dreiecke mit gemeinsamem Knoten (Tangentialkoordinaten zur glatten Normalen):
    // phi_j - phi_t = g.r + 1/2 r^T H r; Gewichte fuer den Gradienten (O(h^2)) und fuer Delta_G = tr H (O(h)).
    struct Ring { std::vector<std::size_t> nb; std::vector<Vec3> gw; std::vector<real> lw; };
    std::vector<Ring> ring;
    std::vector<real> meancurv;                                                // H = tr S / 2
    explicit SurfaceFV(const TriangleMesh& m);
    const TriangleMesh& mesh;
    // grad_G phi (tangential im Dreieck), div_G v (v tangential je Dreieck)
    void grad(const std::vector<cplx>& phi, std::vector<std::array<cplx, 3>>& g) const;
    void div(const std::vector<std::array<cplx, 3>>& v, std::vector<cplx>& dv) const;
    // tangentialer Dirac-Operator D F = sum_c grad_G(F_c) * e_c je Dreieck (Werte); optional D_S F (mit S gewichtet)
    void dirac(const std::vector<Multivector>& F, std::vector<Multivector>& DF, std::vector<Multivector>* DSF = nullptr) const;
    // mit der quadratischen Anpassung: D F, D_S F und D D F = sum_c [Delta_G F_c - (S grad F_c) n] e_c (punktweise konsistent)
    void dirac_fit(const std::vector<Multivector>& F, std::vector<Multivector>& DF, std::vector<Multivector>& DSF,
                   std::vector<Multivector>& DDF) const;
    void laplace_fit(const std::vector<cplx>& phi, std::vector<cplx>& lap) const;
};

// Sprungform 1. Ordnung (v0.14), Dirac-Form 1. bzw. 2. Ordnung mit zusammengesetzten Gradienten (v0.15), Dirac-Form
// 2. Ordnung mit quadratischer Anpassung und geschlossener Form B^2 F = -k^2 F - DDF - 2H n (ikF - DF) (v0.16)
enum class ThinLayerModel { Jump1, Dirac1, Dirac2, Dirac2Fit };

class ThinLayerTransmissionOperator {
public:
    ThinLayerTransmissionOperator(const TriangleMesh& m, const BoundaryOperator& E_inner, const BoundaryOperator& E_outer,
                                  const Medium& inner, const Medium& outer, const std::vector<Coating>& coatings,
                                  real inner_fraction, real omega, ThinLayerModel model = ThinLayerModel::Jump1);
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
    ThinLayerModel model_;
    // Dirac-Form: Medien von aussen nach innen (2, Schichten, 1), Grenzen nu_l, Transmissionsabbildungen je Dreieck
    std::vector<Medium> stack_; std::vector<real> bound_; std::vector<std::vector<Mat8>> Jstep_; std::vector<Mat8> Jsm_;
    void propagate(std::vector<Multivector>& F, const Medium& m, real nu0, real nu1) const;
    void apply_dirac(const std::vector<cplx>& x, std::vector<cplx>& y) const;
};

// Koerper mit duennen Schichten: Flaeche (Referenzflaeche), Kernmedium, Schichten von innen nach aussen, Anteil der
// Schichten innerhalb der Referenzflaeche (0: Referenz = Kernoberflaeche, empfohlen fuer Metallkerne)
struct ThinBody { TriangleMesh surface; Medium core; std::vector<Coating> coatings; real inner_fraction = 0.0; };

// Streuproblem: ein oder mehrere Koerper mit duennen Schichten (je eine Flaeche). Aussen wirkt E_2 auf der Vereinigung
// aller Flaechen, innen ein blockdiagonaler Operator (ein Cauchy-Operator je Kern), J_eff je Koerper.
class ThinLayerScatteringProblem {
public:
    ThinLayerScatteringProblem(const TriangleMesh& surface, const Medium& core, const std::vector<Coating>& coatings,
                               real omega, Medium outer = {}, real inner_fraction = 0.0, HMatrixParams hp = {}, EntryParams ep = {},
                               ThinLayerModel model = ThinLayerModel::Dirac2Fit);
    ThinLayerScatteringProblem(const std::vector<ThinBody>& bodies, real omega, Medium outer = {}, HMatrixParams hp = {},
                               EntryParams ep = {}, ThinLayerModel model = ThinLayerModel::Dirac2Fit);
    LayeredResult solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o = {}) const;
    const TriangleMesh& mesh() const { return all_.all; }
    std::size_t body_begin(std::size_t b) const { return all_.body_begin[b]; }
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;         // T_eff
    void precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    double hmatrix_bytes() const;
private:
    void build(const std::vector<ThinBody>& bodies, HMatrixParams hp, EntryParams ep, ThinLayerModel model);
    Medium outer_;
    real omega_;
    MultiBodyMesh all_;
    std::vector<std::unique_ptr<TriangleMesh>> body_mesh_;
    std::vector<std::unique_ptr<KernelEntries>> ents_;
    std::vector<std::unique_ptr<KernelHMatrix>> hms_;
    std::vector<std::unique_ptr<CauchyOperator>> cops_;
    std::vector<std::unique_ptr<ChiralCauchyOperator>> chops_;               // chirale Kerne
    std::unique_ptr<BlockDiagonalOperator> E1_;
    const CauchyOperator* E2_ = nullptr;
    std::vector<std::unique_ptr<ThinLayerTransmissionOperator>> maps_;          // J_eff je Koerper (lokale Nummerierung)
};

}  // namespace cbem
