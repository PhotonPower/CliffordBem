#pragma once
// Beschichtete Grenzflaechen als Zweitor (S-Matrix-Formulierung; v0.19 eine Schicht, v0.20 Mehrfachschichten).
//
// Flaechen Gamma_0 (Kernoberflaeche), Gamma_1, ..., Gamma_L (Aussenflaeche) gleicher Topologie (Parallelflaechen; Dreieck i
// liegt ueber Dreieck i von Gamma_0). Medien m_0 = Kern, m_1 ... m_L = Schichten von innen nach aussen, m_{L+1} = aussen.
// Unbekannt: die Aussenspur X_L (Medium m_{L+1}) auf Gamma_L, die Kernspur v (Medium m_0) auf Gamma_0 und fuer 0 < l < L
// die Spur X_l im Medium m_{l+1} auf Gamma_l. E_2 wirkt nur auf Gamma_L, E_1 nur auf Gamma_0 -- zwei H-Matrizen, unabhaengig
// von der Zahl der Schichten; kein Feld wird durch eine Schicht fortgesetzt, es gibt keinen Cauchy-Operator eines Schichtgebiets.
//
//   Zeile I (auf Gamma_L):  1/2 (1 + E_2) X_L + iota 1/2 (1 - E_1) v = h_inc
//   Zeile l = 1..L (auf Gamma_0 nummeriert):  u_top - u_bot - g_l(s) B_l (u_top + u_bot) = 0,
//     u_top = J_{m_l <- m_{l+1}} X_l,  u_bot = J_{m_1 <- m_0} v (l = 1) bzw. X_{l-1} (l > 1).
//
// Zeile I paart die Bedingung des Aussenraums mit der des Kerns (jede Flaeche an ihre Nachbargebiete zu koppeln degeneriert fuer
// d -> 0). Zeile l ist die Summe der beiden S-Matrix-Bedingungen der Schicht l, Pi_+ u_bot = P Pi_+ u_top (nach unten
// abklingend) und Pi_- u_top = P Pi_- u_bot (nach oben abklingend), P = exp(-gamma d_l), Pi_pm = (1 +- B/gamma)/2,
// gamma^2 = s: Die Projektoren heben sich heraus, es bleibt eine durch tanh stabilisierte Trapezregel fuer d_nu F = B F mit
//   g_l(s) = tanh(d_l sqrt(s)/2)/sqrt(s) = sum_m (4/d_l) / (s + sigma_m),  sigma_m = ((2m+1) pi / d_l)^2,
// ganz in s bis auf weit entfernte Pole (kein Verzweigungsschnitt); s = -K_l^2 - Delta_G (Laplace-Beltrami auf Gamma_0 je
// Komponente), einige Pole mit duennbesetzten Resolventen, Rest linear in s. d << h: g -> d/2 - d^3 s/24 (lokal); d >> h: die
// Flaechen entkoppeln fuer feine Strukturen. Jede Schicht l wird auf ihrer eigenen unteren Flaeche Gamma_{l-1} ausgewertet
// (eigene Stencils, Formoperator, Laplace-Beltrami): B_l = n (iK_l - D + d_l/2 D_S) in ihrer Mitte (Kruemmung bis O(d_l^2)).
// Eine Entwicklung um die Kernoberflaeche (nu D_S mit nu = Abstand der Schichtmitte vom Kern) ist fuer weit aussen liegende
// Schichten zu grob (chirale Schicht 4,5 nm ueber einer 20-nm-Goldkugel: CD um 75 % falsch). Das Gleichungssystem ist blockbidiagonal und wird als Ganzes geloest: keine Produkte von
// Transfermatrizen, jede Zeile beschraenkt (Stabilitaet wie beim Redheffer-Sternprodukt). Alle d_l = 0: jede Zeile wird zur
// Stetigkeit, mit J_{a<-b} J_{b<-c} = J_{a<-c} folgt v = J_{m_0 <- m_{L+1}} X_L und Zeile I wird T_1. Chirale Medien ueber
// K = k+ P+ + k- P-, auch im Aussenmedium (v0.21). Fuer eine Schicht identisch mit v0.19.
#include <memory>
#include <vector>
#include "cbem/problems/thin_layer_problem.hpp"

namespace cbem {

struct TwoPortOptions {
    int poles = -1;             // Zahl der Pole je Schicht; -1: automatisch (sigma_M >= 10 s_max)
    real inner_tol = 1e-10;     // GMRES-Toleranz der Resolventen
    bool curvature = true;      // B in der Schichtmitte (nu D_S)
    // Unterteilung dicker Schichten: jede Teilschicht hoechstens split / kappa_max dick (kappa_max: groesste Hauptkruemmung
    // der Kernflaeche); Zwischenflaechen durch lineare Interpolation der Knoten (exakt fuer konzentrische Kugeln und
    // Parallelflaechen). Der Kruemmungsfehler einer Schicht waechst mit (d/a)^2; die Unterteilung kostet nur lokale
    // Unbekannte, keine weiteren H-Matrizen. 0: keine Unterteilung.
    real split = 0.1;
};

// Ein Koerper fuer das Zweitor: Flaechen Gamma_0 (Kern) ... Gamma_L (aussen) gleicher Topologie, Kernmedium, Schichten von
// innen nach aussen
struct TwoPortBody { std::vector<TriangleMesh> surfaces; Medium core; std::vector<Coating> layers; };

// Mehrere Koerper (v0.22): E_2 wirkt auf der Vereinigung aller Aussenflaechen, E_1 blockdiagonal auf den Kernen; jeder
// Koerper hat seinen eigenen Schichtstapel mit eigenen Zweitor-Zeilen. Anordnung der Unbekannten: [X_L aller Koerper |
// v aller Koerper | Zwischenspuren Koerper 0 | Koerper 1 | ...]; fuer einen Koerper identisch mit v0.20.
class TwoPortLayerProblem {
public:
    // eine Schicht: core_surface = Gamma_0, outer_surface = Gamma_1, d = Schichtdicke
    TwoPortLayerProblem(const TriangleMesh& core_surface, const TriangleMesh& outer_surface, const Medium& core, const Medium& layer,
                        real d, real omega, Medium outer = {}, HMatrixParams hp = {}, EntryParams ep = {}, TwoPortOptions opt = {});
    // Mehrfachschichten: surfaces = Gamma_0 ... Gamma_L (gleiche Topologie), layers von innen nach aussen (Dicke, Medium)
    TwoPortLayerProblem(const std::vector<TriangleMesh>& surfaces, const Medium& core, const std::vector<Coating>& layers,
                        real omega, Medium outer = {}, HMatrixParams hp = {}, EntryParams ep = {}, TwoPortOptions opt = {});
    // mehrere Koerper
    TwoPortLayerProblem(const std::vector<TwoPortBody>& bodies, real omega, Medium outer = {}, HMatrixParams hp = {},
                        EntryParams ep = {}, TwoPortOptions opt = {});
    // Parallelflaechen aus der Kernoberflaeche (offset_surface, kumulative Dicken; Pruefung auf Faltung/Durchdringung)
    static std::vector<TriangleMesh> layer_surfaces(const TriangleMesh& core_surface, const std::vector<Coating>& layers);
    LayeredResult solve_plane_wave(const Vec3& dir, const CVec3& p, const SolveOptions& o = {}) const;
    // beliebige einfallende Spur b auf outer_mesh() (z. B. project_dipole); r.h: Gesamtvektor, Aussenspur = erste 8 N (v0.26)
    LayeredResult solve_rhs(const std::vector<cplx>& b_outer, const SolveOptions& o = {}) const;
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    void precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    std::size_t size() const { return size_; }
    std::size_t bodies() const { return B_.size(); }
    std::size_t layers(std::size_t b = 0) const { return B_[b]->L; }            // nach der Unterteilung
    std::size_t input_layers(std::size_t b = 0) const { return B_[b]->L_in; }
    int poles(std::size_t l = 0, std::size_t b = 0) const { return B_[b]->lay[l].M; }
    double mean_inner_iterations() const { return inner_calls_ ? double(inner_its_) / inner_calls_ : 0.0; }
    const TriangleMesh& outer_mesh() const { return outer_.all; }             // Vereinigung der Aussenflaechen
private:
    struct Layer {
        Medium m; real d = 0, numid = 0;             // numid: Abstand der Mitte von der unteren Flaeche Gamma_{l-1}
        int M = 0; std::vector<real> sigma; real tail0 = 0, tail1 = 0;
        std::vector<Mat8> Jtop, Jtops;              // X_l (Medium m_{l+1}) -> Medium m_l, Facetten- bzw. glatte Normalen
    };
    struct Body {
        std::vector<TriangleMesh> S; Medium core; std::vector<Layer> lay;
        std::size_t N = 0, L = 0, L_in = 0, pre = 0;               // pre: erster Index in der Vereinigung
        std::size_t offXL = 0, offV = 0; std::vector<std::size_t> offInt;   // Positionen (Koeffizienten) im Vektor
        std::vector<std::unique_ptr<SurfaceFV>> fv;                // je Flaeche Gamma_0 ... Gamma_{L-1}
        std::vector<Mat8> Jb, Jbs;                                 // Kernspur v -> Medium m_1
        std::vector<std::vector<real>> sq;                         // sqrt|tau| je Flaeche
        std::vector<std::vector<cplx>> Pinv;                       // Vorkonditionierer je Dreieck (8(L+1))^2
        std::size_t pos(std::size_t j) const { return j == L ? offXL : j == 0 ? offV : offInt[j - 1]; }   // Position von X_j
        std::size_t row(std::size_t l) const { return l == 0 ? offXL : l == 1 ? offV : offInt[l - 2]; }    // Zeile I bzw. Schicht l
    };
    void build(std::vector<TwoPortBody> bodies, HMatrixParams hp, EntryParams ep);
    void setup_body(Body& Bd);
    void apply_g(const Body& Bd, std::size_t l, const std::vector<Multivector>& z, std::vector<Multivector>& gz) const;
    void resolvent(const SurfaceFV& fv, const std::vector<Multivector>& z, cplx shift, std::vector<Multivector>& x) const;
    std::vector<std::unique_ptr<Body>> B_;
    std::size_t Ntot_ = 0, size_ = 0;
    Medium outer_m_;
    real omega_;
    TwoPortOptions opt_;
    MultiBodyMesh outer_;
    std::vector<std::unique_ptr<KernelEntries>> ents_;
    std::vector<std::unique_ptr<KernelHMatrix>> hms_;
    std::vector<std::unique_ptr<CauchyOperator>> cops_;
    std::vector<std::unique_ptr<ChiralCauchyOperator>> chops_;
    std::unique_ptr<BlockDiagonalOperator> E1bd_;
    const BoundaryOperator* E1_ = nullptr;           // Kerne (blockdiagonal)
    const BoundaryOperator* E2op_ = nullptr;          // Aussenoperator auf der Vereinigung (chiral: P+ E_{k+} + P- E_{k-})
    mutable long inner_its_ = 0, inner_calls_ = 0;
};

}  // namespace cbem
