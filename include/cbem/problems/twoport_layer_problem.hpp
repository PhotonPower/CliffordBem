#pragma once
// Beschichtete Grenzflaeche als Zweitor (S-Matrix-Formulierung, Prototyp v0.19).
//
// Kernoberflaeche Gamma_b und Aussenflaeche Gamma_a der Schicht (Medium c, Dicke d) mit gleicher Topologie (Parallelflaeche;
// Dreieck i von Gamma_a liegt ueber Dreieck i von Gamma_b). Unbekannt sind die Aussenspur h_a (Medium 2) auf Gamma_a und die
// Kernspur v_b (Medium 1) auf Gamma_b. Kein Feld wird durch die Schicht fortgesetzt, es gibt keinen Cauchy-Operator des
// Schichtgebiets.
//
//   Zeile I  (auf Gamma_a):  1/2 (1 + E_2) h_a + iota 1/2 (1 - E_1) v_b = h_inc
//   Zeile II (auf Gamma_b):  u_a - u_b - g(s) B_m (u_a + u_b) = 0,   u_a = J_{c<-2} h_a,  u_b = J_{c<-1} v_b.
//
// Zeile I paart die Bedingung des Aussenraums mit der des Kerns (nicht jede Flaeche mit ihren Nachbargebieten: diese Paarung
// degeneriert fuer d -> 0). Zeile II ist die Summe der beiden S-Matrix-Bedingungen Pi_+ u_b = P Pi_+ u_a (nach unten
// abklingend) und Pi_- u_a = P Pi_- u_b (nach oben abklingend) mit P = exp(-gamma d), Pi_pm = (1 +- B/gamma)/2, gamma^2 = s:
// Die Projektoren heben sich heraus und es
// bleibt (u_a - u_b) = g(s) B (u_a + u_b) mit g(s) = tanh(d sqrt(s)/2)/sqrt(s) -- eine durch tanh stabilisierte Trapezregel
// fuer d_nu F = B F. g ist in s ganz bis auf Pole bei s = -((2m+1) pi/d)^2 (kein Verzweigungsschnitt) und exakt
//   g(s) = sum_m (4/d) / (s + sigma_m),  sigma_m = ((2m+1) pi / d)^2,
// hier mit s = -K_c^2 - Delta_G (Laplace-Beltrami je Komponente, quadratische Anpassung), M Pole mit duennbesetzten Loesungen
// (-Delta_G + sigma_m - k_c^2)^{-1} und dem Rest linear in s. d << h: g -> d/2 - d^3 s/24 (lokal); d >> h: g ~ 1/sqrt(s)
// (die Flaechen entkoppeln fuer feine Strukturen). B_m = n (iK_c - D + d/2 D_S) ist B in der Schichtmitte (Kruemmung bis O(d^2)).
// d = 0: Zeile II erzwingt v_b = J_{1<-2} h_a exakt, Zeile I wird T_1. Chirale Schicht und chiraler Kern ueber K = k+P+ + k-P-.
#include <memory>
#include <vector>
#include "cbem/problems/thin_layer_problem.hpp"

namespace cbem {

struct TwoPortOptions {
    int poles = -1;             // Zahl der Pole von g; -1: automatisch (sigma_M >= 10 s_max)
    real inner_tol = 1e-10;     // GMRES-Toleranz der Resolventen
    bool curvature = true;      // B in der Schichtmitte (d/2 D_S)
};

class TwoPortLayerProblem {
public:
    // core_surface = Gamma_b, outer_surface = Gamma_a (gleiche Topologie), d = Schichtdicke
    TwoPortLayerProblem(const TriangleMesh& core_surface, const TriangleMesh& outer_surface, const Medium& core, const Medium& layer,
                        real d, real omega, Medium outer = {}, HMatrixParams hp = {}, EntryParams ep = {}, TwoPortOptions opt = {});
    LayeredResult solve_plane_wave(const Vec3& dir, const CVec3& p, const SolveOptions& o = {}) const;
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    void precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    std::size_t size() const { return 16 * N_; }
    int poles() const { return M_; }
    double mean_inner_iterations() const { return inner_calls_ ? double(inner_its_) / inner_calls_ : 0.0; }
    const TriangleMesh& outer_mesh() const { return ma_; }
private:
    void apply_g(const std::vector<Multivector>& z, std::vector<Multivector>& gz) const;   // g(s) z (Werte)
    void resolvent(const std::vector<Multivector>& z, cplx shift, std::vector<Multivector>& x) const;  // (-Delta + shift)^{-1} z
    TriangleMesh mb_, ma_;
    std::size_t N_;
    Medium core_, layer_, outer_;
    real d_, omega_;
    TwoPortOptions opt_;
    std::unique_ptr<SurfaceFV> fv_;
    std::unique_ptr<KernelEntries> K2_, K1p_, K1m_;
    std::unique_ptr<KernelHMatrix> H2_, H1p_, H1m_;
    std::unique_ptr<CauchyOperator> E2_, E1p_, E1m_;
    std::unique_ptr<ChiralCauchyOperator> E1ch_;
    const BoundaryOperator* E1_ = nullptr;
    std::vector<Mat8> Jc2_, Jc1_, Jc2s_, Jc1s_;     // Facetten- bzw. glatte Normalen
    std::vector<real> rho_;                          // sqrt(|tau_a| / |tau_b|)
    std::vector<std::array<cplx, 256>> Pinv_;        // 16x16-Vorkonditionierer je Dreieck
    int M_ = 0; std::vector<real> sigma_; real tail0_ = 0, tail1_ = 0;
    mutable long inner_its_ = 0, inner_calls_ = 0;
};

}  // namespace cbem
