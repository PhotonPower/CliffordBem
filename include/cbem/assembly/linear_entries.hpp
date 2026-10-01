#pragma once
// Eintragsauswertung des Cauchy-Operators E_k fuer unstetig lineare Dichten (v0.45, Stufe 2a der gekruemmten Elemente):
//   K^{ab}_c(i, j) = int_{tau_i} int_{tau_j} lambda_a(x) Phi_k(x - y)_c lambda_b(y) dS_y dS_x,   a, b = 0..2, c = 0..3,
// mit den baryzentrischen Koordinaten lambda_a zur Ecke T[i][a] (Skalar + Vektor wie KernelEntries). Gleiche Quadratur wie
// KernelEntries: Fernpaare Gauss 7 x 7; nahe getrennte Paare mit adaptiver aeusserer Regel und analytischem Innenintegral
// (triangle_integrals_linear); benachbarte Paare mit Sauter-Schwab bzw. halbanalytisch. Im Selbstterm ist der singulaere
// Anteil Phi_0 = z/(4 pi r^3) fuer lineare Dichten nicht null (fuer konstante verschwindet er aus Antisymmetrie); er wird
// antisymmetrisiert integriert: int int Phi_0(x - y) [lambda_a(x) lambda_b(y) - lambda_a(y) lambda_b(x)] / 2, absolut
// integrierbar (Sauter-Schwab; bei gestreckten Elementen halbanalytisch mit dem Hauptwert des Innenintegrals). Summe ueber a, b
// = KernelEntries (sum lambda = 1): auf Maschinengenauigkeit, nur gestreckte Selbstterme auf den Quadraturfehler der
// halbanalytischen Regel (1e-5 bei sa_order = 14, faellt mit sa_order).
//
// Fuer die H-Matrix und die Operatoren: je Element orthonormierte Basis psi_a = sum_k S_ak lambda_k (S = G^{-1/2}, G die
// Gram-Matrix der lambda), Indizes I = 3 t + a. Damit ist die Massenmatrix die Identitaet, und T_1 behaelt seine Form.
#include <array>
#include <cmath>
#include <vector>

#include "cbem/assembly/kernel_entries.hpp"

namespace cbem {

using LinearBlock = std::array<KernelComp, 9>;   // [a * 3 + b]: lambda_a (bzw. psi_a) auf tau_i, lambda_b (psi_b) auf tau_j

// Orthonormierte Basis eines Dreiecks der Flaeche A: psi_a = sum_k S[a * 3 + k] lambda_k, S = G^{-1/2} mit der Gram-Matrix
// G = (A/12)(1 + 1 1^T) (Eigenwerte A/3 auf der Konstanten, A/12 senkrecht dazu)
inline std::array<real, 9> psi_matrix(real A) {
    const real s1 = std::sqrt(12.0 / A), s3 = std::sqrt(3.0 / A);
    std::array<real, 9> S{};
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) S[a * 3 + b] = (a == b ? s1 : 0.0) + (s3 - s1) / 3.0;
    return S;
}

class LinearKernelEntries {
public:
    static constexpr int block_size = 3;
    LinearKernelEntries(const TriangleMesh& mesh, cplx k, EntryParams prm = {});
    const TriangleMesh& mesh() const { return m_; }
    cplx wavenumber() const { return k_; }
    std::size_t elements() const { return m_.size(); }
    std::size_t size() const { return 3 * m_.size(); }

    // Elementpaare in der baryzentrischen Basis lambda
    LinearBlock lambda_far(std::size_t i, std::size_t j) const;
    LinearBlock lambda_exact(std::size_t i, std::size_t j) const;   // ungecacht: Nahfeld wo noetig
    LinearBlock lambda_sauter_schwab(std::size_t i, std::size_t j, Adjacency a) const;
    LinearBlock lambda_semi_analytic(std::size_t i, std::size_t j, Adjacency a) const;
    LinearBlock lambda_near(std::size_t i, std::size_t j) const;

    // Elementpaare in der orthonormierten Basis psi (gecacht fuer Nahpaare)
    LinearBlock block(std::size_t i, std::size_t j) const;           // genau
    LinearBlock block_far(std::size_t i, std::size_t j) const;       // Gauss-Fernfeldregel
    // einzelne Eintraege der H-Matrix-Schnittstelle (Indizes I = 3 t + a), psi-Basis
    KernelComp exact(std::size_t I, std::size_t J) const { return block(I / 3, J / 3)[(I % 3) * 3 + J % 3]; }
    KernelComp far(std::size_t I, std::size_t J) const { return block_far(I / 3, J / 3)[(I % 3) * 3 + J % 3]; }

    bool is_near(std::size_t i, std::size_t j) const;
    Adjacency adjacency(std::size_t i, std::size_t j) const;
    real aspect(std::size_t t) const { return m_.hmax[t] * m_.hmax[t] / (2 * m_.area[t]); }
    const std::array<real, 9>& S(std::size_t t) const { return S_[t]; }   // psi_a = sum_k S[a * 3 + k] lambda_k
    real lambda(std::size_t t, int k, const Vec3& x) const { return alpha_[3 * t + k] + dot(beta_[3 * t + k], x); }
    // Punkte und Groessen je Index fuer den Clusterbaum: (Schwerpunkt + Ecke a)/2, Elementgroesse h_max
    std::vector<Vec3> index_points() const;
    std::vector<real> index_sizes() const;
    std::size_t near_pairs() const { return n_near_; }
    double near_seconds() const { return t_near_; }

private:
    const TriangleMesh& m_;
    cplx k_;
    EntryParams prm_;
    QuadRule r7_, rn_;
    MeshQuadrature q7_, qn_;
    PairRule ss_[4];
    std::vector<real> alpha_;   // lambda_k(x) = alpha + beta . x, Index 3 t + k
    std::vector<Vec3> beta_;
    std::vector<std::array<real, 9>> S_;
    std::vector<std::vector<std::pair<std::size_t, LinearBlock>>> cache_;   // psi-Basis, je Zeile nach j sortiert
    bool cached_ = false;
    std::size_t n_near_ = 0;
    double t_near_ = 0;
    LinearBlock to_psi(std::size_t i, std::size_t j, const LinearBlock& L) const;
    void outer_point(const Vec3& x, real w, const real lo[3], std::size_t inner, bool swap, LinearBlock& K) const;
    void near_adaptive(const std::array<Vec3, 3>& outer, std::size_t outer_t, std::size_t inner, int depth, bool swap, LinearBlock& K) const;
    void build_near_cache();
};

}  // namespace cbem
