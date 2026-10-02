#pragma once
// Eintragsauswertung des Cauchy-Operators auf quadratischen (gekruemmten) Elementen mit unstetig linearen Dichten
// (v0.47, Stufe 2b). Die Normale variiert im Element und steht deshalb im Integral:
//   K^{ab}(i, j) = int_{tau_i} int_{tau_j} psi_a(x) Phi_k(x - y) n(y) psi_b(y) dS_y dS_x   (Multivektor),
// mit Phi_k(z) n = s n + v (z . n) + v (z ^ n): sieben Komponenten in der Blade-Reihenfolge
//   [1, e1, e2, e3, e12, e13, e23]   (Skalar v z.n; Vektor s n; Bivektor v z^n).
// Fuer ebene Elemente: Komponenten = (K_vec . n, K_0 n, K_vec ^ n) der vier Komponenten von KernelEntries.
//
// Quadratur im Parameterraum (dS = |X_u x X_v| du dv):
//  - Fernpaare: Dunavant 7 x 7;
//  - benachbarte Paare: Sauter-Schwab in Parameterkoordinaten. Selbstterm: Phi_0(x - y) n(y) = K_a + K_s mit dem
//    antisymmetrischen K_a = Phi_0 (n(x) + n(y))/2 (Hauptwert, mit [psi_a(x) psi_b(y) - psi_a(y) psi_b(x)]/2 absolut
//    integrierbar) und dem schwach singulaeren K_s = Phi_0 (n(y) - n(x))/2;
//  - nahe getrennte Paare: doppelt adaptiv (aeusseres Teilstueck klein gegen den Abstand zum inneren Element, inneres
//    Teilstueck klein gegen den Abstand zum aeusseren Punkt), je Dunavant 7. Analytische Innenintegrale gibt es auf
//    gekruemmten Elementen nicht; die Kriterien sind an der ebenen Gegenprobe gegen LinearKernelEntries eingestellt.
// Basis je Element orthonormal bezueglich der gekruemmten Flaeche: psi = S lambda, S = L^{-1} (curved_psi_matrices).
#include <array>
#include <vector>

#include "cbem/assembly/kernel_entries.hpp"
#include "cbem/geometry/quadratic_mesh.hpp"

namespace cbem {

constexpr int kCurvedComps = 7;
using CurvedComp = std::array<cplx, kCurvedComps>;
using CurvedBlock = std::array<CurvedComp, 9>;   // [a * 3 + b]
// Blades der sieben Komponenten (Bitmasken)
constexpr int kCurvedBlade[kCurvedComps] = {0, 1, 2, 4, 3, 5, 6};

struct CurvedNearParams {
    // doppelt adaptives Verfahren (subtract = false): relativer Fehler etwa 2e-6 (Selbstkonvergenz: 0,5/0,25 -> 2e-5,
    // 0,3/0,15 -> 2e-6, 0,2/0,1 -> 2e-7);
    // zum Vergleich: die analytische Nahquadratur der ebenen Elemente liegt bei 3,5e-6 (Rest mit 7 Gauss-Punkten)
    real outer_ratio = 0.3;    // aeusseres Teilstueck: Umkreisradius < outer_ratio * Abstand zum inneren Element
    real inner_ratio = 0.15;   // inneres Teilstueck: Umkreisradius < inner_ratio * Abstand zum aeusseren Punkt
    int outer_depth = 12, inner_depth = 16;
    // Singularitaetssubtraktion (v0.49): vom Innenintegral wird der singulaere Kern ueber dem Tangentialdreieck am Fusspunkt
    // u* des aeusseren Punktes analytisch abgezogen (triangle_integrals_linear, Normale n(u*), Jacobi-Determinante J(u*));
    // der Rest verhaelt sich bei u* wie 1/r statt 1/r^3 und wird mit dem groeberen Kriterium correction_ratio integriert.
    // Begrenzend war der Grad der Blattregel, nicht die Singularitaet (v0.52): Der Rest ist auf der Skala des Abstands glatt,
    // eine Regel hohen Grades mit lockerem Kriterium ist genauer und billiger als Dunavant 7 mit feiner Unterteilung.
    // Voreinstellung (v0.52) Gauss 5 x 5, Kriterien 1,5/1,5: Fehler 0,6-2,6e-6 (v0.49 mit Dunavant 7 und 0,3/0,3: 5-10e-6),
    // 3-mal schneller. Schnell: Gauss 4 x 4, 1,0/1,0 (Fehler 2,5-4e-5 wie v0.49 mit 0,5/0,5, gleich schnell).
    // false: doppelt adaptiv wie in v0.47 (Referenz, outer_ratio/inner_ratio).
    bool subtract = true;
    real subtract_outer_ratio = 1.5;
    real correction_ratio = 1.5;
    // Regeln an den Blaettern der aeusseren Integration und der Korrektur: 0 = Dunavant 7 (Grad 5), n >= 2 =
    // QuadRule::conical(n) (n x n Punkte, Grad 2n - 2)
    int outer_rule = 5, correction_rule = 5;
};

class CurvedKernelEntries {
public:
    static constexpr int block_size = 3;
    CurvedKernelEntries(const QuadraticMesh& mesh, cplx k, EntryParams prm = {}, CurvedNearParams np = {});
    const QuadraticMesh& mesh() const { return m_; }
    cplx wavenumber() const { return k_; }
    std::size_t elements() const { return m_.size(); }
    std::size_t size() const { return 3 * m_.size(); }

    // Elementpaare in der Basis lambda (baryzentrische Parameter)
    CurvedBlock lambda_far(std::size_t i, std::size_t j) const;
    CurvedBlock lambda_exact(std::size_t i, std::size_t j) const;
    CurvedBlock lambda_sauter_schwab(std::size_t i, std::size_t j, Adjacency a) const;
    CurvedBlock lambda_near(std::size_t i, std::size_t j) const;
    // orthonormierte Basis psi (Nahpaare gecacht)
    CurvedBlock block(std::size_t i, std::size_t j) const;
    CurvedBlock block_far(std::size_t i, std::size_t j) const;

    bool is_near(std::size_t i, std::size_t j) const;
    Adjacency adjacency(std::size_t i, std::size_t j) const;
    const std::array<real, 9>& S(std::size_t t) const { return S_[t]; }
    std::vector<Vec3> index_points() const;
    std::size_t near_pairs() const { return n_near_; }
    double near_seconds() const { return t_near_; }

private:
    const QuadraticMesh& m_;
    cplx k_;
    EntryParams prm_;
    CurvedNearParams np_;
    QuadRule r7_, ro_, rc_;   // Dunavant 7; Blattregeln der aeusseren Integration und der Korrektur
    CurvedQuadrature q7_;
    PairRule ss_[4];
    std::vector<std::array<real, 9>> S_;
    std::vector<std::array<real, 3>> psiw_;   // w psi_a an den Punkten von q7_ je Element (block_far)
    // Geometrie je Element als Polynom in u = lambda_1, v = lambda_2 (v0.51): X = A + B u + C v + D u^2 + E u v + F v^2,
    // inline ausgewertet (QuadraticMesh::X/jacobian sammeln die Knoten je Aufruf ueber die Konnektivitaet)
    struct Poly { Vec3 A, B, C, D, E, F; };
    std::vector<Poly> poly_;
    Vec3 gX(std::size_t t, const std::array<real, 3>& l) const {
        const Poly& g = poly_[t]; const real u = l[1], v = l[2];
        return g.A + (g.B + g.D * u + g.E * v) * u + (g.C + g.F * v) * v;
    }
    void gFrame(std::size_t t, const std::array<real, 3>& l, Vec3& Xu, Vec3& Xv) const {
        const Poly& g = poly_[t]; const real u = l[1], v = l[2];
        Xu = g.B + g.D * (2 * u) + g.E * v; Xv = g.C + g.E * u + g.F * (2 * v);
    }
    real gJ(std::size_t t, const std::array<real, 3>& l, Vec3* n = nullptr) const {
        Vec3 Xu, Xv; gFrame(t, l, Xu, Xv);
        const Vec3 c = cross(Xu, Xv); const real J = norm(c);
        if (n) *n = c / J;
        return J;
    }
    std::vector<std::vector<std::pair<std::size_t, CurvedBlock>>> cache_;
    bool cached_ = false;
    std::size_t n_near_ = 0;
    double t_near_ = 0;
    CurvedBlock to_psi(std::size_t i, std::size_t j, const CurvedBlock& L) const;
    void inner(const Vec3& x, std::size_t j, const std::array<std::array<real, 3>, 3>& tri, real aref, int depth,
               std::array<CurvedComp, 3>& acc) const;
    struct Tangent { Vec3 X0, Xu, Xv, n; real J; std::array<real, 3> lam; };   // affine Taylor-Abbildung am Fusspunkt
    Tangent tangent_at(const Vec3& x, std::size_t j) const;
    void inner_subtracted(const Vec3& x, std::size_t j, std::array<CurvedComp, 3>& acc) const;
    void correction(const Vec3& x, std::size_t j, const Tangent& T, const std::array<std::array<real, 3>, 3>& tri, real aref, int depth,
                    std::array<CurvedComp, 3>& acc) const;
    void correction_point(const Vec3& x, std::size_t j, const Tangent& T, const std::array<real, 3>& l, real w,
                          std::array<CurvedComp, 3>& acc) const;
    void outer(std::size_t i, std::size_t j, const std::array<std::array<real, 3>, 3>& tri, real aref, int depth, CurvedBlock& K) const;
    real distance_to_element(const Vec3& x, std::size_t j) const;
    void build_near_cache();
};

}  // namespace cbem
