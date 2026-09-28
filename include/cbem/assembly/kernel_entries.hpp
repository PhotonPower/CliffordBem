#pragma once
// Eintragsauswertung fuer den Cauchy-Randoperator E_k auf stueckweise konstanten Dichten:
// Kernkomponenten K(i,j) = int_{tau_i} int_{tau_j} Phi_k(x - y) dS_y dS_x  (Skalar + Vektor, 4 Werte).
// Fernpaare: Gauss 7x7. Nahpaare: Phi_0 und 1/r analytisch (inneres Integral), Rest mit Gauss;
// Selbstterm des Phi_0-Anteils exakt 0 (Antisymmetrie). Diese Klasse ist die Schnittstelle zur
// H-Matrix (AP 3) und spaeter zu Loesern (AP 2).
#include <vector>
#include "cbem/geometry/mesh.hpp"
#include "cbem/geometry/quadrature.hpp"
#include "cbem/geometry/sauter_schwab.hpp"

namespace cbem {

struct EntryParams {
    real near_factor = 2.5;   // Nahpaar, wenn Schwerpunktabstand < near_factor * max(h_i, h_j)
    int near_subdivision = 4; // aeussere Regel fuer Nahpaare: sub^2 * 7 Punkte
    bool sauter_schwab = true; // Paare mit gemeinsamer Ecke/Kante/Flaeche mit Sauter-Schwab (konforme Netze)
    int ss_order = 5;          // Gauss-Punkte je Richtung der Sauter-Schwab-Regeln
    bool cache_near = true;    // alle Nahpaare einmal vorab berechnen und speichern (Symmetrie K_ji = (s, -v))
    bool adaptive_outer = true;// aeussere Regel adaptiv: Halbierung der laengsten Kante, bis jedes Teildreieck
                               // klein gegen seinen Abstand zum inneren Dreieck ist (sonst feste sub^2*7-Regel)
    real adapt_ratio = 0.5;    // Kriterium: Umkreisradius des Teilstuecks < adapt_ratio * Abstand
    int adapt_depth = 14;      // maximale Halbierungstiefe
    bool adapt_to_boundary = false; // liegt das Teilstueck ganz auf einer Seite der Ebene des inneren Dreiecks, zaehlt der
                               // Abstand zu dessen Randkanten statt zum Dreieck: das analytische Innenintegral ist dort
                               // glatt bis auf den Rand (duenne Schichten: parallele Flaechen im Abstand d << h)
    int sa_order = 14;         // halbanalytische Regel fuer benachbarte Paare: Gauss-Punkte je Richtung
    real ss_max_aspect = 1.6;  // Sauter-Schwab nur, wenn beide Dreiecke Seitenverhaeltnis <= ss_max_aspect;
                               // sonst halbanalytisch (Innenintegral exakt, Aussenregel zur Singularitaet gradiert)
};

class KernelEntries {
public:
    KernelEntries(const TriangleMesh& mesh, cplx k, EntryParams prm = {});
    const TriangleMesh& mesh() const { return m_; }
    cplx wavenumber() const { return k_; }
    bool is_near(std::size_t i, std::size_t j) const;
    KernelComp far(std::size_t i, std::size_t j) const;      // nur fuer getrennte Paare genau
    KernelComp near(std::size_t i, std::size_t j) const;     // Nahfeldbehandlung (analytisches Innenintegral)
    KernelComp sauter_schwab(std::size_t i, std::size_t j, Adjacency a) const;
    KernelComp semi_analytic(std::size_t i, std::size_t j, Adjacency a) const;   // benachbarte Paare, beliebige Form
    real aspect(std::size_t t) const { return m_.hmax[t] * m_.hmax[t] / (2 * m_.area[t]); }
    Adjacency adjacency(std::size_t i, std::size_t j) const;
    KernelComp exact(std::size_t i, std::size_t j) const;
    std::size_t near_pairs() const { return n_near_; }
    double near_seconds() const { return t_near_; }
private:
    const TriangleMesh& m_;
    cplx k_;
    EntryParams prm_;
    MeshQuadrature q7_, q2_, qn_;
    std::vector<real> rho_;    // Umkreisradius um den Schwerpunkt
    PairRule ss_[4];
    KernelComp exact_uncached(std::size_t i, std::size_t j) const;
    void outer_point(const Vec3& x, real w, std::size_t inner, bool skip_phi0, KernelComp& K) const;
    void near_adaptive(const std::array<Vec3, 3>& outer, std::size_t inner, int depth, bool swap, bool self, KernelComp& K) const;
    void build_near_cache();
    std::vector<std::vector<std::pair<std::size_t, KernelComp>>> cache_;   // je Zeile nach j sortiert
    bool cached_ = false; std::size_t n_near_ = 0; double t_near_ = 0;
};

}  // namespace cbem
