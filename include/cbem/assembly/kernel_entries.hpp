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
    Adjacency adjacency(std::size_t i, std::size_t j) const;
    KernelComp exact(std::size_t i, std::size_t j) const;
private:
    const TriangleMesh& m_;
    cplx k_;
    EntryParams prm_;
    MeshQuadrature q7_, qn_;
    PairRule ss_[4];
};

}  // namespace cbem
