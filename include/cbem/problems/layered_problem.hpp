#pragma once
// Beschichtete Grenzflaechen und geschichtete Koerper (Kern-Schale, Mehrfachschichten, mehrere Koerper).
//
// Geometrie als Grenzflaechengraph: geschlossene, nach aussen orientierte Flaechen Gamma_s, jede trennt ein
// Innengebiet in(s) von einem Aussengebiet out(s); Gebiet 0 ist der Aussenraum. Unbekannte je Flaeche ist die
// Spur h_s des Feldes auf der Aussenseite (8 Komponenten je Dreieck); die Innenspur ist J_s h_s mit
// J_s = transmission_map(n, medium[in(s)], medium[out(s)]).
//
// Jedes Gebiet R mit Rand dR = Vereinigung seiner Flaechen hat den Cauchy-Operator E_R (Wellenzahl des Mediums,
// Flaechen mit ihren Aussennormalen). Mit sigma_s = +1, wenn R innen an Gamma_s liegt, sonst -1, und der Spur
// u_s = J_s h_s (R innen) bzw. h_s (R aussen) lautet die Cauchy-Bedingung des Gebiets (Normalen nach aussen aus R)
//   u - E_R(sigma u) = 0   (Aussenraum: fuer die Streuspur; das einfallende Feld liefert die rechte Seite).
// Zeile s des Systems ist die halbe Summe der Bedingungen beider Nachbargebiete:
//   1/2 (h_s - [E_out(sigma u)]_s) + 1/2 (J_s h_s - [E_in(sigma u)]_s) = h_inc,s   (h_inc nur an Flaechen zum Aussenraum).
// Fuer einen homogenen Koerper ist das genau T_1 = E_2^+ + E_1^- J; an jeder Flaeche hat das System dieselbe
// lokale Struktur (punktweise Vorkonditionierung 2 (1 + J_s)^{-1}). Eindeutigkeit ist fuer verschachtelte Gebiete
// nicht bewiesen (AP 1 behandelt eine Grenzflaeche); geprueft wird gegen die Aden-Kerker-Loesung.
#include <memory>
#include <vector>
#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/problems/scattering_problem.hpp"

namespace cbem {

struct LayeredGeometry {
    std::vector<TriangleMesh> surfaces;         // geschlossene, nach aussen orientierte Flaechen
    std::vector<int> inside, outside;           // Gebiet innen / aussen je Flaeche
    std::vector<Medium> region_medium;          // Medium je Gebiet; Gebiet 0 = Aussenraum
    explicit LayeredGeometry(Medium exterior = {}) : region_medium{exterior} {}
    int add_region(const Medium& m) { region_medium.push_back(m); return static_cast<int>(region_medium.size()) - 1; }
    int add_surface(const TriangleMesh& m, int in, int out) {
        surfaces.push_back(m); inside.push_back(in); outside.push_back(out); return static_cast<int>(surfaces.size()) - 1;
    }
    std::size_t triangles() const { std::size_t n = 0; for (auto& s : surfaces) n += s.size(); return n; }
};

// Schicht (Beschichtung) mit Dicke und Medium
struct Coating { real thickness; Medium medium; };

// Homogener Koerper im Gebiet parent (Standard: Aussenraum); gibt das Gebiet des Koerpers zurueck.
int add_body(LayeredGeometry& g, const TriangleMesh& surface, const Medium& inner, int parent = 0);

// Verschachtelte Flaechen, von aussen nach innen: media[l] ist das Medium innerhalb surfaces[l] (und ausserhalb
// surfaces[l+1]). Gibt die Gebiete der Schichten zurueck (letztes = Kern).
std::vector<int> add_layered_body(LayeredGeometry& g, const std::vector<TriangleMesh>& surfaces_outer_first,
                                  const std::vector<Medium>& media, int parent = 0);

// Beschichteter Koerper: coatings von innen (am Kern) nach aussen. outward = true: surface ist die Kernoberflaeche,
// die Schichten wachsen nach aussen (Parallelflaechen im Abstand d_1, d_1 + d_2, ...). outward = false: surface ist die
// Aussenflaeche des beschichteten Koerpers, die Schichten liegen innerhalb (verdraengen Kernmaterial, z. B. Oxidation).
std::vector<int> add_coated_body(LayeredGeometry& g, const TriangleMesh& surface, const Medium& core,
                                 const std::vector<Coating>& coatings, bool outward = true, int parent = 0);

// Operator des Systems (siehe oben); lebt nur zusammen mit den Gebietsoperatoren von LayeredScatteringProblem
class LayeredTransmissionOperator {
public:
    struct Region {
        const BoundaryOperator* E;               // Cauchy-Operator auf dem Gebietsrand (lokale Nummerierung)
        std::vector<std::size_t> tri;            // globaler Dreiecksindex je lokalem Dreieck
        std::vector<char> is_inner;              // 1: Gebiet liegt innen an dieser Flaeche (Spur J h, sigma = +1)
    };
    LayeredTransmissionOperator(std::vector<Region> regions, std::vector<Mat8> J);
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    void precondition(const std::vector<cplx>& x, std::vector<cplx>& y) const;   // 2 (1 + J)^{-1} je Dreieck
    std::size_t size() const { return 8 * J_.size(); }
    const std::vector<Mat8>& J() const { return J_; }
private:
    std::vector<Region> R_;
    std::vector<Mat8> J_, P_;
};

struct LayeredResult {
    real sigma_ext = 0;       // Extinktionsquerschnitt (optisches Theorem)
    cplx forward = 0;         // Vorwaertsamplitude S(0) = -i k conj(p).E_inf(d) / |p|^2 (Bohren-Huffman: S_1(0) = S_2(0))
    int iterations = 0; real residual = 0;
    std::vector<cplx> h;      // Aussenspuren aller Flaechen (Gesamtfeld)
};

// Voreinstellung der Eintragsauswertung fuer geschichtete Koerper: adaptive Nahfeldregel mit Randabstand
// (adapt_to_boundary), sonst wie EntryParams{} (siehe docs/results_coated.md)
inline EntryParams layered_entry_params() { EntryParams ep; ep.adapt_to_boundary = true; return ep; }

class LayeredScatteringProblem {
public:
    LayeredScatteringProblem(const LayeredGeometry& g, real omega, HMatrixParams hp = {}, EntryParams ep = layered_entry_params());
    LayeredResult solve_plane_wave(const Vec3& d, const CVec3& p, const SolveOptions& o = {}) const;
    const TriangleMesh& mesh() const { return all_.all; }             // alle Flaechen (globale Nummerierung)
    std::size_t surface_begin(std::size_t s) const { return all_.body_begin[s]; }
    const LayeredTransmissionOperator& T() const { return *T_; }
    double hmatrix_bytes() const;
    std::size_t near_pairs() const;
    double near_seconds() const;
private:
    LayeredGeometry g_;
    real omega_;
    MultiBodyMesh all_;
    std::vector<std::unique_ptr<TriangleMesh>> region_mesh_;
    std::vector<std::unique_ptr<KernelEntries>> ents_;
    std::vector<std::unique_ptr<KernelHMatrix>> hms_;
    std::vector<std::unique_ptr<CauchyOperator>> cops_;
    std::vector<std::unique_ptr<ChiralCauchyOperator>> chops_;
    std::unique_ptr<LayeredTransmissionOperator> T_;
    std::vector<std::size_t> ext_tri_;          // globale Dreiecke am Aussenraum (in der Reihenfolge von region_mesh_[0])
};

}  // namespace cbem
