#pragma once
// Nahfeld im Aussenraum (v0.24).
//
// Das Streufeld an einem Punkt x ausserhalb aller Koerper ist das Cauchy-Integral der Streuspur h_s = h - h_inc auf der
// Flaeche zum Aussenraum:  F_s(x) = sign * sum_tau int_tau Phi_k(x - y) dS_y  n_tau u_tau  (u_tau Wert der stueckweise
// konstanten Dichte). Quadratur: fern die 7-Punkt-Regel; nahe (Abstand < 3 h) die Zerlegung Phi_k = Phi_0 + (-ik/(4 pi r))
// + Rest mit analytischen Dreiecksintegralen fuer die singulaeren Anteile (triangle_integrals) und der 7-Punkt-Regel fuer den
// glatten Rest -- genau bis dicht an die Oberflaeche. Chirales Aussenmedium: je Helizitaet mit k_pm und P_pm h_s.
// E = vec(F) / sqrt(eps), H aus dem Bivektoranteil (I sqrt(mu) H). Punkte innerhalb eines Koerpers (Windungszahl der
// Aussenflaeche > 1/2) werden markiert; dort liefert die Aussendarstellung nicht das physikalische Feld. Ebenso Punkte naeher als
// 2 % der mittleren Elementgroesse am Netz (Kantenspitzen der stueckweise konstanten Dichte; ab etwa 3 % der Elementgroesse
// trifft das Nahfeld die Mie-Loesung auf wenige Prozent).
#include <vector>
#include "cbem/hmatrix/aca.hpp"
#include "cbem/hmatrix/cluster_tree.hpp"
#include "cbem/hmatrix/hmatrix.hpp"
#include "cbem/sources/chiral_incidence.hpp"
#include "cbem/sources/incident_field.hpp"

namespace cbem {

struct NearFieldPoint {
    CVec3 E{}, H{};              // Gesamtfeld (einfallend + gestreut)
    bool inside = false;         // innerhalb eines Koerpers (Feld nicht definiert)
    bool too_close = false;      // naeher als 2 % der Elementgroesse am Netz: stueckweise konstante Dichten springen an den
                                 // Kanten, das Cauchy-Integral hat dort logarithmische Spitzen (Wert unzuverlaessig)
    real enhancement = 0;        // |E|^2 / |E_0|^2
    real chirality = 0;          // Im(conj(E).H) / |Im(conj(E_0).H_0)| (optische Chiralitaet relativ zur zirkularen ebenen Welle)
};

// Streufeld F_s an Punkten (Multivektor), Streuspur hs auf der Flaeche m, Wellenzahl k (ein Helizitaetskanal); direkte Summation
std::vector<Multivector> scattered_field(const TriangleMesh& m, const std::vector<cplx>& hs, cplx k, const std::vector<Vec3>& pts);

// Nahfeldoperator als rechteckige H-Matrix (v0.25): Zeilen = Auswertepunkte, Spalten = Dreiecke, Eintraege
// K(i, t) = int_tau Phi_k(x_i - y) dS (Skalar + Vektor, exakt wie in scattered_field). Zulaessige Bloecke (Punkt- und
// Dreieckscluster mit min(diam) <= eta dist, dist > sep_factor h, |k| diam <= max_kdiam) per ACA gemeinsam ueber die vier
// Komponenten, sonst dicht. Anwendung F_i = sum_c sum_t K_c(i, t) Z_{t,c} mit Z_{t,0} = n u, Z_{t,a} = e_a n u.
// Aufwand O((M + N) log) statt O(M N); lohnt ab einigen tausend Punkten.
class NearFieldOperator {
public:
    NearFieldOperator(const TriangleMesh& m, cplx k, const std::vector<Vec3>& pts, HMatrixParams prm = {});
    std::vector<Multivector> apply(const std::vector<cplx>& hs) const;   // F_s an den Punkten
    const HStats& stats() const { return st_; }
private:
    using Comp = std::array<cplx, 4>;
    struct Dense { std::vector<std::size_t> R, C; std::vector<Comp> K; };
    struct LR { std::vector<std::size_t> R, C; LowRank f; };
    void partition(int t, int s, std::vector<std::pair<int, int>>& adm, std::vector<std::pair<int, int>>& inadm) const;
    const TriangleMesh& m_;
    cplx k_;
    HMatrixParams prm_;
    ClusterTree rows_, cols_;
    std::vector<Dense> dense_;
    std::vector<LR> lr_;
    HStats st_;
};

// Optionen der Nahfeldauswertung: H-Matrix ab 'hmatrix_min_points' Punkten (0: immer, sehr gross: nie), ACA-Toleranz eps.
// eps = 1e-4 ergibt einen Fehler von etwa 1e-5 bezogen auf max|F| (weit unter dem Diskretisierungsfehler von etwa 1e-2);
// 28 800 Punkte x 5 760 Dreiecke: 17,6 s statt 78,8 s direkt, Anwendung allein 0,9 s (fuer wiederholte Auswertung).
struct NearFieldOptions { std::size_t hmatrix_min_points = 2000; real eps = 1e-4; };

// Gesamtes Nahfeld im Aussenraum: outer = Flaeche(n) zum Aussenraum, h = Gesamtspur darauf (wie im Loeser), Aussenmedium m.
// Die chirale Normierung bezieht sich auf die zirkulare ebene Welle derselben Richtung im selben Medium (auch bei linearer p).
std::vector<NearFieldPoint> exterior_near_field(const TriangleMesh& outer, const std::vector<cplx>& h, const Medium& m, real omega,
                                                const Vec3& d, const CVec3& p, const std::vector<Vec3>& pts, const NearFieldOptions& opt = {});

// Allgemeines einfallendes Feld (v0.37): b = Projektion des einfallenden Feldes auf outer (z. B. inc.project(outer, m) oder
// project_dipole), h_s = h - b; Gesamtfeld = Streufeld + inc.eval(x); Bezugsgroessen fuer Verstaerkung und Chiralitaet aus inc
std::vector<NearFieldPoint> exterior_near_field(const TriangleMesh& outer, const std::vector<cplx>& h, const std::vector<cplx>& b,
                                                const Medium& m, real omega, const IncidentField& inc, const std::vector<Vec3>& pts,
                                                const NearFieldOptions& o = {});

}  // namespace cbem
