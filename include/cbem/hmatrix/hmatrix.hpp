#pragma once
// H-Matrix der vier Kernkomponenten K_c(i,j) (c = 0: Skalar, 1..3: Vektor) des Cauchy-Operators.
// Unzulaessige Blaetter: dicht (mit Nahfeld). Zulaessige Bloecke: ACA
//   Joint:         eine ACA auf der gestapelten Matrix |tau| x 4|sigma| (gemeinsame Zeilenfaktoren),
//   Componentwise: vier getrennte ACAs.
// Anwendung: Y_i += sum_c sum_j K_c(i,j) Z_{j,c}  mit Z: (N x 4 x 8) -> Y: (N x 8).
#include <vector>
#include "cbem/assembly/kernel_entries.hpp"
#include "cbem/assembly/linear_entries.hpp"
#include "cbem/hmatrix/aca.hpp"
#include "cbem/hmatrix/cluster_tree.hpp"
#include "cbem/clifford/multivector.hpp"

namespace cbem {

// Joint: eine skalare ACA auf den gestapelten Komponenten (gemeinsame Zeilenfaktoren);
// Componentwise: vier getrennte skalare ACAs;
// Multivector: Kreuzapproximation ueber der Algebra, K ~ sum_k u_k(i) w_k(j) (geometrisches Produkt) mit
//   Multivektor-Pivots K_ij^{-1}; Faktoren sind Multivektoren (8 Komponenten), ohne Nachkompression.
enum class AcaMode { Joint, Componentwise, Multivector };

struct HMatrixParams {
    real eps = 1e-4;             // relative ACA-Toleranz je Block
    real eta = 1.0;              // Zulaessigkeit min(diam) <= eta dist
    std::size_t leaf = 32;       // maximale Blattgroesse
    real sep_factor = 3.0;       // zusaetzlich dist > sep_factor * h_max (Gauss-Fernfeldgenauigkeit); 0 = aus
    bool exact_in_lowrank = false;  // ACA mit exakten Eintraegen (Nahfeld wo noetig) statt reiner Gauss-Fernfeldregel;
                                    // erlaubt sep_factor = 0 bei gestreckten Elementen
    real max_kdiam = 20.0;       // |k| diam <= max_kdiam
    AcaMode mode = AcaMode::Joint;
    bool single_precision = true;  // Eintraege nach dem Aufbau als complex<float> gespeichert (v0.56): halber Speicher; Fehler
                                   // auf sigma_ext 1e-12 bis 2e-8, 30- bis 2400-mal kleiner als die ACA-Kompression mit eps = 1e-4.
                                   // Aufbau und Produkt rechnen weiter in double (nur KernelHMatrix/CurvedHMatrix, nicht Multivector)
    bool aca_plus = true;        // ACA+ statt teilpivotisierter ACA (Modi Joint und Separate, Nahfeld; v0.30): gleiche Raenge,
                                 // etwa 10 % genauer, 2-6 % mehr Aufbauzeit, robust gegen unbemerkt zu kleinen Rang
};

struct HStats {
    std::size_t n_dense = 0, n_lowrank = 0;
    std::size_t entries_dense = 0, entries_lowrank = 0;   // gespeicherte komplexe Zahlen
    double mean_rank = 0; std::size_t max_rank = 0;
    double seconds = 0;
    double entry_bytes = 16;                                // 8 bei single_precision
    double bytes() const { return entry_bytes * (entries_dense + entries_lowrank); }
};

class KernelHMatrix {
public:
    KernelHMatrix(const KernelEntries& entries, HMatrixParams prm = {});
    // unstetig lineare Dichten (v0.45): Indizes I = 3 t + a (psi-Basis), Clusterbaum ueber Elementen mit je drei
    // zusammenhaengenden Indizes; dichte Bloecke und ACA-Zeilen/-Spalten je Elementpaar einmal ausgewertet; nur AcaMode::Joint
    KernelHMatrix(const LinearKernelEntries& entries, HMatrixParams prm = {});
    void apply(const std::vector<cplx>& Z, std::vector<cplx>& Y) const;   // Z: N*4*8, Y: N*8 (wird addiert)
    const HStats& stats() const { return st_; }
    std::size_t size() const { return N_; }
private:
    // Eintraege in double (K, f) oder nach dem Aufbau in einfacher Genauigkeit (Kf, ff; HMatrixParams::single_precision)
    struct FloatFactor { std::size_t r = 0; std::vector<std::complex<float>> U, V; };
    struct Dense { std::vector<std::size_t> R, C; std::vector<KernelComp> K; std::vector<std::complex<float>> Kf; };
    struct LR { std::vector<std::size_t> R, C; std::vector<LowRank> f;       // Joint: f.size()==1 (V: 4C x r)
                std::vector<FloatFactor> ff;                                  // wie f in einfacher Genauigkeit
                std::vector<Multivector> mu, mw; std::size_t mrank = 0; };     // Multivector: u (|R| x r), w (|C| x r), spaltenweise
    void partition(int t, int s, std::vector<std::pair<int, int>>& adm, std::vector<std::pair<int, int>>& inadm) const;
    template <class E> void build(const E& entries);
    std::size_t N_;
    HMatrixParams prm_;
    ClusterTree tree_;
    cplx k_;
    std::vector<Dense> dense_;
    std::vector<LR> lr_;
    HStats st_;
};

}  // namespace cbem
