#pragma once
// H-Matrix der sieben Kernkomponenten gekruemmter Elemente (v0.47, Stufe 2b): Indizes I = 3 t + a (psi-Basis), Clusterbaum
// ueber Elementen mit je drei zusammenhaengenden Indizes (wie KernelHMatrix fuer lineare Dichten), Joint-ACA ueber die
// gestapelten Komponenten, dichte Bloecke und ACA-Zeilen/-Spalten je Elementpaar einmal ausgewertet.
//   Y_I += sum_c sum_J K_c(I, J) Z_(J, c),   Z: N x 7 x 8 (Z_(J, c) = e_c x_J mit den Blades kCurvedBlade), Y: N x 8.
#include <vector>

#include "cbem/assembly/curved_entries.hpp"
#include "cbem/hmatrix/aca.hpp"
#include "cbem/hmatrix/cluster_tree.hpp"
#include "cbem/hmatrix/hmatrix.hpp"

namespace cbem {

// Voreinstellung fuer gekruemmte Elemente (v0.56): ACA-Toleranz 1e-5 statt 1e-4. Seit den genaueren Quadraturen (v0.52,
// v0.53) dominierte eps = 1e-4 den Fehler (sigma_ext Gold 1280: 1,5e-6, mit 1e-5 1,6e-7) bei etwa 5 % mehr Zeit und 15 %
// mehr Speicher (results_curved.md).
inline HMatrixParams curved_hmatrix_params() { HMatrixParams p; p.eps = 1e-5; return p; }

class CurvedHMatrix {
public:
    CurvedHMatrix(const CurvedKernelEntries& entries, HMatrixParams prm = curved_hmatrix_params());
    std::size_t size() const { return N_; }
    void apply(const std::vector<cplx>& Z, std::vector<cplx>& Y) const;
    const HStats& stats() const { return st_; }
private:
    // Eintraege in double (K, f) oder nach dem Aufbau in einfacher Genauigkeit (Kf, Uf, Vf; HMatrixParams::single_precision)
    struct Dense { std::vector<std::size_t> R, C; std::vector<CurvedComp> K; std::vector<std::complex<float>> Kf; };
    struct LR { std::vector<std::size_t> R, C; LowRank f; std::size_t r = 0; std::vector<std::complex<float>> Uf, Vf; };   // V: (7 n) x r
    std::size_t N_;
    HMatrixParams prm_;
    ClusterTree tree_;
    cplx k_;
    std::vector<Dense> dense_;
    std::vector<LR> lr_;
    HStats st_;
    void partition(int t, int s, std::vector<std::pair<int, int>>& adm, std::vector<std::pair<int, int>>& inadm) const;
};

}  // namespace cbem
