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

class CurvedHMatrix {
public:
    CurvedHMatrix(const CurvedKernelEntries& entries, HMatrixParams prm = {});
    std::size_t size() const { return N_; }
    void apply(const std::vector<cplx>& Z, std::vector<cplx>& Y) const;
    const HStats& stats() const { return st_; }
private:
    struct Dense { std::vector<std::size_t> R, C; std::vector<CurvedComp> K; };
    struct LR { std::vector<std::size_t> R, C; LowRank f; };   // V: (7 n) x r, komponentenweise gestapelt
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
