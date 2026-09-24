#pragma once
// H-Matrix der vier Kernkomponenten K_c(i,j) (c = 0: Skalar, 1..3: Vektor) des Cauchy-Operators.
// Unzulaessige Blaetter: dicht (mit Nahfeld). Zulaessige Bloecke: ACA
//   Joint:         eine ACA auf der gestapelten Matrix |tau| x 4|sigma| (gemeinsame Zeilenfaktoren),
//   Componentwise: vier getrennte ACAs.
// Anwendung: Y_i += sum_c sum_j K_c(i,j) Z_{j,c}  mit Z: (N x 4 x 8) -> Y: (N x 8).
#include <vector>
#include "cbem/assembly/kernel_entries.hpp"
#include "cbem/hmatrix/aca.hpp"
#include "cbem/hmatrix/cluster_tree.hpp"

namespace cbem {

enum class AcaMode { Joint, Componentwise };

struct HMatrixParams {
    real eps = 1e-4;             // relative ACA-Toleranz je Block
    real eta = 1.0;              // Zulaessigkeit min(diam) <= eta dist
    std::size_t leaf = 32;       // maximale Blattgroesse
    real sep_factor = 3.0;       // zusaetzlich dist > sep_factor * h_max (Gauss-Fernfeldgenauigkeit)
    real max_kdiam = 20.0;       // |k| diam <= max_kdiam
    AcaMode mode = AcaMode::Joint;
};

struct HStats {
    std::size_t n_dense = 0, n_lowrank = 0;
    std::size_t entries_dense = 0, entries_lowrank = 0;   // gespeicherte komplexe Zahlen
    double mean_rank = 0; std::size_t max_rank = 0;
    double seconds = 0;
    double bytes() const { return 16.0 * (entries_dense + entries_lowrank); }
};

class KernelHMatrix {
public:
    KernelHMatrix(const KernelEntries& entries, HMatrixParams prm = {});
    void apply(const std::vector<cplx>& Z, std::vector<cplx>& Y) const;   // Z: N*4*8, Y: N*8 (wird addiert)
    const HStats& stats() const { return st_; }
    std::size_t size() const { return N_; }
private:
    struct Dense { std::vector<std::size_t> R, C; std::vector<KernelComp> K; };
    struct LR { std::vector<std::size_t> R, C; std::vector<LowRank> f; };   // Joint: f.size()==1 (V: 4C x r)
    void partition(int t, int s, std::vector<std::pair<int, int>>& adm, std::vector<std::pair<int, int>>& inadm) const;
    std::size_t N_;
    HMatrixParams prm_;
    ClusterTree tree_;
    cplx k_;
    std::vector<Dense> dense_;
    std::vector<LR> lr_;
    HStats st_;
};

}  // namespace cbem
