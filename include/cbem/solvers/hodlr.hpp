#pragma once
// Hierarchische Faktorisierung des Systems T_1 im HODLR-Format (schwache Zulaessigkeit): Auf jeder Stufe des
// geometrischen Clusterbaums werden die beiden Nebendiagonalbloecke niedrigrangig approximiert
// (Block-ACA mit 8x8-Pivots, Nachkompression), die Diagonalbloecke rekursiv zerlegt, die Blaetter dicht per LU.
// Die Inverse wird ueber die Woodbury-Formel angewandt:
//   A = [[A11, U1 V1^T], [U2 V2^T, A22]] = D + W Z^T,
//   A^{-1} b = D^{-1} b - Y S^{-1} Z^T D^{-1} b,  Y = D^{-1} W,  S = I + Z^T Y.
// Kosten O(N r^2 log^2 N). Mit grober Toleranz als Vorkonditionierer fuer GMRES gedacht; mit feiner Toleranz
// ein direkter Loeser. (Klassische H-LU mit starker Zulaessigkeit ist nicht umgesetzt.)
#include <functional>
#include <memory>
#include <vector>
#include "cbem/clifford/multivector.hpp"
#include "cbem/hmatrix/aca.hpp"
#include "cbem/hmatrix/cluster_tree.hpp"
#include "cbem/linalg/dense.hpp"

namespace cbem {

using EntryFn = std::function<Mat8(std::size_t i, std::size_t j)>;   // 8x8-Block T_ij des Systems

// Block-ACA fuer eine Matrix aus 8x8-Bloecken (Zeilen-/Spaltendreiecke R, C): A ~ U V^T
LowRank block_aca(const EntryFn& A, const std::vector<std::size_t>& R, const std::vector<std::size_t>& C, real eps, std::size_t max_rank = 2000);

struct HodlrParams { std::size_t leaf = 64; real eps = 1e-3; };

class HodlrSolver {
public:
    HodlrSolver(const TriangleMesh& m, const EntryFn& A, HodlrParams p = {});
    void apply(const std::vector<cplx>& b, std::vector<cplx>& x) const;          // x ~ T^{-1} b
    double seconds() const { return sec_; }
    double bytes() const { return 16.0 * entries_; }
    std::size_t max_rank() const { return max_rank_; }
private:
    struct Node {
        std::size_t begin = 0, end = 0;          // Bereich in perm (Dreiecke), Unbekannte 8*(begin..end)
        int child[2] = {-1, -1};
        Matrix LU; std::vector<std::size_t> piv; // Blatt
        LowRank B12, B21;                        // T[t1,t2] ~ B12.U B12.V^T,  T[t2,t1] ~ B21.U B21.V^T
        Matrix Y1, Y2;                           // A11^{-1} B12.U, A22^{-1} B21.U
        Matrix S; std::vector<std::size_t> Spiv; // LU von S
    };
    int build(const ClusterTree& t, int cnode, const EntryFn& A, const TriangleMesh& m);
    void solve(int node, Matrix& B, std::size_t row0) const;     // B: Zeilen = Unbekannte des Knotens (ab row0 in B)
    std::vector<Node> nodes_;
    std::vector<std::size_t> perm_;
    HodlrParams p_;
    double sec_ = 0; std::size_t entries_ = 0, max_rank_ = 0;
};

}  // namespace cbem
