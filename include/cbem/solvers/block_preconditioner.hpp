#pragma once
// Blockvorkonditionierung: exakte Inverse der Diagonalbloecke von T_1 auf Dreiecksgruppen, sonst punktweise
// 2(1+J)^{-1}. Gruppen nach geometrischen Merkmalen (Kanten, Ecken; AP 2.5) oder nach Clustern (allgemeine
// Geometrien, mehrere Koerper). Der innere Block wird ueber eine Funktion geliefert, damit chirale und
// mehrere Koerper (blockdiagonaler Innenoperator) ohne Sonderfaelle behandelt werden.
#include <functional>
#include <utility>
#include <vector>
#include "cbem/linalg/dense.hpp"
#include "cbem/operators/transmission_operator.hpp"

namespace cbem {

struct FeatureSet {
    std::vector<Vec3> vertices;                          // Ecken
    std::vector<std::pair<Vec3, Vec3>> edges;            // Kanten als Strecken
    static FeatureSet cube(real a = 1.0, const Vec3& center = {0, 0, 0});   // Wuerfel center + [-a,a]^3
    void append(const FeatureSet& o) { vertices.insert(vertices.end(), o.vertices.begin(), o.vertices.end()); edges.insert(edges.end(), o.edges.begin(), o.edges.end()); }
};

// Gruppen: Dreiecke mit Schwerpunkt-Abstand < R zu einer Ecke (Vorrang) bzw. zur naechsten Kante
std::vector<std::vector<std::size_t>> group_by_features(const TriangleMesh& m, const FeatureSet& f, real R);
// Gruppen: Blaetter eines geometrischen Clusterbaums mit hoechstens max_size Dreiecken (ueberdecken alles)
std::vector<std::vector<std::size_t>> group_by_clusters(const TriangleMesh& m, std::size_t max_size);

using BlockFn = std::function<Matrix(const std::vector<std::size_t>&)>;   // E_1 auf B x B (8|B| x 8|B|)

class BlockPreconditioner {
public:
    BlockPreconditioner(const TriangleMesh& m, const KernelEntries& E_inner, const KernelEntries& E_outer,
                        const TransmissionOperator& T, const std::vector<std::vector<std::size_t>>& groups);
    BlockPreconditioner(const TriangleMesh& m, const BlockFn& inner_block, const KernelEntries& E_outer,
                        const TransmissionOperator& T, const std::vector<std::vector<std::size_t>>& groups);
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    std::size_t max_block() const { return max_block_; }
    double fraction() const { return frac_; }       // Anteil der Unbekannten in Bloecken
    double seconds() const { return sec_; }
private:
    void build(const TriangleMesh& m, const BlockFn& inner, const KernelEntries& E2);
    const TransmissionOperator& T_;
    std::vector<std::vector<std::size_t>> groups_;
    std::vector<Matrix> lu_;
    std::vector<std::vector<std::size_t>> piv_;
    std::size_t max_block_ = 0; double frac_ = 0, sec_ = 0;
};

}  // namespace cbem
