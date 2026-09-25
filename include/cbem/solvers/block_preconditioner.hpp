#pragma once
// Kanten-/Eck-Blockvorkonditionierung (AP 2.5): exakte Inverse der Diagonalbloecke von T_1 auf
// Dreiecksgruppen nahe geometrischer Merkmale (Kanten, Ecken); sonst punktweise 2(1+J)^{-1}.
#include <utility>
#include <vector>
#include "cbem/linalg/dense.hpp"
#include "cbem/operators/transmission_operator.hpp"

namespace cbem {

struct FeatureSet {
    std::vector<Vec3> vertices;                          // Ecken
    std::vector<std::pair<Vec3, Vec3>> edges;            // Kanten als Strecken
    static FeatureSet cube(real a = 1.0);                // Wuerfel [-a,a]^3
};

// Gruppen: Dreiecke mit Schwerpunkt-Abstand < R zu einer Ecke (Vorrang) bzw. zur naechsten Kante
std::vector<std::vector<std::size_t>> group_by_features(const TriangleMesh& m, const FeatureSet& f, real R);

class BlockPreconditioner {
public:
    BlockPreconditioner(const TriangleMesh& m, const KernelEntries& E_inner, const KernelEntries& E_outer,
                        const TransmissionOperator& T, const std::vector<std::vector<std::size_t>>& groups);
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const;
    std::size_t max_block() const { return max_block_; }
    double fraction() const { return frac_; }       // Anteil der Unbekannten in Bloecken
    double seconds() const { return sec_; }
private:
    const TransmissionOperator& T_;
    std::vector<std::vector<std::size_t>> groups_;
    std::vector<Matrix> lu_;
    std::vector<std::vector<std::size_t>> piv_;
    std::size_t max_block_ = 0; double frac_ = 0, sec_ = 0;
};

}  // namespace cbem
