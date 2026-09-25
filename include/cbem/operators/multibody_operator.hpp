#pragma once
// Innerer Randoperator fuer mehrere getrennte Koerper: E_1 = diag(E_1^(b)). Jeder Koerper hat seinen eigenen
// Innenraum (Medium, Wellenzahl) und damit seinen eigenen Cauchy-Operator auf dem eigenen Teilnetz; die
// Dreiecke von Koerper b liegen im Gesamtvektor im Bereich [begin[b], begin[b+1]).
// (Der aeussere Operator E_2 wirkt dagegen auf der Vereinigung aller Raender.)
#include <vector>
#include "cbem/operators/cauchy_operator.hpp"

namespace cbem {

class BlockDiagonalOperator : public BoundaryOperator {
public:
    BlockDiagonalOperator(std::vector<const BoundaryOperator*> blocks, std::vector<std::size_t> begin)
        : blocks_(std::move(blocks)), begin_(std::move(begin)) {}
    void apply(const std::vector<cplx>& x, std::vector<cplx>& y) const override {
        y.assign(x.size(), cplx(0));
        std::vector<cplx> xb, yb;
        for (std::size_t b = 0; b < blocks_.size(); ++b) {
            const std::size_t o = 8 * begin_[b], n = 8 * (begin_[b + 1] - begin_[b]);
            xb.assign(x.begin() + o, x.begin() + o + n);
            blocks_[b]->apply(xb, yb);
            std::copy(yb.begin(), yb.end(), y.begin() + o);
        }
    }
    std::size_t size() const override { return 8 * begin_.back(); }
private:
    std::vector<const BoundaryOperator*> blocks_;
    std::vector<std::size_t> begin_;
};

}  // namespace cbem
