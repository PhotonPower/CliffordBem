#include "cbem/hmatrix/curved_hmatrix.hpp"

#include <algorithm>
#include <chrono>

namespace cbem {

namespace {
constexpr int C7 = kCurvedComps;

// acc[0..15] (8 komplexe Zahlen, Real- und Imaginaerteil abwechselnd) += k z[0..7] in reeller Arithmetik: dieselbe
// Rechnung wie std::complex (ac - bd, ad + bc), aber ohne die NaN-Pruefung mit __muldc3-Rueckfall, die GCC nach jeder
// komplexen Multiplikation erzeugt. acc ist ein lokales Feld (in Registern), Ergebnisse bitgleich zur alten Schleife (v0.54).
inline void load8(real* acc, const cplx* y) { const real* p = reinterpret_cast<const real*>(y); for (int i = 0; i < 16; ++i) acc[i] = p[i]; }
inline void store8(cplx* y, const real* acc) { real* p = reinterpret_cast<real*>(y); for (int i = 0; i < 16; ++i) p[i] = acc[i]; }
inline void axpy8(real* yr, cplx k, const cplx* z) {
    const real* zr = reinterpret_cast<const real*>(z);
    const real a = k.real(), b = k.imag();
    for (int r = 0; r < 8; ++r) {
        const real c = zr[2 * r], d = zr[2 * r + 1];
        yr[2 * r] += a * c - b * d;
        yr[2 * r + 1] += a * d + b * c;
    }
}
}

CurvedHMatrix::CurvedHMatrix(const CurvedKernelEntries& E, HMatrixParams prm)
    : N_(E.size()), prm_(prm), tree_(E.mesh().flat, std::max<std::size_t>(1, prm.leaf / 3)), k_(E.wavenumber()) {
    auto t0 = std::chrono::steady_clock::now();
    std::vector<std::size_t> perm; perm.reserve(3 * tree_.perm.size());
    for (std::size_t t : tree_.perm) for (std::size_t a = 0; a < 3; ++a) perm.push_back(3 * t + a);
    tree_.perm = std::move(perm);
    for (auto& nd : tree_.nodes) { nd.begin *= 3; nd.end *= 3; }
    std::vector<std::pair<int, int>> adm, inadm;
    partition(0, 0, adm, inadm);
    auto idx = [&](int node) {
        const ClusterNode& c = tree_.nodes[node];
        return std::vector<std::size_t>(tree_.perm.begin() + c.begin, tree_.perm.begin() + c.end);
    };
    dense_.resize(inadm.size());
    CBEM_OMP(omp parallel for schedule(dynamic))
    for (long b = 0; b < static_cast<long>(inadm.size()); ++b) {
        Dense& D = dense_[b];
        D.R = idx(inadm[b].first); D.C = idx(inadm[b].second);
        const std::size_t nc = D.C.size();
        D.K.resize(D.R.size() * nc);
        for (std::size_t a = 0; a < D.R.size(); a += 3)
            for (std::size_t c = 0; c < nc; c += 3) {
                const CurvedBlock Bk = E.block(D.R[a] / 3, D.C[c] / 3);
                for (int p = 0; p < 3; ++p) for (int q = 0; q < 3; ++q) D.K[(a + p) * nc + c + q] = Bk[p * 3 + q];
            }
    }
    lr_.resize(adm.size());
    CBEM_OMP(omp parallel for schedule(dynamic))
    for (long b = 0; b < static_cast<long>(adm.size()); ++b) {
        LR& B = lr_[b];
        B.R = idx(adm[b].first); B.C = idx(adm[b].second);
        const std::size_t m = B.R.size(), n = B.C.size();
        const bool ex = prm_.exact_in_lowrank;
        auto blk = [&](std::size_t ti, std::size_t tj) { return ex ? E.block(ti, tj) : E.block_far(ti, tj); };
        // Elementpaare, die schon im Zwischenspeicher der anderen Richtung stehen, werden von dort kopiert statt neu berechnet
        std::vector<std::vector<cplx>> rowc(m / 3), colc(n / 3);
        RowFn row = [&](std::size_t i, cplx* out) {
            std::vector<cplx>& rc = rowc[i / 3];
            if (rc.empty()) {
                rc.assign(3 * C7 * n, cplx(0));
                for (std::size_t j = 0; j < n; j += 3) {
                    const std::vector<cplx>& cc = colc[j / 3];
                    if (!cc.empty()) {
                        for (int p = 0; p < 3; ++p) for (int q = 0; q < 3; ++q) for (int c = 0; c < C7; ++c)
                            rc[p * C7 * n + c * n + j + q] = cc[(q * C7 + c) * m + i - i % 3 + p];
                        continue;
                    }
                    const CurvedBlock Bk = blk(B.R[i] / 3, B.C[j] / 3);
                    for (int p = 0; p < 3; ++p) for (int q = 0; q < 3; ++q) for (int c = 0; c < C7; ++c)
                        rc[p * C7 * n + c * n + j + q] = Bk[p * 3 + q][c];
                }
            }
            std::copy(rc.begin() + (i % 3) * C7 * n, rc.begin() + (i % 3 + 1) * C7 * n, out);
        };
        ColFn col = [&](std::size_t J, cplx* out) {
            const std::size_t c = J / n, j = J % n;
            std::vector<cplx>& cc = colc[j / 3];
            if (cc.empty()) {
                cc.assign(3 * C7 * m, cplx(0));
                for (std::size_t i = 0; i < m; i += 3) {
                    const std::vector<cplx>& rc = rowc[i / 3];
                    if (!rc.empty()) {
                        for (int p = 0; p < 3; ++p) for (int q = 0; q < 3; ++q) for (int c7 = 0; c7 < C7; ++c7)
                            cc[(q * C7 + c7) * m + i + p] = rc[p * C7 * n + c7 * n + j - j % 3 + q];
                        continue;
                    }
                    const CurvedBlock Bk = blk(B.R[i] / 3, B.C[j] / 3);
                    for (int p = 0; p < 3; ++p) for (int q = 0; q < 3; ++q) for (int c7 = 0; c7 < C7; ++c7)
                        cc[(q * C7 + c7) * m + i + p] = Bk[p * 3 + q][c7];
                }
            }
            std::copy(cc.begin() + ((j % 3) * C7 + c) * m, cc.begin() + ((j % 3) * C7 + c + 1) * m, out);
        };
        B.f = aca_select(prm_.aca_plus, row, col, m, C7 * n, prm_.eps);
        recompress(B.f, prm_.eps);
    }
    st_.n_dense = dense_.size(); st_.n_lowrank = lr_.size();
    std::size_t rsum = 0;
    for (auto& D : dense_) st_.entries_dense += C7 * D.K.size();
    for (auto& B : lr_) { st_.entries_lowrank += B.f.storage(); rsum += B.f.rank(); st_.max_rank = std::max(st_.max_rank, B.f.rank()); }
    st_.mean_rank = lr_.empty() ? 0 : double(rsum) / lr_.size();
    st_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void CurvedHMatrix::partition(int t, int s, std::vector<std::pair<int, int>>& adm, std::vector<std::pair<int, int>>& inadm) const {
    const ClusterNode& a = tree_.nodes[t]; const ClusterNode& b = tree_.nodes[s];
    real dist = box_distance(a, b), h = std::max(a.hmax, b.hmax);
    bool ok = std::min(a.diam, b.diam) <= prm_.eta * dist && dist > 0 && dist > prm_.sep_factor * h &&
              std::abs(k_) * std::max(a.diam, b.diam) <= prm_.max_kdiam;
    if (ok) { adm.emplace_back(t, s); return; }
    if (a.leaf() || b.leaf()) { inadm.emplace_back(t, s); return; }
    for (int x : a.child) for (int y : b.child) partition(x, y, adm, inadm);
}

void CurvedHMatrix::apply(const std::vector<cplx>& Z, std::vector<cplx>& Y) const {
    auto dense_block = [&](const Dense& D, std::vector<cplx>& Yo) {
        const std::size_t nc = D.C.size();
        for (std::size_t a = 0; a < D.R.size(); ++a) {
            cplx* y = &Yo[D.R[a] * 8];
            real acc[16]; load8(acc, y);
            for (std::size_t c = 0; c < nc; ++c) {
                const CurvedComp& K = D.K[a * nc + c];
                const cplx* z = &Z[D.C[c] * C7 * 8];
                for (int q = 0; q < C7; ++q) axpy8(acc, K[q], z + q * 8);
            }
            store8(y, acc);
        }
    };
    auto lr_block = [&](const LR& B, std::vector<cplx>& Yo, std::vector<cplx>& tmp) {
        const std::size_t m = B.R.size(), n = B.C.size(), r = B.f.rank();
        tmp.assign(r * 8, cplx(0));
        for (std::size_t k = 0; k < r; ++k) {
            const cplx* v = B.f.V.col(k);
            real acc[16]; load8(acc, &tmp[k * 8]);
            for (int c = 0; c < C7; ++c)
                for (std::size_t j = 0; j < n; ++j) axpy8(acc, v[c * n + j], &Z[B.C[j] * C7 * 8 + c * 8]);
            store8(&tmp[k * 8], acc);
        }
        for (std::size_t a = 0; a < m; ++a) {
            cplx* y = &Yo[B.R[a] * 8];
            real acc[16]; load8(acc, y);
            for (std::size_t k = 0; k < r; ++k) axpy8(acc, B.f.U(a, k), &tmp[k * 8]);
            store8(y, acc);
        }
    };
    const long nd = static_cast<long>(dense_.size()), nb = nd + static_cast<long>(lr_.size());
    if (omp_threads() == 1) {
        std::vector<cplx> tmp;
        for (long b = 0; b < nb; ++b) { if (b < nd) dense_block(dense_[b], Y); else lr_block(lr_[b - nd], Y, tmp); }
        return;
    }
    CBEM_OMP(omp parallel)
    {
        std::vector<cplx> Yl(Y.size(), cplx(0)), tmp;
        CBEM_OMP(omp for schedule(dynamic) nowait)
        for (long b = 0; b < nb; ++b) { if (b < nd) dense_block(dense_[b], Yl); else lr_block(lr_[b - nd], Yl, tmp); }
        CBEM_OMP(omp critical)
        for (std::size_t i = 0; i < Y.size(); ++i) Y[i] += Yl[i];
    }
}

}  // namespace cbem
