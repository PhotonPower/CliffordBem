#include "cbem/hmatrix/curved_hmatrix.hpp"
#include "cbem/linalg/axpy8.hpp"

#include <algorithm>
#include <chrono>

namespace cbem {

namespace {
constexpr int C7 = kCurvedComps;

// Y += K Z fuer einen dichten Block; K: |R| x |C| Eintraege zu je 7 Komponenten (T = cplx oder complex<float>)
template <class T>
void dense_apply(const std::vector<std::size_t>& R, const std::vector<std::size_t>& C, const T* K, const std::vector<cplx>& Z, std::vector<cplx>& Yo) {
    const std::size_t nc = C.size();
    std::vector<cplx> buf;
    for (std::size_t a = 0; a < R.size(); ++a) {
        cplx* y = &Yo[R[a] * 8];
        real acc[16]; load8(acc, y);
        const cplx* row = row_as_double(K + a * nc * C7, nc * C7, buf);
        for (std::size_t c = 0; c < nc; ++c) {
            const cplx* k = row + c * C7;
            const cplx* z = &Z[C[c] * C7 * 8];
            for (int q = 0; q < C7; ++q) axpy8(acc, k[q], z + q * 8);
        }
        store8(y, acc);
    }
}

// Y += U V^T Z fuer einen niedrigrangigen Block; U: m x r, V: (7 n) x r, spaltenweise
template <class T>
void lr_apply(const std::vector<std::size_t>& R, const std::vector<std::size_t>& C, std::size_t r, const T* U, const T* V,
              const std::vector<cplx>& Z, std::vector<cplx>& Yo, std::vector<cplx>& tmp) {
    const std::size_t m = R.size(), n = C.size();
    tmp.assign(r * 8, cplx(0));
    for (std::size_t k = 0; k < r; ++k) {
        const T* v = V + k * C7 * n;
        real acc[16]; load8(acc, &tmp[k * 8]);
        for (int c = 0; c < C7; ++c)
            for (std::size_t j = 0; j < n; ++j) axpy8(acc, v[c * n + j], &Z[C[j] * C7 * 8 + c * 8]);
        store8(&tmp[k * 8], acc);
    }
    for (std::size_t a = 0; a < m; ++a) {
        cplx* y = &Yo[R[a] * 8];
        real acc[16]; load8(acc, y);
        for (std::size_t k = 0; k < r; ++k) axpy8(acc, U[k * m + a], &tmp[k * 8]);
        store8(y, acc);
    }
}

std::vector<std::complex<float>> to_float(const cplx* p, std::size_t n) {
    std::vector<std::complex<float>> f(n);
    for (std::size_t i = 0; i < n; ++i) f[i] = std::complex<float>(static_cast<float>(p[i].real()), static_cast<float>(p[i].imag()));
    return f;
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
    for (auto& B : lr_) B.r = B.f.rank();
    if (prm_.single_precision) {                                          // einfache Genauigkeit; double-Speicher freigeben
        for (auto& D : dense_) { D.Kf = to_float(D.K.empty() ? nullptr : D.K[0].data(), C7 * D.K.size()); std::vector<CurvedComp>().swap(D.K); }
        for (auto& B : lr_) { B.Uf = to_float(B.f.U.a.data(), B.f.U.a.size()); B.Vf = to_float(B.f.V.a.data(), B.f.V.a.size()); B.f = LowRank{}; }
        st_.entry_bytes = 8;
    }
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
    const bool sp = prm_.single_precision;
    auto dense_block = [&](const Dense& D, std::vector<cplx>& Yo) {
        if (sp) dense_apply(D.R, D.C, D.Kf.data(), Z, Yo);
        else if (!D.K.empty()) dense_apply(D.R, D.C, D.K[0].data(), Z, Yo);
    };
    auto lr_block = [&](const LR& B, std::vector<cplx>& Yo, std::vector<cplx>& tmp) {
        if (sp) lr_apply(B.R, B.C, B.r, B.Uf.data(), B.Vf.data(), Z, Yo, tmp);
        else lr_apply(B.R, B.C, B.r, B.f.U.a.data(), B.f.V.a.data(), Z, Yo, tmp);
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
