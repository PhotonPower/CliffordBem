#include "cbem/hmatrix/hmatrix.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#ifdef CBEM_USE_OPENMP
#include <omp.h>
#endif

namespace cbem {

KernelHMatrix::KernelHMatrix(const KernelEntries& E, HMatrixParams prm)
    : N_(E.mesh().size()), prm_(prm), tree_(E.mesh(), prm.leaf), k_(E.wavenumber()) {
    auto t0 = std::chrono::steady_clock::now();
    std::vector<std::pair<int, int>> adm, inadm;
    partition(0, 0, adm, inadm);
    dense_.resize(inadm.size()); lr_.resize(adm.size());
#ifdef CBEM_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (long b = 0; b < static_cast<long>(inadm.size()); ++b) {
        Dense& D = dense_[b];
        D.R = tree_.indices(tree_.nodes[inadm[b].first]); D.C = tree_.indices(tree_.nodes[inadm[b].second]);
        D.K.resize(D.R.size() * D.C.size());
        for (std::size_t a = 0; a < D.R.size(); ++a)
            for (std::size_t c = 0; c < D.C.size(); ++c) D.K[a * D.C.size() + c] = E.exact(D.R[a], D.C[c]);
    }
#ifdef CBEM_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (long b = 0; b < static_cast<long>(adm.size()); ++b) {
        LR& B = lr_[b];
        B.R = tree_.indices(tree_.nodes[adm[b].first]); B.C = tree_.indices(tree_.nodes[adm[b].second]);
        const std::size_t m = B.R.size(), n = B.C.size();
        if (prm_.mode == AcaMode::Joint) {
            RowFn row = [&](std::size_t i, cplx* out) {
                for (std::size_t j = 0; j < n; ++j) { KernelComp K = E.far(B.R[i], B.C[j]); for (int c = 0; c < 4; ++c) out[c * n + j] = K[c]; }
            };
            ColFn col = [&](std::size_t J, cplx* out) {
                std::size_t c = J / n, j = J % n;
                for (std::size_t i = 0; i < m; ++i) out[i] = E.far(B.R[i], B.C[j])[c];
            };
            LowRank f = aca_partial(row, col, m, 4 * n, prm_.eps); recompress(f, prm_.eps);
            B.f.push_back(std::move(f));
        } else {
            for (int c = 0; c < 4; ++c) {
                RowFn row = [&, c](std::size_t i, cplx* out) { for (std::size_t j = 0; j < n; ++j) out[j] = E.far(B.R[i], B.C[j])[c]; };
                ColFn col = [&, c](std::size_t j, cplx* out) { for (std::size_t i = 0; i < m; ++i) out[i] = E.far(B.R[i], B.C[j])[c]; };
                LowRank f = aca_partial(row, col, m, n, prm_.eps); recompress(f, prm_.eps);
                B.f.push_back(std::move(f));
            }
        }
    }
    st_.n_dense = dense_.size(); st_.n_lowrank = lr_.size();
    std::size_t rsum = 0, rcnt = 0;
    for (auto& D : dense_) st_.entries_dense += 4 * D.K.size();
    for (auto& B : lr_) for (auto& f : B.f) { st_.entries_lowrank += f.storage(); rsum += f.rank(); ++rcnt; st_.max_rank = std::max(st_.max_rank, f.rank()); }
    st_.mean_rank = rcnt ? double(rsum) / rcnt : 0;
    st_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void KernelHMatrix::partition(int t, int s, std::vector<std::pair<int, int>>& adm, std::vector<std::pair<int, int>>& inadm) const {
    const ClusterNode& a = tree_.nodes[t]; const ClusterNode& b = tree_.nodes[s];
    real dist = box_distance(a, b), h = std::max(a.hmax, b.hmax);
    bool ok = std::min(a.diam, b.diam) <= prm_.eta * dist && dist > prm_.sep_factor * h &&
              std::abs(k_) * std::max(a.diam, b.diam) <= prm_.max_kdiam;
    if (ok) { adm.emplace_back(t, s); return; }
    if (a.leaf() || b.leaf()) { inadm.emplace_back(t, s); return; }
    for (int x : a.child) for (int y : b.child) partition(x, y, adm, inadm);
}

void KernelHMatrix::apply(const std::vector<cplx>& Z, std::vector<cplx>& Y) const {
    // Dichte Bloecke
    for (auto& D : dense_) {
        const std::size_t n = D.C.size();
        for (std::size_t a = 0; a < D.R.size(); ++a) {
            cplx* y = &Y[D.R[a] * 8];
            for (std::size_t c = 0; c < n; ++c) {
                const KernelComp& K = D.K[a * n + c]; const cplx* z = &Z[D.C[c] * 32];
                for (int q = 0; q < 4; ++q) for (int r = 0; r < 8; ++r) y[r] += K[q] * z[q * 8 + r];
            }
        }
    }
    // Niedrigrangbloecke
    std::vector<cplx> tmp;
    for (auto& B : lr_) {
        const std::size_t m = B.R.size(), n = B.C.size();
        for (std::size_t fi = 0; fi < B.f.size(); ++fi) {
            const LowRank& f = B.f[fi]; const std::size_t r = f.rank();
            tmp.assign(r * 8, cplx(0));
            // tmp = V^T zs,  zs: Zeilen (c, j) bei Joint, (j) bei Componentwise mit c = fi
            for (std::size_t k = 0; k < r; ++k) {
                const cplx* v = f.V.col(k);
                if (B.f.size() == 1) {
                    for (int c = 0; c < 4; ++c) for (std::size_t j = 0; j < n; ++j) {
                        const cplx* z = &Z[B.C[j] * 32 + c * 8]; cplx vv = v[c * n + j];
                        for (int q = 0; q < 8; ++q) tmp[k * 8 + q] += vv * z[q];
                    }
                } else {
                    for (std::size_t j = 0; j < n; ++j) {
                        const cplx* z = &Z[B.C[j] * 32 + fi * 8]; cplx vv = v[j];
                        for (int q = 0; q < 8; ++q) tmp[k * 8 + q] += vv * z[q];
                    }
                }
            }
            for (std::size_t a = 0; a < m; ++a) {
                cplx* y = &Y[B.R[a] * 8];
                for (std::size_t k = 0; k < r; ++k) { cplx u = f.U(a, k); for (int q = 0; q < 8; ++q) y[q] += u * tmp[k * 8 + q]; }
            }
        }
    }
}

}  // namespace cbem
