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
        const bool ex = prm_.exact_in_lowrank;
        auto ent = [&](std::size_t i, std::size_t j) { return ex ? E.exact(i, j) : E.far(i, j); };
        if (prm_.mode == AcaMode::Multivector) {
            // Kreuzapproximation ueber Cl3(C): Eintrag = Paravektor s + v
            auto entry = [&](std::size_t a, std::size_t b) { KernelComp K = ent(B.R[a], B.C[b]); Multivector M; M.c[0] = K[0]; M.c[1] = K[1]; M.c[2] = K[2]; M.c[4] = K[3]; return M; };
            std::vector<Multivector>& U = B.mu; std::vector<Multivector>& W = B.mw;   // Spalte k: U[k*m + a], W[k*n + b]
            std::vector<char> used(m, 0); std::size_t i = 0, r = 0; real S2 = 0;
            std::vector<Multivector> row(n), col(m);
            const std::size_t rmax = std::min(m, n);
            while (r < rmax) {
                for (std::size_t b = 0; b < n; ++b) { Multivector v = entry(i, b); for (std::size_t k = 0; k < r; ++k) v = v - U[k * m + i] * W[k * n + b]; row[b] = v; }
                used[i] = 1;
                // Pivot: groesste Norm unter den invertierbaren Eintraegen der Restzeile
                std::size_t j = n; real best = 0; Multivector Pinv;
                std::vector<std::size_t> order(n); for (std::size_t b = 0; b < n; ++b) order[b] = b;
                std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return mv_norm2(row[x]) > mv_norm2(row[y]); });
                for (std::size_t q = 0; q < std::min<std::size_t>(n, 8); ++q) {
                    bool ok; Multivector inv = mv_inverse(row[order[q]], &ok);
                    if (ok) { j = order[q]; best = mv_norm2(row[j]); Pinv = inv; break; }
                }
                if (j == n || best <= 1e-300) {
                    auto it = std::find(used.begin(), used.end(), 0); if (it == used.end()) break; i = it - used.begin(); continue;
                }
                for (std::size_t a = 0; a < m; ++a) { Multivector v = entry(a, j); for (std::size_t k = 0; k < r; ++k) v = v - U[k * m + a] * W[k * n + j]; col[a] = v; }
                real nu = 0, nw = 0;
                for (std::size_t a = 0; a < m; ++a) { U.push_back(col[a] * Pinv); nu += mv_norm2(U.back()); }
                for (std::size_t b = 0; b < n; ++b) { W.push_back(row[b]); nw += mv_norm2(row[b]); }
                ++r; S2 += nu * nw;
                if (std::sqrt(nu * nw) <= prm_.eps * std::sqrt(S2)) break;
                std::size_t inext = m; real bu = -1;
                for (std::size_t a = 0; a < m; ++a) if (!used[a]) { real t = mv_norm2(U[(r - 1) * m + a]); if (t > bu) { bu = t; inext = a; } }
                if (inext == m) break;
                i = inext;
            }
            B.mrank = r;
        } else if (prm_.mode == AcaMode::Joint) {
            RowFn row = [&](std::size_t i, cplx* out) {
                for (std::size_t j = 0; j < n; ++j) { KernelComp K = ent(B.R[i], B.C[j]); for (int c = 0; c < 4; ++c) out[c * n + j] = K[c]; }
            };
            ColFn col = [&](std::size_t J, cplx* out) {
                std::size_t c = J / n, j = J % n;
                for (std::size_t i = 0; i < m; ++i) out[i] = ent(B.R[i], B.C[j])[c];
            };
            LowRank f = aca_select(prm_.aca_plus, row, col, m, 4 * n, prm_.eps); recompress(f, prm_.eps);
            B.f.push_back(std::move(f));
        } else {
            for (int c = 0; c < 4; ++c) {
                RowFn row = [&, c](std::size_t i, cplx* out) { for (std::size_t j = 0; j < n; ++j) out[j] = ent(B.R[i], B.C[j])[c]; };
                ColFn col = [&, c](std::size_t j, cplx* out) { for (std::size_t i = 0; i < m; ++i) out[i] = ent(B.R[i], B.C[j])[c]; };
                LowRank f = aca_select(prm_.aca_plus, row, col, m, n, prm_.eps); recompress(f, prm_.eps);
                B.f.push_back(std::move(f));
            }
        }
    }
    st_.n_dense = dense_.size(); st_.n_lowrank = lr_.size();
    std::size_t rsum = 0, rcnt = 0;
    for (auto& D : dense_) st_.entries_dense += 4 * D.K.size();
    for (auto& B : lr_) {
        for (auto& f : B.f) { st_.entries_lowrank += f.storage(); rsum += f.rank(); ++rcnt; st_.max_rank = std::max(st_.max_rank, f.rank()); }
        if (prm_.mode == AcaMode::Multivector) { st_.entries_lowrank += 8 * (B.mu.size() + B.mw.size()); rsum += B.mrank; ++rcnt; st_.max_rank = std::max(st_.max_rank, B.mrank); }
    }
    st_.mean_rank = rcnt ? double(rsum) / rcnt : 0;
    st_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void KernelHMatrix::partition(int t, int s, std::vector<std::pair<int, int>>& adm, std::vector<std::pair<int, int>>& inadm) const {
    const ClusterNode& a = tree_.nodes[t]; const ClusterNode& b = tree_.nodes[s];
    real dist = box_distance(a, b), h = std::max(a.hmax, b.hmax);
    bool ok = std::min(a.diam, b.diam) <= prm_.eta * dist && dist > 0 && dist > prm_.sep_factor * h &&
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
        if (!B.mu.empty() || (prm_.mode == AcaMode::Multivector)) {
            // y_i += sum_k u_k(i) t_k,  t_k = sum_j w_k(j) z_j,  z_j = n_j x_j / sqrt|tau_j| = Z[j][c = 0]
            const std::size_t m = B.R.size(), n = B.C.size();
            for (std::size_t k = 0; k < B.mrank; ++k) {
                Multivector t;
                for (std::size_t b = 0; b < n; ++b) { Multivector z; for (int q = 0; q < 8; ++q) z.c[q] = Z[B.C[b] * 32 + q]; t = t + B.mw[k * n + b] * z; }
                for (std::size_t a = 0; a < m; ++a) { Multivector y = B.mu[k * m + a] * t; cplx* yy = &Y[B.R[a] * 8]; for (int q = 0; q < 8; ++q) yy[q] += y.c[q]; }
            }
            continue;
        }
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
