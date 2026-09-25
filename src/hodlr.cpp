#include "cbem/solvers/hodlr.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>

namespace cbem {

LowRank block_aca(const EntryFn& A, const std::vector<std::size_t>& R, const std::vector<std::size_t>& C, real eps, std::size_t max_rank) {
    const std::size_t m = 8 * R.size(), n = 8 * C.size();
    std::vector<Matrix> Us, Vs;                  // je Schritt: U (m x q), V (n x q), Beitrag U V^T
    std::vector<char> used(R.size(), 0);
    std::size_t i = 0, rank = 0; real norm2 = 0;
    Matrix row(8, n), col(m, 8);
    for (std::size_t it = 0; it < R.size() && rank < max_rank; ++it) {
        // Blockzeile i (8 x n) minus bisherige Naeherung
        for (std::size_t b = 0; b < C.size(); ++b) { Mat8 E = A(R[i], C[b]); for (int r = 0; r < 8; ++r) for (int q = 0; q < 8; ++q) row(r, 8 * b + q) = E[r * 8 + q]; }
        for (std::size_t k = 0; k < Us.size(); ++k)
            for (std::size_t q = 0; q < Us[k].cols; ++q)
                for (int r = 0; r < 8; ++r) { cplx u = Us[k](8 * i + r, q); if (u == cplx(0)) continue; for (std::size_t c = 0; c < n; ++c) row(r, c) -= u * Vs[k](c, q); }
        used[i] = 1;
        // Pivot-Blockspalte: groesste Frobeniusnorm
        std::size_t jb = 0; real best = -1;
        for (std::size_t b = 0; b < C.size(); ++b) { real s = 0; for (int r = 0; r < 8; ++r) for (int q = 0; q < 8; ++q) s += std::norm(row(r, 8 * b + q)); if (s > best) { best = s; jb = b; } }
        if (best <= 1e-300) {
            auto it2 = std::find(used.begin(), used.end(), 0);
            if (it2 == used.end()) break;
            i = it2 - used.begin();
            continue;
        }
        // Blockspalte jb minus Naeherung
        for (std::size_t a = 0; a < R.size(); ++a) { Mat8 E = A(R[a], C[jb]); for (int r = 0; r < 8; ++r) for (int q = 0; q < 8; ++q) col(8 * a + r, q) = E[r * 8 + q]; }
        for (std::size_t k = 0; k < Us.size(); ++k)
            for (std::size_t q = 0; q < Us[k].cols; ++q)
                for (int c = 0; c < 8; ++c) { cplx v = Vs[k](8 * jb + c, q); if (v == cplx(0)) continue; for (std::size_t r = 0; r < m; ++r) col(r, c) -= Us[k](r, q) * v; }
        // Pivotblock P = row(:, jb) (8x8), Pseudoinverse ueber SVD
        Matrix P(8, 8); for (int r = 0; r < 8; ++r) for (int q = 0; q < 8; ++q) P(r, q) = row(r, 8 * jb + q);
        Matrix W, Z; std::vector<real> s; svd_jacobi(P, W, s, Z);
        std::size_t q = 0; while (q < 8 && s[q] > 1e-12 * s[0]) ++q;
        if (q == 0) {
            auto it2 = std::find(used.begin(), used.end(), 0);
            if (it2 == used.end()) break;
            i = it2 - used.begin();
            continue;
        }
        // A ~ col P^+ row = (col Z S^-1) (W^H row);  U = col Z S^-1 (m x q), V^T = W^H row  ->  V = row^T conj(W)
        Matrix U(m, q), V(n, q);
        for (std::size_t k = 0; k < q; ++k) {
            for (std::size_t r = 0; r < m; ++r) { cplx t = 0; for (int c = 0; c < 8; ++c) t += col(r, c) * Z(c, k); U(r, k) = t / s[k]; }
            for (std::size_t c = 0; c < n; ++c) { cplx t = 0; for (int r = 0; r < 8; ++r) t += std::conj(W(r, k)) * row(r, c); V(c, k) = t; }
        }
        real nu = frobenius(U), nv = frobenius(V); norm2 += nu * nu * nv * nv;
        Us.push_back(std::move(U)); Vs.push_back(std::move(V)); rank += q;
        if (nu * nv <= eps * std::sqrt(norm2)) break;
        // naechste Pivotzeile: groesste Norm in der neuen Spalte unter den unbenutzten Dreiecken
        std::size_t inext = R.size(); real bu = -1;
        for (std::size_t a = 0; a < R.size(); ++a) {
            if (used[a]) continue;
            real t = 0;
            for (int r = 0; r < 8; ++r) for (std::size_t k = 0; k < q; ++k) t += std::norm(Us.back()(8 * a + r, k));
            if (t > bu) { bu = t; inext = a; }
        }
        if (inext == R.size()) break;
        i = inext;
    }
    LowRank lr; lr.U = Matrix(m, rank); lr.V = Matrix(n, rank);
    std::size_t off = 0;
    for (std::size_t k = 0; k < Us.size(); ++k) {
        for (std::size_t q = 0; q < Us[k].cols; ++q) { std::copy(Us[k].col(q), Us[k].col(q) + m, lr.U.col(off + q)); std::copy(Vs[k].col(q), Vs[k].col(q) + n, lr.V.col(off + q)); }
        off += Us[k].cols;
    }
    if (rank > 0) recompress(lr, eps);
    return lr;
}

HodlrSolver::HodlrSolver(const TriangleMesh& m, const EntryFn& A, HodlrParams p) : p_(p) {
    auto t0 = std::chrono::steady_clock::now();
    ClusterTree t(m, p.leaf); perm_ = t.perm;
    nodes_.reserve(t.nodes.size());
    build(t, 0, A, m);
    sec_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

namespace {
// C = A^T B (A: k x r, B: k x c)
Matrix atb(const Matrix& A, const Matrix& B) {
    Matrix C(A.cols, B.cols);
    for (std::size_t i = 0; i < A.cols; ++i) for (std::size_t j = 0; j < B.cols; ++j) { cplx s = 0; const cplx* a = A.col(i); const cplx* b = B.col(j); for (std::size_t r = 0; r < A.rows; ++r) s += a[r] * b[r]; C(i, j) = s; }
    return C;
}
}  // namespace

int HodlrSolver::build(const ClusterTree& t, int cn, const EntryFn& A, const TriangleMesh& m) {
    const ClusterNode& c = t.nodes[cn];
    int id = static_cast<int>(nodes_.size()); nodes_.emplace_back();
    nodes_[id].begin = c.begin; nodes_[id].end = c.end;
    if (c.leaf()) {
        const std::size_t n = c.size(); Matrix D(8 * n, 8 * n);
        for (std::size_t a = 0; a < n; ++a) for (std::size_t b = 0; b < n; ++b) {
            Mat8 E = A(perm_[c.begin + a], perm_[c.begin + b]);
            for (int r = 0; r < 8; ++r) for (int q = 0; q < 8; ++q) D(8 * a + r, 8 * b + q) = E[r * 8 + q];
        }
        lu_factor(D, nodes_[id].piv); nodes_[id].LU = std::move(D); entries_ += 64 * n * n;
        return id;
    }
    int l = build(t, c.child[0], A, m), r = build(t, c.child[1], A, m);
    Node& N = nodes_[id]; N.child[0] = l; N.child[1] = r;
    const ClusterNode& c1 = t.nodes[c.child[0]]; const ClusterNode& c2 = t.nodes[c.child[1]];
    std::vector<std::size_t> I1(perm_.begin() + c1.begin, perm_.begin() + c1.end), I2(perm_.begin() + c2.begin, perm_.begin() + c2.end);
    N.B12 = block_aca(A, I1, I2, p_.eps); N.B21 = block_aca(A, I2, I1, p_.eps);
    max_rank_ = std::max({max_rank_, N.B12.rank(), N.B21.rank()});
    entries_ += N.B12.storage() + N.B21.storage();
    // Y1 = A11^{-1} U12, Y2 = A22^{-1} U21
    Matrix Y1 = N.B12.U, Y2 = N.B21.U;
    solve(l, Y1, 0); solve(r, Y2, 0);
    const std::size_t r1 = N.B12.rank(), r2 = N.B21.rank();
    // S = [[I, V12^T Y2], [V21^T Y1, I]]
    Matrix S(r1 + r2, r1 + r2);
    for (std::size_t k = 0; k < r1 + r2; ++k) S(k, k) = 1.0;
    Matrix S12 = atb(N.B12.V, Y2), S21 = atb(N.B21.V, Y1);
    for (std::size_t a = 0; a < r1; ++a) for (std::size_t b = 0; b < r2; ++b) S(a, r1 + b) = S12(a, b);
    for (std::size_t a = 0; a < r2; ++a) for (std::size_t b = 0; b < r1; ++b) S(r1 + a, b) = S21(a, b);
    if (r1 + r2 > 0) lu_factor(S, N.Spiv);
    N.S = std::move(S); N.Y1 = std::move(Y1); N.Y2 = std::move(Y2);
    entries_ += N.Y1.a.size() + N.Y2.a.size() + N.S.a.size();
    return id;
}

// B: Spalten sind rechte Seiten; die Zeilen row0 .. row0 + 8*(end-begin) gehoeren zum Knoten
void HodlrSolver::solve(int id, Matrix& B, std::size_t row0) const {
    const Node& N = nodes_[id];
    const std::size_t n = 8 * (N.end - N.begin);
    if (N.child[0] < 0) {
        std::vector<cplx> tmp(n);
        for (std::size_t j = 0; j < B.cols; ++j) { cplx* c = B.col(j) + row0; std::copy(c, c + n, tmp.begin()); lu_solve(N.LU, N.piv, tmp.data()); std::copy(tmp.begin(), tmp.end(), c); }
        return;
    }
    const Node& L = nodes_[N.child[0]];
    const std::size_t n1 = 8 * (L.end - L.begin);
    solve(N.child[0], B, row0); solve(N.child[1], B, row0 + n1);
    const std::size_t r1 = N.B12.rank(), r2 = N.B21.rank();
    if (r1 + r2 == 0) return;
    // t = Z^T x = [V12^T x2; V21^T x1];  t <- S^{-1} t;  x -= Y t  (x1 -= Y1 t1, x2 -= Y2 t2)
    std::vector<cplx> tv(r1 + r2);
    for (std::size_t j = 0; j < B.cols; ++j) {
        cplx* x = B.col(j) + row0;
        for (std::size_t k = 0; k < r1; ++k) { cplx s = 0; const cplx* v = N.B12.V.col(k); for (std::size_t q = 0; q < n - n1; ++q) s += v[q] * x[n1 + q]; tv[k] = s; }
        for (std::size_t k = 0; k < r2; ++k) { cplx s = 0; const cplx* v = N.B21.V.col(k); for (std::size_t q = 0; q < n1; ++q) s += v[q] * x[q]; tv[r1 + k] = s; }
        lu_solve(N.S, N.Spiv, tv.data());
        for (std::size_t k = 0; k < r1; ++k) { const cplx* y = N.Y1.col(k); cplx t = tv[k]; for (std::size_t q = 0; q < n1; ++q) x[q] -= y[q] * t; }
        for (std::size_t k = 0; k < r2; ++k) { const cplx* y = N.Y2.col(k); cplx t = tv[r1 + k]; for (std::size_t q = 0; q < n - n1; ++q) x[n1 + q] -= y[q] * t; }
    }
}

void HodlrSolver::apply(const std::vector<cplx>& b, std::vector<cplx>& x) const {
    const std::size_t N = perm_.size();
    Matrix B(8 * N, 1);
    for (std::size_t a = 0; a < N; ++a) for (int q = 0; q < 8; ++q) B(8 * a + q, 0) = b[8 * perm_[a] + q];
    solve(0, B, 0);
    x.resize(8 * N);
    for (std::size_t a = 0; a < N; ++a) for (int q = 0; q < 8; ++q) x[8 * perm_[a] + q] = B(8 * a + q, 0);
}

}  // namespace cbem
