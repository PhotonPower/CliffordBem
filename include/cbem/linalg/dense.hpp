#pragma once
// Kleine dichte komplexe Matrizen (spaltenweise) mit QR (CGS2) und einseitiger Jacobi-SVD.
// Ohne externe Abhaengigkeiten; ausgelegt fuer die Nachkompression von Niedrigrangfaktoren.
#include <vector>
#include "cbem/core/types.hpp"

namespace cbem {

struct Matrix {
    std::size_t rows = 0, cols = 0;
    std::vector<cplx> a;                               // spaltenweise: a[c*rows + r]
    Matrix() = default;
    Matrix(std::size_t r, std::size_t c) : rows(r), cols(c), a(r * c, cplx(0)) {}
    cplx& operator()(std::size_t r, std::size_t c) { return a[c * rows + r]; }
    const cplx& operator()(std::size_t r, std::size_t c) const { return a[c * rows + r]; }
    cplx* col(std::size_t c) { return &a[c * rows]; }
    const cplx* col(std::size_t c) const { return &a[c * rows]; }
};

// A = Q R  (Q: m x n mit orthonormalen Spalten, R: n x n obere Dreiecksmatrix), m >= n
void qr_cgs2(const Matrix& A, Matrix& Q, Matrix& R);
// A = U diag(s) V^H  (einseitige Jacobi-Rotationen), A: m x n, m >= n; s absteigend sortiert
void svd_jacobi(const Matrix& A, Matrix& U, std::vector<real>& s, Matrix& V, real tol = 1e-14, int max_sweeps = 60);
real frobenius(const Matrix& A);

}  // namespace cbem
