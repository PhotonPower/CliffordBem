#include "cbem/clifford/multivector.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    Multivector I = Multivector::blade(7), e1 = Multivector::blade(1), e2 = Multivector::blade(2);
    Multivector II = I * I; CHECK(std::abs(II.c[0] + 1.0) < 1e-15, "I^2 = %g", II.c[0].real());
    Multivector e12 = e1 * e2; CHECK(std::abs(e12.c[3] - 1.0) < 1e-15, "e1 e2 = e12");
    Multivector e21 = e2 * e1; CHECK(std::abs(e21.c[3] + 1.0) < 1e-15, "e2 e1 = -e12");
    Multivector iI = I * cplx(0, 1); Multivector s = iI * iI; CHECK(std::abs(s.c[0] - 1.0) < 1e-15, "(iI)^2 = 1");
    for (int b = 0; b < 8; ++b) { Multivector x = I * Multivector::blade(b) - Multivector::blade(b) * I; for (auto v : x.c) CHECK(std::abs(v) < 1e-15, "I zentral"); }
    // Linksmultiplikationsmatrix gegen Produkt
    Multivector a, x; for (int b = 0; b < 8; ++b) { a.c[b] = cplx(0.3 * b - 1, 0.1 * b); x.c[b] = cplx(1.0 / (b + 1), -0.2 * b); }
    Mat8 L = a.left_matrix(); cplx y[8]; apply(L, x.c.data(), y); Multivector ax = a * x;
    for (int b = 0; b < 8; ++b) CHECK(std::abs(y[b] - ax.c[b]) < 1e-13, "left_matrix");
    REPORT();
}
