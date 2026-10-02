#pragma once
// Innere Schleifen der H-Matrix-Produkte (v0.54/v0.55): acc[0..15] (8 komplexe Zahlen, Real- und Imaginaerteil
// abwechselnd) += k z[0..7] in reeller Arithmetik. Dieselbe Rechnung wie std::complex (ac - bd, ad + bc), aber ohne die
// NaN-Pruefung mit __muldc3-Rueckfall, die GCC sonst nach jeder komplexen Multiplikation erzeugt; acc ist ein lokales Feld
// (in Registern). Die Reihenfolge der Additionen bleibt, die Ergebnisse sind bitgleich zur Schleife mit std::complex.
#include "cbem/core/types.hpp"

namespace cbem {

inline void load8(real* acc, const cplx* y) {
    const real* p = reinterpret_cast<const real*>(y);
    for (int i = 0; i < 16; ++i) acc[i] = p[i];
}

inline void store8(cplx* y, const real* acc) {
    real* p = reinterpret_cast<real*>(y);
    for (int i = 0; i < 16; ++i) p[i] = acc[i];
}

inline void axpy8(real* acc, cplx k, const cplx* z) {
    const real* zr = reinterpret_cast<const real*>(z);
    const real a = k.real(), b = k.imag();
    for (int r = 0; r < 8; ++r) {
        const real c = zr[2 * r], d = zr[2 * r + 1];
        acc[2 * r] += a * c - b * d;
        acc[2 * r + 1] += a * d + b * c;
    }
}

}  // namespace cbem
