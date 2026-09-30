#pragma once
// Adaptive Kreuzapproximation mit partieller Pivotisierung und QR/SVD-Nachkompression.
#include <functional>
#include "cbem/linalg/dense.hpp"

namespace cbem {

struct LowRank {                       // A ~ U V^T  (U: m x r, V: n x r)
    Matrix U, V;
    std::size_t rank() const { return U.cols; }
    std::size_t storage() const { return U.a.size() + V.a.size(); }
};

using RowFn = std::function<void(std::size_t i, cplx* out)>;   // Zeile i (Laenge n)
using ColFn = std::function<void(std::size_t j, cplx* out)>;   // Spalte j (Laenge m)

LowRank aca_partial(const RowFn& row, const ColFn& col, std::size_t m, std::size_t n, real eps, std::size_t rmax = 300);
// ACA+ (Grasedyck 2005; v0.30): zusaetzlich eine Referenzzeile und eine Referenzspalte, deren Reste nach jedem Schritt
// mitgefuehrt werden. Pivot aus dem groesseren Rest der beiden Referenzen; eine Referenz wird ersetzt, sobald sie Pivot war
// oder ihr Rest verschwindet. Abbruch erst, wenn das neue Kreuz |u||v| <= eps |S| ist UND die Referenzreste, auf den Block
// hochgerechnet, ebenfalls unter eps |S| liegen. So wird der ganze Block abgetastet; die teilpivotisierte ACA kann Teile nie
// sehen (Zeilen und Spalten, auf denen die bisherigen Kreuze verschwinden) und bricht dann unbemerkt mit zu kleinem Rang ab.
LowRank aca_plus(const RowFn& row, const ColFn& col, std::size_t m, std::size_t n, real eps, std::size_t rmax = 300);
// Auswahl nach Parameter
inline LowRank aca_select(bool plus, const RowFn& row, const ColFn& col, std::size_t m, std::size_t n, real eps, std::size_t rmax = 300) {
    return plus ? aca_plus(row, col, m, n, eps, rmax) : aca_partial(row, col, m, n, eps, rmax);
}
void recompress(LowRank& lr, real eps);

}  // namespace cbem
