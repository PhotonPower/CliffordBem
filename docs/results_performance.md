# Rechenzeit auf einem Kern: komplexe Arithmetik (v0.54, v0.55)

**Ausgangslage.** Bei 1 280 gekrümmten Elementen kostete das Lösen so viel wie der Aufbau. Je Iteration wird die
H-Matrix angewendet, und das Produkt las 481 MB mit 1,1 GB/s, bei einer Speicherbandbreite von 7,2 GB/s auf einem Kern
(Windows, MinGW g++ 16, `-O3` ohne `-march`). Die Messung steht in `results_curved.md` (v0.54).

**Ursache.** GCC rechnet `std::complex<double>`-Multiplikationen nach C99 Anhang G: Nach der Formel (ac − bd, ad + bc)
prüft es, ob das Ergebnis NaN + i·NaN ist, und springt dann in `__muldc3`, das Unendlichkeiten rettet. Für endliche Zahlen
ist das Ergebnis die Formel selbst; die Prüfung kostet aber eine Verzweigung je Multiplikation und verhindert die
Vektorisierung. Die Bibliothek enthielt 1 146 solche Aufrufstellen (`objdump -dr datei.obj | grep -c __muldc3`), u. a.
`twoport_layer_problem` 166, `dipole` 132, `recycling_gmres` 74, `KernelHMatrix` 45, `gmres` 27, HODLR 23.

**Änderungen.**
- v0.54, `CurvedHMatrix::apply`: innere Schleife y += k z reell in lokalen Akkumulatoren (`include/cbem/linalg/axpy8.hpp`);
  ohne die Prüfung, und das Ergebnis bleibt in Registern statt nach jeder Operation in den Speicher zu gehen.
- v0.55, alle Pfade: `-fcx-fortran-rules` für die Bibliothek (CMake-Option `CBEM_FAST_COMPLEX`, Voreinstellung an, nur
  wenn der Compiler den Schalter kennt): Multiplikation ohne NaN-Rückfall, Division mit Bereichsreduktion. Danach keine
  Aufrufstelle von `__muldc3` mehr. Dazu die Akkumulatoren aus v0.54 in `KernelHMatrix::apply` (konstante und lineare
  Dichten).

**Bitgleichheit.** Für endliche Zahlen rechnen beide Wege dieselbe Formel in derselben Reihenfolge. Gemessen (ein Kern,
tol = 10⁻¹⁰): σ_ext auf 17 Stellen und die Iterationszahlen sind in allen Fällen gleich.

| Fall (Gold, ε = −11 + 1,2i, ωa = 0,5) | σ_ext (alle Stände gleich) | Iterationen |
|---|---|---:|
| konstant, 5 120 Elemente | 1,8618421659175515 | 49 |
| konstant, 5 120 Elemente, HODLR | 1,8618421659940361 | 49 |
| linear, 1 280 Elemente | 1,8238575140081954 | 52 |
| gekrümmt, 320 Elemente | 1,8534737347487384 | 47 |

**Rechenzeit** (ein Kern):

| Fall | Aufbau v0.54 → v0.55 | Lösen v0.54 → v0.55 |
|---|---:|---:|
| konstant, 5 120 Elemente | 50,8 → 45,5 s | 26,4 → 18,8 s |
| konstant, 5 120 Elemente, HODLR | 463 → 347 s | 36,0 → 26,3 s |
| linear, 1 280 Elemente | 15,2 → 14,7 s | 21,3 → 15,3 s |
| gekrümmt, 320 Elemente | 3,9 → 3,9 s | 3,7 → 3,9 s (v0.54 schon umgestellt) |

Der Schalter allein brachte 7–25 % (Aufbau HODLR −25 %, Lösen linear −17 %), die Akkumulatoren im Produkt den Rest.
Die gekrümmten Elemente waren in v0.54 schon umgestellt: Lösen Gold 1 280 38,7 → 22,3 s.

**`-march=native`.** Mit `-DCBEM_NATIVE=ON` (AVX2, FMA) halbiert sich das H-Matrix-Produkt noch einmal, die dichten
Blöcke laufen dann nahe der Speicherbandbreite, und der Aufbau wird 15 % schneller. Die Ergebnisse sind dann wegen FMA
nicht mehr bitgleich. Die Option bleibt in der Voreinstellung aus, damit Programme und Wheels auf jedem Rechner laufen.

**Hinweise.**
- Laufzeiten immer für Aufbau *und* Lösen messen.
- Neue heiße Schleifen mit `objdump -dr` auf `__muldc3` prüfen, falls ohne `CBEM_FAST_COMPLEX` gebaut wird
  (z. B. andere Compiler).
- Nächster Schritt für das Produkt: Speicherung in einfacher Genauigkeit (complex64), weil die dichten Blöcke mit
  `-march=native` an der Speicherbandbreite liegen.
