# Rechenzeit und Speicher auf einem Kern: komplexe Arithmetik (v0.54, v0.55), einfache Genauigkeit (v0.56)

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
- Speicherung in einfacher Genauigkeit: folgender Abschnitt (v0.56).

## H-Matrix in einfacher Genauigkeit (v0.56)

**Messung der Genauigkeit** (Prototyp: Einträge nach dem Aufbau auf `complex<float>` gerundet, Rechnung in double;
tol = 10⁻¹⁰). Verglichen mit dem Fehler der ACA-Kompression (ε = 10⁻⁴ gegen 10⁻⁶):

| Fall (Gold) | float gegen double | ACA 10⁻⁴ gegen 10⁻⁶ |
|---|---:|---:|
| konstant, 1 280 Elemente | 1,6·10⁻⁹ | 5,1·10⁻⁷ |
| linear, 1 280 Elemente | 2·10⁻¹² | 4,8·10⁻⁶ |
| gekrümmt, 320 Elemente (Glas) | 2,3·10⁻⁸ | 4,7·10⁻⁸ |
| gekrümmt, 320 Elemente | 1,7·10⁻⁹ | 1,7·10⁻⁶ |
| gekrümmt, 1 280 Elemente | 3,7·10⁻⁹ | 1,5·10⁻⁶ |

Die Rundung ist 2- bis 2 400-mal kleiner als die Kompression; die Iterationszahlen ändern sich nicht.

**Änderung.** `HMatrixParams::single_precision` (Voreinstellung an): `KernelHMatrix` und `CurvedHMatrix` speichern dichte
Blöcke und die Faktoren U, V nach dem Aufbau als `complex<float>` und geben den double-Speicher frei; Aufbau und Produkt
rechnen in double (je Blockzeile werden die Einträge am Stück umgewandelt). `HStats::entry_bytes` = 8. Der Modus
`AcaMode::Multivector` behält seine Faktoren in double. Die eingebaute Variante gibt exakt die Werte des Prototyps.

**Wirkung** (ein Kern):

| Fall | Speicher double → float | s je Iteration, portabel | s je Iteration, `-march=native` |
|---|---:|---:|---:|
| konstant Gold 1 280 | 116 → 58 MB | 0,059 → 0,058 | 0,039 → 0,034 |
| linear Gold 1 280 | 566 → 283 MB | 0,286 → 0,288 | 0,187 → 0,163 |
| gekrümmt Gold 1 280 | 953 → 477 MB | 0,466 → 0,492 | 0,211 → 0,272 |

Der Speicher halbiert sich überall. Bei vier Kernkomponenten wird das Produkt mit `-march=native` 13 % schneller. Im
gekrümmten Pfad (sieben Komponenten) ist das dichte Produkt auch nativ durch die Rechnung begrenzt, nicht durch die
Bandbreite; dort kostet die Umwandlung 5 % (portabel) bzw. 20 % (nativ) Zeit des Produkts (H-Matrix allein: 0,101 → 0,122 s
nativ). Der halbe Speicher wiegt das auf; wer die Zeit braucht, setzt `single_precision = false`.

**Nebenbefund.** Mit `-march=native` weicht σ_ext vom portablen Bau bis 5·10⁻⁶ ab: FMA ändert die Rundung, und die
Pivotwahl der ACA verstärkt das bis zur Größe ihres Kompressionsfehlers (ε = 10⁻⁴). Float und double stimmen innerhalb
eines Baus auf 10⁻⁹ überein.

**Prüfung.** `test_hmatrix`: float gegen double für Joint und Componentwise 2,6·10⁻⁸, Speicher genau halbiert;
`test_curved` (8) ebenso für `CurvedHMatrix`.
