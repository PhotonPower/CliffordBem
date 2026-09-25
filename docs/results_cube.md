# Ergebnisse: Würfel mit Kanten-/Eck-Blockvorkonditionierung (C++-Kern v0.3)

## 1. Kompression auf gradierten Netzen: Kriterien für gestreckte Elemente

Das Tensornetz `make_cube_graded(L)` ist zu den Kanten geometrisch gradiert; die Elemente an
den Kanten sind stark gestreckt (lange Kante bis 0,5, kurze bis 2⁻ᴸ). Mit der Zulässigkeit
dist > 3 h_max waren diese Netze dicht-dominiert. Neu: `HMatrixParams::exact_in_lowrank` (ACA mit
exakten Einträgen, d. h. Nahfeldbehandlung wo nötig) erlaubt `sep_factor = 0`; die Zulässigkeit
hängt dann nur noch von η ab.

| L | Dreiecke | alt: dist > 3 h_max | neu: sep 0, exakte ACA | davon dicht | Fehler (neu) |
|---:|---:|---:|---:|---:|---:|
| 3 | 768 | 35 MB | 20 MB | 10 MB | 3,8·10⁻⁵ |
| 5 | 1 728 | 129 MB | 61 MB | 36 MB | 6,6·10⁻⁵ |
| 7 | 3 072 | 293 MB | 127 MB | 74 MB | 6,5·10⁻⁵ |
| 9 | 4 800 | – | 216 MB | 100 MB | 3,6·10⁻⁵ |

Der Speicher sinkt auf 43–57 %; der gradierte Würfel braucht jetzt weniger als der
gleichmäßige mit ähnlicher Dreieckszahl (4 800 Dreiecke gradiert: 216 MB; 6 912 gleichmäßig:
488 MB). Der Fehler bleibt unter ε = 10⁻⁴. Die Aufbauzeit bleibt hoch (L = 7: 64 s je Wellenzahl,
L = 9: 127 s), weil das Nahfeldkriterium D < 2,5 h_max bei gestreckten Elementen viele Paare
analytisch behandelt. Das ist der nächste Optimierungspunkt.

## 2. Streuung und Vorkonditionierung

`scatter_cube`: T₁ mit beiden Operatoren als H-Matrizen (neue Kriterien, ε = 10⁻⁴), ebene Welle,
ωa = 0,5 (a = halbe Kantenlänge), GMRES (Neustart 300, Toleranz 10⁻⁶), ohne Symmetriezerlegung.
„Blöcke“: exakte Inverse der Diagonalblöcke von T₁ auf allen Dreiecken mit Schwerpunktabstand < R
zu einer Ecke bzw. Kante (`BlockPreconditioner`, `group_by_features`), sonst 2(1+J)⁻¹.

### Gold, ε₁ = −11 + 1,2i

| L | Unbekannte | punktweise | Blöcke R = 0,125 | Blöcke R = 0,25 | σ_ext |
|---:|---:|---:|---:|---:|---:|
| 3 | 6 144 | 47 | 40 | 34 | 19,78547 |
| 5 | 13 824 | 54 | 41 | 34 | 19,55851 |
| 7 | 24 576 | 59 | 41 | – | 19,46273 |

### κ = −3 + 0,3i (innerhalb der kritischen Kurve der Kanten)

| L | Unbekannte | punktweise | Blöcke R = 0,25 | σ_ext |
|---:|---:|---:|---:|---:|
| 3 | 6 144 | 59 | 39 | 19,16417 |
| 5 | 13 824 | 80 | 40 | 19,08545 |
| 7 | 24 576 | 96 | 40 | 18,77837 |

- **Iterationen:** Mit Blöcken konstant (40–41 bzw. 34), punktweise wachsend (47 → 59 bzw.
  59 → 96). Das bestätigt den Prototyp (AP 2.5) im C++-Kern ohne Symmetriezerlegung.
- **Lösung:** Punktweise und blockweise vorkonditioniert stimmen auf 6–7 Stellen überein.
- **Genauigkeit:** Für Gold fällt σ_ext monoton (Differenzen 0,23; 0,10). Für −3 + 0,3i
  beschleunigt sich die Abnahme (0,08; 0,31), passend zum Resonanzfenster aus AP 1.
- **Blockanteil:** Auf dem Tensornetz liegen 75–94 % der Unbekannten in den Blöcken (größter
  Block bis 1 680 Unbekannte, Aufbau bis 57 s). Für große Streuer sind anisotrope Kantennetze
  nötig, bei denen der Anteil mit der Kantenlänge statt mit der Fläche wächst, sowie eine
  H-LU der Blöcke (AP 3.4).
