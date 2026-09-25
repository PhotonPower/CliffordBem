# Ergebnisse: Würfel mit Kanten-/Eck-Blockvorkonditionierung (C++-Kern v0.3, Nahfeld korrigiert in v0.5)

## 0. Korrektur in v0.5: Nahfeld auf gestreckten Elementen

Sauter-Schwab mit 5 Punkten je Richtung ist auf fast gleichseitigen Dreiecken sehr genau (Kugel,
Seitenverhältnis ≤ 1,45), auf gestreckten Dreiecken aber nicht. Gemessen gegen eine unabhängige
Referenz am gradierten Würfel (Seitenverhältnis h_max²/(2A) bis 16):

| Paar | Sauter-Schwab Ordnung 5 | Ordnung 16 | halbanalytisch, 14 Punkte |
|---|---:|---:|---:|
| gleiches Dreieck | 1,4·10⁻² | 2,3·10⁻⁵ | 8,9·10⁻⁵ |
| gemeinsame Kante | 6,7·10⁻³ | 6,9·10⁻⁶ | 4,4·10⁻⁵ |
| gemeinsame Ecke | 5,3·10⁻⁴ | 1,4·10⁻⁵ | 1,4·10⁻⁸ |

Bei gemeinsamer Kante hängt die Konvergenz zudem stark von der Orientierung der Kante ab (in der
ungünstigen Orientierung noch 38 % Fehler bei Ordnung 8). Ab v0.5 gilt deshalb:
- Sauter-Schwab (Ordnung 5) nur, wenn beide Dreiecke ein Seitenverhältnis ≤ 1,6 haben; sonst die
  **halbanalytische Regel**: inneres Integral exakt (Wilton), äußere Regel in Koordinaten, die
  kubisch zur gemeinsamen Kante bzw. Ecke gradiert sind (14 × 14 Punkte).
- Nahe, nicht benachbarte Paare: äußere Regel adaptiv (Halbierung der längsten Kante, bis der
  Umkreisradius jedes Teilstücks < 0,5 × Abstand ist). Abweichung zu feiner Referenz 5,5·10⁻⁶.
- Alle Nahpaare werden einmal vorab berechnet und gespeichert (mit der Symmetrie K_ji = (s, −v)).
- Robuste Kantenintegrale (asinh-Form) für Punkte nahe einer Kantengeraden.

Folgen: Die Würfelwerte aus v0.3/v0.4 enthielten Nahfeldfehler im Prozentbereich einzelner
Einträge. Die Iterationszahlen ändern sich praktisch nicht; σ_ext ändert sich bei L = 7 um 0,3 %
(Gold). Der Aufbau wird deutlich schneller (L = 7, beide Wellenzahlen: 45 s statt 125 s; Blöcke
13 s statt 36 s). Die Tabellen unten sind mit v0.5 neu gerechnet; die alten Werte liegen in
`results/cube_gold_v03.csv` und `results/cube_m3_v03.csv`.

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

| L | Unbekannte | punktweise | Blöcke R = 0,125 | σ_ext (v0.5) | σ_ext (v0.3) |
|---:|---:|---:|---:|---:|---:|
| 3 | 6 144 | 48 | 41 | 19,76256 | 19,78547 |
| 5 | 13 824 | 54 | 41 | 19,58598 | 19,55851 |
| 7 | 24 576 | 60 | 41 | 19,51713 | 19,46273 |

### κ = −3 + 0,3i (innerhalb der kritischen Kurve der Kanten)

| L | Unbekannte | punktweise | Blöcke R = 0,25 | σ_ext (v0.5) | σ_ext (v0.3) |
|---:|---:|---:|---:|---:|---:|
| 3 | 6 144 | 60 | 40 | 19,20251 | 19,16417 |
| 5 | 13 824 | 80 | 40 | 19,09556 | 19,08545 |
| 7 | 24 576 | 97 | 40 | 18,76510 | 18,77837 |

- **Iterationen:** Mit Blöcken konstant (40–41 bzw. 34), punktweise wachsend (47 → 59 bzw.
  59 → 96). Das bestätigt den Prototyp (AP 2.5) im C++-Kern ohne Symmetriezerlegung.
- **Lösung:** Punktweise und blockweise vorkonditioniert stimmen auf 6–7 Stellen überein.
- **Genauigkeit:** Für Gold fällt σ_ext monoton und konvergiert (Differenzen 0,18; 0,07; Quotient
  0,39). Für −3 + 0,3i beschleunigt sich die Abnahme auch mit korrigiertem Nahfeld (0,11; 0,33):
  Das ist kein Quadraturfehler, sondern passt zum Resonanzfenster aus AP 1.
- **Blockanteil:** Auf dem Tensornetz liegen 75–94 % der Unbekannten in den Blöcken (größter
  Block bis 1 680 Unbekannte, Aufbau bis 57 s). Für große Streuer sind anisotrope Kantennetze
  nötig, bei denen der Anteil mit der Kantenlänge statt mit der Fläche wächst, sowie eine
  H-LU der Blöcke (AP 3.4).

## 3. Aufwand des Nahfelds

Beim gradierten Würfel gilt jedes gestreckte Dreieck als nah zu vielen anderen (Kriterium
Schwerpunktabstand < 2,5 h_max; 518 Nahpaare je Dreieck bei L = 7). Das Kriterium lässt sich nicht
lockern: Die Gauß-7×7-Regel hat bei D/h_max ∈ [1,5; 2) bis 1,5·10⁻⁴ Fehler, erst ab 2,5 unter
3·10⁻⁶. Ein nicht benachbartes Nahpaar kostet 21 µs (Gauß 7×7: 1,7 µs). Die Vorberechnung dauert
bei L = 7 (3 072 Dreiecke) 19 s, bei L = 9 (4 800 Dreiecke) 93 s je Wellenzahl auf einem Kern; sie
ist mit OpenMP parallelisiert. L = 9 ist deshalb im Container (1 Kern, 300 s) nicht mehr
rechenbar.
