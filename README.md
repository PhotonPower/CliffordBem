# CliffordBem

Rechenkern für eine koordinatenfreie Randelementmethode (BEM) der Nano-Optik in der
komplexifizierten Clifford-Algebra Cl₃(ℂ). Grundlage ist die Dirac-Formulierung der
Maxwell-Gleichungen mit dem Faraday-Multivektor **F** = √ε **E** + I √μ **H** und der
resonanzfreien Transmissionsgleichung T₁ = E₂⁺ + E₁⁻ J (siehe `docs/papers`).

Stand 0.2: **Cauchy-Randoperator E_k** mit H-Matrix-Kompression (ACA) und **Streulöser** für die
resonanzfreie Gleichung T₁ h = h_inc (Transmissionsabbildung J mit Standardwahl, ebene Welle,
Fernfeld/Extinktion, GMRES mit punktweiser Vorkonditionierung). Validiert gegen Mie bis
92 160 Unbekannte (`docs/results_scattering.md`).

Stand 0.7: **Mehrkörperprobleme und Gmsh-Import.** Beliebig viele getrennte Körper mit eigenem Medium
(auch chiral), `ScatteringProblem` als zentrale Einstiegsklasse, Netzimport aus Gmsh (ASCII 2.2/4.1) und
Beispielgeometrien. Beispiel: Zirkulardichroismus eines Born-Kuhn-Dimers aus Gold-Stäben
(`docs/results_multibody.md`).

Stand 0.6: **Chirale Medien und Zirkulardichroismus.** Pasteur-Medien im Inneren (Helizitätszerlegung,
chirale Transmissionsabbildung aus AP 1), zirkular polarisierte Anregung, chirale Mie-Referenz. Q₊, Q₋ und
CD treffen die Mie-Lösung nach Extrapolation auf 10⁻⁴ (`docs/results_chiral.md`).

Stand 0.5: **Nahfeld für gestreckte Elemente.** Sauter-Schwab ist auf gestreckten Dreiecken
ungenau (bis 1,4 % bei Seitenverhältnis 16); dort wird jetzt eine halbanalytische, zur gemeinsamen
Kante/Ecke gradierte Regel verwendet, für nahe nicht benachbarte Paare eine adaptive Außenregel.
Alle Nahpaare werden vorab berechnet und gespeichert. Die Würfelergebnisse sind neu gerechnet
(`docs/results_cube.md`).

Stand 0.4: **Sauter-Schwab-Quadratur** für Dreieckspaare mit gemeinsamer Fläche, Kante oder
Ecke. Damit konvergiert auch Gold mit Ordnung 2 gegen Mie (Extrapolation auf 2·10⁻⁵ statt
1,4·10⁻³ zuvor).

Stand 0.3: **Kanten-/Eck-Blockvorkonditionierung** (`BlockPreconditioner`) und **Kriterien für
gestreckte Elemente** (ACA mit exakten Einträgen, `sep_factor = 0`). Am gradierten Würfel sind
die GMRES-Iterationen damit unabhängig von der Kantenauflösung, und der Speicher sinkt auf
43–57 % (`docs/results_cube.md`). Die Architektur ist auf den
Ausbau zu einem vollständigen BEM-Löser ausgelegt (siehe `docs/ARCHITECTURE.md`).

## Bauen und Testen

Voraussetzungen: C++17-Compiler (getestet: g++ 13), CMake ≥ 3.16, optional OpenMP.
Keine weiteren Abhängigkeiten.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release     # Optionen: -DCBEM_OPENMP=ON -DCBEM_NATIVE=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Tests:
- `test_multivector`: Algebra (I² = −1, (iI)² = 1, Zentralität von I, Linksmultiplikationsmatrix).
- `test_triangle_integrals`: analytische Wilton-Integrale gegen feine Quadratur (Fehler < 10⁻¹⁰).
- `test_dense`: QR (CGS2) und einseitige Jacobi-SVD für komplexe Matrizen.
- `test_aca`: Kreuzapproximation eines getrennten Blocks, Fehler ∼ ε.
- `test_hmatrix`: H-Matrix-Produkt gegen dichtes Produkt, Fehler < 10 ε (beide ACA-Varianten).
- `test_plemelj`: Spurtrennung E h = ±h für innere und äußere Dirac-Lösungen (Konvergenz erster Ordnung).
- `test_gmres`: GMRES mit und ohne Vorkonditionierung, mit Neustart.
- `test_scattering`: Kugelstreuung, identisch mit der dichten Lösung des Python-Prototyps.
- `test_block_precond`: Blockvorkonditionierung (gleiche Lösung, weniger Iterationen; voller Block = exakte Inverse).
- `test_sauter_schwab`: Gewichtssummen, Exaktheit für glatte Integranden, exponentielle Konvergenz für
  die singulären Kernintegrale (gleiches Dreieck, gemeinsame Kante, gemeinsame Ecke).
- `test_nearfield`: Nahfeldregeln auf dem gradierten Würfel gegen feine Referenz; Nahfeld-Cache.
- `test_chiral`: chirale Transmissionsabbildung gegen den Prototyp, χ = 0 gleich achiral, Symmetrie σ_s(χ) = σ_{−s}(−χ).
- `test_multibody`: großer Abstand = Summe der Einzelquerschnitte, zwei unabhängige Formulierungen, Vertauschungssymmetrie.
- `test_gmsh`: Gmsh-Rundreise (2.2), Format 4.1, Orientierung, identische Streulösung.

Die Kerneinträge stimmen mit dem Python-Prototyp (`prototype/ap2`) auf 10⁻¹³ überein.

## Streurechnung an der Kugel

```bash
./build/scatter_sphere --n 8,12,16 --omega 1.0 --eps1 2.25,0 --csv results/scatter_glass.csv
./build/scatter_sphere --n 8,12,16 --omega 0.5 --eps1 -11,1.2 --csv results/scatter_gold.csv
cd tools && python3 analyze_scattering.py ../results/scatter_glass.csv   # Vergleich mit Mie
```

| Fall | Unbekannte | Q_ext | Mie | GMRES |
|---|---:|---:|---:|---:|
| Glas, ωa = 1 | 92 160 | 0,214745 (extrapoliert 0,215110) | 0,215098 | 10 |
| Gold, ωa = 0,5 | 92 160 | 0,591164 (extrapoliert 0,589967) | 0,590018 | 24 |

## Beliebige Geometrien aus Gmsh

```bash
pip install gmsh
python3 tools/make_geometries.py bornkuhn 0.12 examples/bornkuhn.msh --angle 60
./build/scatter_mesh --mesh examples/bornkuhn.msh --omega 0.5 --media "-11,1.2" --pol circ
```

Ein Medium je physikalischer Gruppe (`--media "eps1;eps2;..."`, `--chi "chi1;chi2;..."`). Für Kugel-Dimere
gibt es `scatter_multi`.

## Chirale Kugel (Zirkulardichroismus)

```bash
./build/scatter_chiral --n 8,12,16 --omega 1.0 --eps1 2.25,0 --chi 0.2,0 --csv results/chiral_lossless.csv
cd tools && python3 analyze_chiral.py ../results/chiral_lossless.csv     # Vergleich mit chiraler Mie-Lösung
```

## Würfel mit Kanten-/Eck-Blockvorkonditionierung

```bash
./build/scatter_cube --mesh graded --L 3,5,7 --omega 0.5 --eps1 -11,1.2 --R 0.125,0.25 --csv results/cube_gold.csv
```

| Gold, gradierter Würfel | L = 3 | L = 5 | L = 7 |
|---|---:|---:|---:|
| GMRES punktweise | 48 | 54 | 60 |
| GMRES mit Blöcken (R = 0,125) | 41 | 41 | 41 |

## Kompressions-Benchmark

```bash
./build/bench_compression --geometry sphere --n 8,16,24,32 --k 1.5 --eps 1e-4 --mode joint --csv results/results_sphere.csv
python3 tools/analyze_compression.py results/results_sphere.csv
```

Optionen: `--geometry sphere|cube|cube_uniform`, `--n` Liste (Kugel: 20 n² Dreiecke, `cube`:
Gradierungsstufen L mit 48 (L+1)² Dreiecken, `cube_uniform`: 12 n² Dreiecke), `--k-per-n c` (k = c·n),
`--k`/`--ki` Wellenzahl, `--eps` ACA-Toleranz, `--mode joint|comp`, `--leaf`, `--eta`, `--sep`,
`--check` Zahl der exakt geprüften Zeilen.

### Ergebnis (Kugel, k = 1,5, ε = 10⁻⁴, ein Kern, g++ -O3)

| Dreiecke | Unbekannte | Speicher gesamt | davon Niedrigrang | MB/(N log²N)·10³ | Aufbau | Mat-Vek | Fehler |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 280 | 10 240 | 55 MB | 25 MB | 0,83 | 6 s | 0,05 s | 1,1·10⁻⁵ |
| 2 880 | 23 040 | 169 MB | 91 MB | 0,93 | 16 s | 0,20 s | 9,0·10⁻⁶ |
| 5 120 | 40 960 | 340 MB | 221 MB | 0,91 | 33 s | 0,34 s | 9,8·10⁻⁶ |
| 8 000 | 64 000 | 662 MB | 342 MB | 1,02 | 52 s | 0,62 s | 7,0·10⁻⁶ |
| 11 520 | 92 160 | 959 MB | 644 MB | 0,95 | 80 s | 0,88 s | 6,9·10⁻⁶ |
| 15 680 | 125 440 | 1 474 MB | 866 MB | 1,01 | 113 s | 1,42 s | 5,9·10⁻⁶ |
| 20 480 | 163 840 | 1 881 MB | 1 388 MB | 0,93 | 156 s | 1,77 s | 8,2·10⁻⁶ |

Der Quotient Speicher/(N log² N) ist ab 2 880 Dreiecken innerhalb von ±7 % konstant; der
Speicher wächst also wie O(N log² N). Die lokalen Exponenten schwanken wegen der diskreten
Blockstruktur zwischen 0,9 und 1,5, im Mittel 1,23. Die Aufbauzeit wächst mit Exponent ≈ 1,15.
Eine dichte Systemmatrix hätte bei 20 480 Dreiecken etwa 430 GB.

Weitere Studien (gleichmäßiger Würfel mit Kanten und Ecken, Frequenzabhängigkeit bis kD = 48,
feste Elementzahl je Wellenlänge, gradierter Würfel) in `docs/results_compression.md`. Kurz:
Kanten und Ecken ändern das O(N log² N)-Wachstum nicht; bei konstanter Auflösung je Wellenlänge
kostet die Frequenz bis kD ≈ 16 nur etwa 10 %; stark gradierte Netze sind derzeit durch die
Nahfeldkriterien dicht-dominiert.

Gemeinsame (`joint`) gegen komponentenweise (`comp`) ACA, Niedrigrang-Speicher bei gleicher
Genauigkeit: 25/31 MB (N = 1 280), 221/273 MB (5 120), 644/793 MB (11 520); das sind etwa
19 % weniger, bei halber Aufbauzeit.

## Aufbau

```
include/cbem/core        Grundtypen (Vec3, komplexe Zahlen, Kernkomponenten)
include/cbem/clifford    Cl3(C): Multivektoren, geometrisches Produkt, Linksmultiplikation
include/cbem/geometry    Dreiecksnetze (Kugel, Würfel, Mehrkörper), Gmsh-Import, Quadratur, Sauter-Schwab
include/cbem/kernel      Dirac-Fundamentallösung, analytische Dreiecksintegrale
include/cbem/assembly    Eintragsauswertung (Fern-/Nahfeld) = Schnittstelle zu H-Matrix und Lösern
include/cbem/linalg      kleine dichte Matrizen: QR, Jacobi-SVD
include/cbem/hmatrix     Clusterbaum, Blockpartition, ACA, H-Matrix
include/cbem/operators   Randoperatoren (Cauchy-Operator E_k, chiraler Cauchy-Operator, Transmissionsoperator T_1)
include/cbem/solvers     GMRES
include/cbem/sources     ebene Wellen, Fernfeld, Extinktion
include/cbem/problems    ScatteringProblem: ein oder mehrere Körper, alle Operatoren, Lösung für ebene Wellen
apps/                    Benchmarks
tests/                   Tests (CTest)
tools/                   Auswertungsskripte
results/                 Benchmark-Ergebnisse (CSV)
prototype/               Python-Prototypen aus AP 1-3 (Referenz und Validierung)
docs/                    Architektur, Arbeitspapiere (AP 1-3), Zusammenfassung, Antrag
```

## Status und Ausbauplan

Siehe `docs/ARCHITECTURE.md`. Kurz: Als Nächstes folgen ein schnelleres Nahfeld für gestreckte
Elemente, anisotrope Kantennetze und H-LU für große Blöcke, danach Block-ACA mit Multivektor-Pivots, Netzimport (Gmsh) und
Python-Anbindung.

Die Beweise und Aussagen in den Arbeitspapieren sind vorläufig und nicht unabhängig geprüft.

## Lizenz

MIT, siehe `LICENSE`.
