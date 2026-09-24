# CliffordBem

Rechenkern für eine koordinatenfreie Randelementmethode (BEM) der Nano-Optik in der
komplexifizierten Clifford-Algebra Cl₃(ℂ). Grundlage ist die Dirac-Formulierung der
Maxwell-Gleichungen mit dem Faraday-Multivektor **F** = √ε **E** + I √μ **H** und der
resonanzfreien Transmissionsgleichung T₁ = E₂⁺ + E₁⁻ J (siehe `docs/papers`).

Stand 0.1: Kern für den **Cauchy-Randoperator E_k** mit H-Matrix-Kompression (ACA), gebaut
für den Asymptotiktest der Kompression (Hypothese H3 des Antrags). Die Architektur ist auf den
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

Die Kerneinträge stimmen mit dem Python-Prototyp (`prototype/ap2`) auf 10⁻¹³ überein.

## Kompressions-Benchmark

```bash
./build/bench_compression --geometry sphere --n 8,16,24,32 --k 1.5 --eps 1e-4 --mode joint --csv results/results_sphere.csv
python3 tools/analyze_compression.py results/results_sphere.csv
```

Optionen: `--geometry sphere|cube`, `--n` Liste (Kugel: 20 n² Dreiecke, Würfel: Gradierungsstufen),
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

Gemeinsame (`joint`) gegen komponentenweise (`comp`) ACA, Niedrigrang-Speicher bei gleicher
Genauigkeit: 25/31 MB (N = 1 280), 221/273 MB (5 120), 644/793 MB (11 520); das sind etwa
19 % weniger, bei halber Aufbauzeit.

## Aufbau

```
include/cbem/core        Grundtypen (Vec3, komplexe Zahlen, Kernkomponenten)
include/cbem/clifford    Cl3(C): Multivektoren, geometrisches Produkt, Linksmultiplikation
include/cbem/geometry    Dreiecksnetze (Kugel, gradierter Würfel), Quadraturregeln
include/cbem/kernel      Dirac-Fundamentallösung, analytische Dreiecksintegrale
include/cbem/assembly    Eintragsauswertung (Fern-/Nahfeld) = Schnittstelle zu H-Matrix und Lösern
include/cbem/linalg      kleine dichte Matrizen: QR, Jacobi-SVD
include/cbem/hmatrix     Clusterbaum, Blockpartition, ACA, H-Matrix
include/cbem/operators   Randoperatoren (Cauchy-Operator E_k)
apps/                    Benchmarks
tests/                   Tests (CTest)
tools/                   Auswertungsskripte
results/                 Benchmark-Ergebnisse (CSV)
prototype/               Python-Prototypen aus AP 1-3 (Referenz und Validierung)
docs/                    Architektur, Arbeitspapiere (AP 1-3), Zusammenfassung, Antrag
```

## Status und Ausbauplan

Siehe `docs/ARCHITECTURE.md`. Kurz: Als Nächstes folgen der Transmissionsoperator T₁ mit
Standardwahl von J, rechte Seiten und Fernfeld, GMRES und die Kanten-/Eck-Blockvorkonditionierung.
Danach Sauter-Schwab-Quadratur, Block-ACA mit Multivektor-Pivots, Netzimport (Gmsh) und
Python-Anbindung.

Die Beweise und Aussagen in den Arbeitspapieren sind vorläufig und nicht unabhängig geprüft.

## Lizenz

MIT, siehe `LICENSE`.
