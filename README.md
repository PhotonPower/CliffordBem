# CliffordBem

Rechenkern für eine koordinatenfreie Randelementmethode (BEM) der Nano-Optik in der komplexifizierten
Clifford-Algebra Cl₃(ℂ). Grundlage ist die Dirac-Formulierung der Maxwell-Gleichungen mit dem
Faraday-Multivektor **F** = √ε **E** + I √μ **H** und der resonanzfreien Transmissionsgleichung
T₁ = E₂⁺ + E₁⁻ J (Theorie: `docs/papers`).

**Stand 0.15.** Galerkin-BEM mit stückweise konstanten Multivektor-Dichten auf ebenen Dreiecken:
H-Matrix-Kompression (ACA), Sauter-Schwab- und halbanalytische Nahfeldquadratur, achirale und chirale Medien,
mehrere Körper, beschichtete Grenzflächen (Kern-Schale, Mehrfachschichten, dünne Oxidschichten; exakt oder als
Dünnschicht-Näherung zweiter Ordnung auf einer Fläche), dispersive
Materialien, Spektren mit Orientierungsmittelung, Phase der Vorwärtsamplitude, Block- und hierarchische
Vorkonditionierung, Gmsh-Import. Versionsgeschichte: `CHANGELOG.md`.

## Was der Kern kann

| Bereich | Umfang | Validierung |
|---|---|---|
| Streuung an glatten Körpern | T₁ mit Standardwahl von J, ebene Wellen (linear/zirkular), Fernfeld, Extinktion | Kugel gegen Mie: Glas und Gold Ordnung 2, extrapoliert 10⁻⁵–10⁻⁴ (bis 92 160 Unbekannte) |
| Chirale Medien | Pasteur-Medien innen, Helizitätszerlegung, chirale J | chirale Mie-Lösung: Q₊, Q₋, CD auf 10⁻⁴ |
| Mehrere Körper | eigenes Medium je Körper, blockdiagonaler Innenoperator | Additivität bei großem Abstand, zwei unabhängige Formulierungen, Enantiomere mit entgegengesetztem CD |
| Beschichtete Grenzflächen | verschachtelte Gebiete (Grenzflächengraph), Parallelflächen, Schichten nach außen/innen, chirale Schichten, Vorwärtsamplitude S(0) mit Phase | beschichtete Kugel gegen Aden–Kerker: Q_ext und S(0) Ordnung 2; dünne Schichten (d/h bis 0,04) als Differenz zur neutralen Rechnung auf 1–5 % |
| Dünnschicht-Näherung | eine Fläche; zweite Ordnung in Dirac-Form (Schritte ∂_ν F = n(ik − D)F auf Parallelflächen, mit Formoperator), glatte Normalen; erste Ordnung als Sprungform | Schichtwirkung gegen Aden–Kerker: 0,1–0,7 % bis d/a = 0,05, −3,7 % bei 0,1 (erste Ordnung: 1,5 d/a) |
| Kanten und Ecken | gradierte Netze, Kanten-/Eckblöcke, Kriterien für gestreckte Elemente | Iterationen unabhängig von der Kantenauflösung (Würfel, Würfel-Dimer) |
| Kompression | ACA (gemeinsam, komponentenweise, Multivektor-Pivots) | Speicher O(N log² N) bis 163 840 Unbekannte; gemeinsame ACA am günstigsten |
| Spektren | Johnson-Christy Au/Ag, Hintergrundmedium, Lebedev-Mittelung | Goldkugel in Wasser gegen Mie (≤ 1,5 % bei 1 280 Dreiecken) |
| Vorkonditionierung | punktweise, Blöcke (Kanten/Ecken, Cluster), HODLR | voller Block bzw. feine HODLR-Toleranz = direkter Löser |
| Geometrie | Kugel, Würfel (gleichmäßig/gradiert), Mehrkörper, Gmsh 2.2/4.1 | Gmsh-Kugel konvergiert gegen Mie |

## Bauen und Testen

Voraussetzungen: C++17-Compiler (getestet: g++ 13), CMake ≥ 3.16, optional OpenMP. Keine weiteren Abhängigkeiten.
Für die Python-Werkzeuge: NumPy, SciPy, Matplotlib, optional `gmsh` (Geometrien; benötigt libGLU).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release     # Optionen: -DCBEM_OPENMP=ON -DCBEM_NATIVE=ON
cmake --build build -j
ctest --test-dir build --output-on-failure         # 18 Tests
```

| Test | prüft |
|---|---|
| `test_multivector` | Algebra (I² = −1, (iI)² = 1, Zentralität von I, Linksmultiplikation) |
| `test_triangle_integrals` | analytische Wilton-Integrale gegen feine Quadratur |
| `test_dense` | QR (CGS2), Jacobi-SVD |
| `test_aca` | Kreuzapproximation eines getrennten Blocks |
| `test_hmatrix` | H-Matrix-Produkt gegen dicht (gemeinsam, komponentenweise, Multivektor) |
| `test_plemelj` | Spurtrennung E h = ±h |
| `test_gmres` | GMRES mit/ohne Vorkonditionierung, Neustart |
| `test_scattering` | Kugelstreuung, Regression gegen den Python-Prototyp |
| `test_block_precond` | Blockvorkonditionierung, auch mehrere Körper mit chiralem Medium |
| `test_sauter_schwab` | Gewichtssummen, glatte Integranden, exponentielle Konvergenz der singulären Integrale |
| `test_nearfield` | Nahfeldregeln auf gestreckten Elementen, Nahfeld-Cache |
| `test_chiral` | chirale J gegen Prototyp, χ = 0 gleich achiral, Symmetrie σ_s(χ) = σ_{−s}(−χ) |
| `test_multibody` | großer Abstand, zwei Formulierungen, Vertauschungssymmetrie |
| `test_gmsh` | Rundreise 2.2, Format 4.1, Orientierung, gleiche Streulösung |
| `test_materials` | Materialtabellen, Lebedev-Momente |
| `test_hodlr` | Systemeinträge gegen Operator, HODLR als direkter Löser und Vorkonditionierer |
| `test_layered` | Parallelflächen, ohne Schicht = T₁, beschichtete Kugel gegen Aden–Kerker, neutrale Schale, chirale Schale |
| `test_thin_layer` | Flächenoperatoren und Formoperator auf der Kugel, ohne Schicht = T₁, neutrale Schicht ohne Wirkung, Schichtwirkung erster und zweiter Ordnung gegen Aden–Kerker |

## Anwendungen

| Programm | Zweck | Beispiel |
|---|---|---|
| `scatter_sphere` | Kugel, T₁, gegen Mie | `--n 8,12,16 --omega 1.0 --eps1 2.25,0` |
| `scatter_chiral` | chirale Kugel, Q₊/Q₋/CD | `--n 8,12 --omega 1.0 --eps1 2.25,0 --chi 0.2,0` |
| `scatter_cube` | gradierter Würfel, Kanten-/Eckblöcke | `--mesh graded --L 3,5 --eps1 -11,1.2 --R 0.125` |
| `scatter_cube_dimer` | Dimer aus gradierten Würfeln | `--L 3,5 --gap 0.5` |
| `scatter_multi` | Kugel-Dimer, auch chiral | `--n 8 --dist 3 --pol circ` |
| `scatter_mesh` | beliebige Gmsh-Geometrie | `--mesh stab.msh --media "-11,1.2" --pol circ --precond hodlr:1e-2` |
| `spectrum` | Spektren, Orientierungsmittelung, Beschichtungen | `--mesh x.msh --unit 25 --materials Ag --nbg 1.33 --lambda 340:520:20 [--coating "2:2.89,0"]` |
| `scatter_coated` | beschichtete Kugel/Gmsh-Körper, gegen Aden–Kerker, neutrale Referenz, Dünnschicht-Näherung | `--n 4,8 --omega 0.5 --core -11,1.2 --coat 0.02,2.25,0 --neutral` bzw. `--thin --bare` (`--thin-model jump\|dirac1\|dirac2`) |
| `bench_compression` | Asymptotik der Kompression | `--geometry sphere --n 8,16,24 --mode joint` |

Vorkonditionierung in `spectrum`, `scatter_mesh`, `scatter_multi`: `--precond point | cluster:G | hodlr:eps[:leaf]`.
Geometrien: `python3 tools/make_geometries.py sphere|rod|bornkuhn|roundcube h out.msh [--radius ρ] [--angle φ]`.
Auswertung: `tools/analyze_scattering.py`, `analyze_chiral.py`, `analyze_compression.py`, `mie_spectrum.py`,
`mie_coated.py` (Aden–Kerker, Mehrfachschichten), `analyze_coated.py`, `analyze_thin.py`, `plot_coated.py`.

## Ergebnisse

Übersicht in `docs/README.md`. Die wichtigsten Befunde:

- **Validierung:** Glas- und Goldkugel konvergieren mit Ordnung 2 gegen Mie; die chirale Kugel gegen die chirale
  Mie-Lösung (`results_scattering.md`, `results_chiral.md`).
- **Nahfeld:** Sauter-Schwab ist auf gestreckten Dreiecken ungenau (bis 1,4 %); dort wird eine halbanalytische,
  zur gemeinsamen Kante gradierte Regel verwendet (`results_cube.md`).
- **Kompression (H3):** Speicher O(N log² N); gemeinsame skalare ACA spart 19 % gegenüber komponentenweise, ACA mit
  Multivektor-Pivots braucht 2,3–2,9-mal mehr, weil der Kern skalar separabel ist (`results_compression.md`,
  `results_aca_multivector.md`).
- **Vorkonditionierung:** Kanten-/Eckblöcke machen die Iterationen unabhängig von der Kantenauflösung; bei glatten
  plasmonischen Körpern und an Resonanzen helfen lokale Blöcke und HODLR nur bei vielen rechten Seiten
  (`results_preconditioning.md`).
- **Resonanzfenster an Ecken:** in 2D mit Galerkin und komplexer Streckung in log r gelöst; in 3D an Kanten nicht
  übertragbar (Kernsingularitäten der Streckung), an Spitzen zulässig; Absorber und angereicherte Eckelemente
  nicht transparent (`results_resonance_window.md`). Abgerundete Kanten machen das Fenster rechenbar; die
  Ergebnisse hängen dann stark vom Rundungsradius ab (`results_roundcube.md`).
- **Beschichtete Grenzflächen:** exakte Schichtformulierung (jede Schichtgrenze eine Fläche), gegen Aden–Kerker mit
  Ordnung 2 in Extinktion und Phase. Bei Schichten dünner als die Elemente trägt die Rechnung einen systematischen
  Fehler, der sich in der Differenz zu einer neutralen Rechnung (Schicht aus Außenmedium, gleiche Netze) heraushebt;
  so ist die Wirkung einer Oxidschicht schon bei d/h ≈ 0,04 auf wenige Prozent genau. Für d ≪ h und d ≪ a gibt es
  zusätzlich eine Dünnschicht-Näherung erster Ordnung auf einer Fläche (Greensche Funktion der Schichtfolge in
  erster Ordnung, Ableitungen auf der Dichte statt im Kern), seit v0.15 in zweiter Ordnung mit Krümmung: Kosten wie ohne
  Schicht, keine neutrale Rechnung, Schichtwirkung bis d/a = 0,05 auf 0,1–0,7 % (Extinktion und Phase), Silberkugel mit
  2 nm Oxid so genau wie die exakte Rechnung. Referenzfläche auf der Metallseite (`results_coated.md`).
- **Anwendungen:** CD-Spektrum eines Born-Kuhn-Dimers, Silberwürfel-Spektren in Abhängigkeit vom Rundungsradius,
  Silberkugel und -würfel mit 2 nm Oxid: Rotverschiebung um 10 bzw. 20 nm, Phasenänderung der Vorwärtsamplitude bis
  0,9 bzw. 0,5 rad (`results_multibody.md`, `results_spectra.md`, `results_roundcube.md`, `results_coated.md`).

## Aufbau

```
include/cbem/core        Grundtypen, Materialmodelle
include/cbem/clifford    Cl3(C): Multivektoren, geometrisches Produkt, Inverse, Linksmultiplikation
include/cbem/geometry    Dreiecksnetze, Mehrkörper, Parallelflächen, Gmsh-Import, Quadratur, Sauter-Schwab
include/cbem/kernel      Dirac-Fundamentallösung, analytische Dreiecksintegrale
include/cbem/assembly    Eintragsauswertung (Fern-/Nahfeld, Nahfeld-Cache) = Schnittstelle zu H-Matrix und Lösern
include/cbem/linalg      dichte Matrizen: QR, Jacobi-SVD, LU
include/cbem/hmatrix     Clusterbaum, ACA-Varianten, H-Matrix
include/cbem/operators   Cauchy-Operator, chiraler und blockdiagonaler Innenoperator, T₁, dichte Blöcke
include/cbem/solvers     GMRES, Block-Vorkonditionierung, HODLR
include/cbem/sources     ebene Wellen, Fernfeld, Extinktion, Lebedev-Richtungen
include/cbem/problems    ScatteringProblem (Einstiegsklasse), LayeredScatteringProblem (beschichtete/geschichtete Körper),
                         ThinLayerScatteringProblem (Dünnschicht-Näherung erster Ordnung)
apps/  tests/  tools/  results/  examples/
data/materials/          Johnson-Christy Au, Ag (refractiveindex.info, CC0)
prototype/               Python-Prototypen: ap1 (Theorie), ap2 (3D-Galerkin, H-Matrix), resonance (Resonanzfenster)
docs/                    Architektur, Ergebnisberichte, Arbeitspapiere AP 1–3, Zusammenfassung, Antrag
```

## Offene Punkte

Siehe `docs/ARCHITECTURE.md`. Wichtigste: reflexionsfreier Abschluss an 3D-Kanten im Resonanzfenster,
Streckung um Spitzen im C++-Kern, Krylov-Recycling für viele rechte Seiten, Substrate, gekrümmte Elemente,
Ansätze höherer Ordnung, parallele Tests auf Mehrkernrechnern, Python-Anbindung; für beschichtete Körper:
Block-/HODLR-Vorkonditionierung, ein Eindeutigkeitsbeweis für verschachtelte Gebiete; für die Dünnschicht-Näherung
konsistente zweite Ableitungen (quadratische Anpassung), Stabilität für d ≳ h, mehrere Körper.

Die Beweise und Aussagen in den Arbeitspapieren sind vorläufig und nicht unabhängig geprüft.

## Lizenz

MIT, siehe `LICENSE`. Materialdaten: refractiveindex.info-Datenbank (CC0).
