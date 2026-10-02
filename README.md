# CliffordBem

Rechenkern für eine koordinatenfreie Randelementmethode (BEM) der Nano-Optik in der komplexifizierten
Clifford-Algebra Cl₃(ℂ). Grundlage ist die Dirac-Formulierung der Maxwell-Gleichungen mit dem
Faraday-Multivektor **F** = √ε **E** + I √μ **H** und der resonanzfreien Transmissionsgleichung
T₁ = E₂⁺ + E₁⁻ J (Theorie: `docs/papers`).

**Stand 0.59.** Galerkin-BEM mit stückweise konstanten (seit v0.45 auch unstetig linearen) Multivektor-Dichten auf ebenen Dreiecken, seit v0.47 auch auf gekrümmten (quadratischen) Elementen:
H-Matrix-Kompression (ACA), Sauter-Schwab- und halbanalytische Nahfeldquadratur, achirale und chirale Medien,
mehrere Körper, beschichtete Grenzflächen (Kern-Schale, Mehrfachschichten, dünne Oxidschichten; exakt oder als
Dünnschicht-Näherung zweiter Ordnung auf einer Fläche, auch chirale Schichten; Zweitor-Formulierung für Mehrfachschichten
beliebiger Dicke), chirale Außenmedien (Teilchen in chiraler Lösung), dispersive
Materialien, Spektren mit Orientierungsmittelung, Nahfeldkarten (Feldverstärkung, optische Chiralität), Dipolanregung, Zerfallsraten und Fluoreszenzverstärkung von Emittern (Fluorophore), Phase der Vorwärtsamplitude, Block- und hierarchische
Vorkonditionierung, Gmsh-Import, Python-Anbindung für Skripte (`cliffordbem`, pybind11). Versionsgeschichte: `CHANGELOG.md`.

## Was der Kern kann

| Bereich | Umfang | Validierung |
|---|---|---|
| Streuung an glatten Körpern | T₁ mit Standardwahl von J, ebene Wellen (linear/zirkular), Fernfeld, Extinktion | Kugel gegen Mie: Glas und Gold Ordnung 2, extrapoliert 10⁻⁵–10⁻⁴ (bis 92 160 Unbekannte) |
| Chirale Medien | Pasteur-Medien innen, Helizitätszerlegung, chirale J | chirale Mie-Lösung: Q₊, Q₋, CD auf 10⁻⁴ |
| Mehrere Körper | eigenes Medium je Körper, blockdiagonaler Innenoperator | Additivität bei großem Abstand, zwei unabhängige Formulierungen, Enantiomere mit entgegengesetztem CD |
| Beschichtete Grenzflächen | verschachtelte Gebiete (Grenzflächengraph), Parallelflächen, Schichten nach außen/innen, chirale Schichten, Vorwärtsamplitude S(0) mit Phase | beschichtete Kugel gegen Aden–Kerker: Q_ext und S(0) Ordnung 2; dünne Schichten (d/h bis 0,04) als Differenz zur neutralen Rechnung auf 1–5 % |
| Dipolanregung | elektrischer Dipol als Quelle (beliebige rechte Seite `solve_rhs`), Zerfallsraten gesamt / strahlend / nicht strahlend, Quantenausbeute, konform verdichtetes Kugelnetz | Goldkugel gegen die Reihenlösung: strahlende Rate 1–3 %, Gesamtrate bei Elementgröße ∼ d/8 am Fußpunkt auf wenige Prozent |
| Nahfeld | Cauchy-Integral der Streuspur im Außenraum, halbanalytische Nahquadratur, rechteckige H-Matrix (Punkte × Dreiecke, ACA), Feldverstärkung und optische Chiralität, Karten (`nearfield`) | Goldkugel gegen Mie: Feldvektor auf 0,3–2 % (n = 12), auch 0,02 Radien vor der Oberfläche; Fernfeldgrenze 7,6·10⁻⁴ |
| Chirales Außenmedium | Helizitätswellen mit k_±, Außenoperator P₊E_{k₊} + P₋E_{k₋}, optisches Theorem je Kanal; in allen Formulierungen; seit v0.59 auch lineare und gekrümmte Elemente | Goldkugel in chiralem Wasser gegen Mie: CD Ordnung 2; Kugel aus dem Außenmedium unsichtbar |
| Zweitor (S-Matrix) | E₂ auf der Außenfläche, E₁ auf dem Kern, je Schicht ein Zweitor u_oben − u_unten = g(s)B(u_oben + u_unten) mit g = tanh(d√s/2)/√s (Partialbrüche, dünnbesetzte Resolventen), jede Schicht auf ihrer eigenen Fläche, automatische Unterteilung dicker Schichten | Mehrfachschichten gegen Mie auf 0,1–0,6 % (n = 12), dicke Schalen bis d/a = 0,5 konvergent, beliebiges d/h; CD chiraler Hüllen, auch im Abstand zum Gold; mehrere Körper (Dimer gegen exakte Methode auf 0,1 %) |
| Dünnschicht-Näherung | eine Fläche je Körper, auch mehrere Körper; zweite Ordnung in Dirac-Form (Schritte ∂_ν F = n(ik − D)F auf Parallelflächen, Formoperator, quadratische Anpassung der zweiten Ableitungen), glatte Normalen; erste Ordnung als Sprungform | Schichtwirkung gegen Aden–Kerker: 0,1–0,7 % bis d/a = 0,05; Dimer und Silberwürfel gegen die exakte Rechnung extrapoliert ≈ 1 %; chirale Schichten: CD gegen chirale Schicht-Mie-Lösung mit Ordnung 2 (0,5–1 % bei n = 12) |
| Kanten und Ecken | gradierte Netze, Kanten-/Eckblöcke, Kriterien für gestreckte Elemente | Iterationen unabhängig von der Kantenauflösung (Würfel, Würfel-Dimer) |
| Kompression | ACA (gemeinsam, komponentenweise, Multivektor-Pivots) | Speicher O(N log² N) bis 163 840 Unbekannte; gemeinsame ACA am günstigsten |
| Spektren | Johnson-Christy Au/Ag, Hintergrundmedium, Lebedev-Mittelung | Goldkugel in Wasser gegen Mie (≤ 1,5 % bei 1 280 Dreiecken) |
| Vorkonditionierung | punktweise, Blöcke (Kanten/Ecken, Cluster), HODLR | voller Block bzw. feine HODLR-Toleranz = direkter Löser |
| Geometrie | Kugel, Würfel (gleichmäßig/gradiert), Mehrkörper, Gmsh 2.2/4.1 (seit v0.48 auch zweiter Ordnung für gekrümmte Elemente) | Gmsh-Kugel konvergiert gegen Mie |
| Lineare Dichten (v0.45) | `LinearScatteringProblem`: unstetig lineare Dichten (24 Unbekannte je Dreieck), analytische Nahfeldintegrale, auch chirale und mehrere Körper; Grundlage gekrümmter Elemente | Einträge summieren sich exakt zu den konstanten (10⁻¹⁵), Kugel wie die unabhängige Vorhersage aus Stufe 1 (`docs/results_curved.md`) |
| Gekrümmte Elemente (v0.47) | `CurvedScatteringProblem`: quadratische Elemente (wie Gmsh `ElementOrder 2`) mit unstetig linearen Dichten, Normale im Integral (7 Kernkomponenten), Transmissionsabbildung als Elementmatrix , Nahfeld und optische Kräfte (v0.58), chirales Außenmedium (v0.59) | Kugel mit 320 Elementen: Fehler gegen Mie −0,012 % (Glas) und −0,007 % (Gold) statt −5,7 % und +7,3 %, Konvergenz etwa O(h⁴) (`docs/results_curved.md`) |
| Python-Anbindung | Paket `cliffordbem` (pybind11): alle Problemklassen, Anregungen, Nahfeld, Dipole, Kräfte, Algebra und die Bausteine (Einträge, H-Matrix, Cauchy-Operator, T₁, GMRES mit Python-Operatoren); eigene einfallende Felder als Python-Klassen (Bessel-, Wirbel-, Stehwellen); NumPy-Arrays, Rechnung ohne GIL | gleiche Ergebnisse wie C++ (Regression, Bausteine bitgleich), Nahfeld und Kräfte gegen Mie (`docs/python.md`) |

## Bauen und Testen

Voraussetzungen: C++17-Compiler (getestet: g++ 13), CMake ≥ 3.16, optional OpenMP. Keine weiteren Abhängigkeiten.
Für die Python-Werkzeuge: NumPy, SciPy, Matplotlib, optional `gmsh` (Geometrien; benötigt libGLU).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release     # Optionen: -DCBEM_OPENMP=ON -DCBEM_NATIVE=ON
cmake --build build -j
ctest --test-dir build --output-on-failure         # 29 Tests (30 mit -DCBEM_BUILD_PYTHON=ON)
```

Python-Anbindung (v0.42, `docs/python.md`): `pip install .` oder `cmake -S . -B build -DCBEM_BUILD_PYTHON=ON` (benötigt
pybind11 und die Python-Entwicklungsdateien; das Paket liegt dann in `build/python`):

```python
import cliffordbem as cb
P = cb.ScatteringProblem(cb.make_icosphere(8), cb.Medium(eps=-11 + 1.2j), omega=0.5)
print(P.solve_plane_wave(d=(0, 0, 1), p=(1, 0, 0)).sigma_ext / cb.pi)   # Q_ext der Goldkugel
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
| `test_gmsh` | Rundreise 2.2, Format 4.1, Orientierung, gleiche Streulösung; zweite Ordnung: echte Gmsh-Datei (Mitten auf der CAD-Kugel), Rundreise mit Umorientierung, Ablehnung erster Ordnung, Streuung gegen Mie |
| `test_materials` | Materialtabellen, Lebedev-Momente |
| `test_hodlr` | Systemeinträge gegen Operator, HODLR als direkter Löser und Vorkonditionierer |
| `test_layered` | Parallelflächen, ohne Schicht = T₁, beschichtete Kugel gegen Aden–Kerker, neutrale Schale, chirale Schale |
| `test_glmt` | Selbstprüfung der GLMT-Referenz (Python) |
| `test_beams` | Strahlen als Maxwell-Lösungen, Leistung, Impulserhaltung bei Dipolanregung, Teilchen im Fokus gegen Dipolnäherung |
| `test_optical_force` | Strahlungsdruck gegen Mie, drei Integrationswege, Dimer: Bindungskräfte und Summe; Dipolnäherung: ebene Welle, kleine und chirale Kugel gegen die volle Rechnung |
| `test_recycling` | Krylov-Recycling und GCRO-DR: gleiche Lösungen wie GMRES, weniger Iterationen bei Ausreißer-Eigenwerten; Eigenlöser |
| `test_dipole` | Dipol: Kugel aus dem Außenmedium (Raten = 1), Goldkugel gegen die Reihenlösung, verdichtetes Netz, Fluoreszenz, chiraler Emitter (Helizitätsdipol, magnetischer Dipol gegen Mie) |
| `test_near_field` | Nahfeld: Fernfeldgrenze, Goldkugel gegen Mie (\|E\|², Chiralität), Markierungen, χ → 0, Zweitor gegen Mie |
| `test_chiral_host` | chirales Außenmedium: Unsichtbarkeit, Ablehnung linearer Polarisation, χ → 0, Spiegelsymmetrie, Goldkugel gegen Mie, Zweitor gegen exakte Methode |
| `test_twoport` | Zweitor: neutrale Schicht, Glasschale d/h = 0,29 und 1,14 gegen Aden–Kerker, chirale Schicht (Symmetrie, CD), Mehrfachschichten (Zerlegung, Oxid + Glas, chirale Schicht auf Abstandshalter), Unterteilung dicker Schichten |
| `test_thin_layer` | Flächenoperatoren, Formoperator und Laplace–Beltrami auf der Kugel, ohne Schicht = T₁, neutrale Schicht ohne Wirkung, Schichtwirkung erster und zweiter Ordnung gegen Aden–Kerker (auch nach innen), mehrere Körper, chirale Schicht (Spiegelsymmetrie, CD gegen Mie, chiraler Kern = T₁) |
| `test_linear` | lineare Dichten: Innenintegrale gegen Brute Force, Summenidentität der Einträge in allen Zweigen, Projektion, Plemelj, Kugel gegen die Vorhersage aus Stufe 1, chirale Symmetrie, Mehrkörper |
| `test_curved` | gekrümmte Elemente: Geometrie (O(h⁴)), ebene Gegenproben der Einträge und des Streuproblems, Plemelj auf der gekrümmten Fläche, Kugel aus Glas und Gold gegen Mie, chirale Symmetrie, Nahquadratur mit Subtraktion, Fernblöcke in der ψ-Basis und schneller Kern |
| `test_curved_near_field` | Nahfeld und Strahlungsdruck auf gekrümmten Elementen gegen Mie (Konvergenzordnung), Fernfeldgrenze, allgemeine Projektion, Markierungen |
| `test_python` | Python-Anbindung (nur mit `-DCBEM_BUILD_PYTHON=ON`): Algebra, Netze, Regression und Bausteine bitgleich zu C++, Mie, chirale Medien, Nahfeld, Kräfte, Dipol, Strahlen, Schichten, Fehlerbehandlung, Threads ohne GIL; eigene einfallende Felder bitgleich zum Kern, Überlagerung, Maxwell-Prüfung |

## Anwendungen

| Programm | Zweck | Beispiel |
|---|---|---|
| `scatter_sphere` | Kugel, T₁, gegen Mie | `--n 8,12,16 --omega 1.0 --eps1 2.25,0` |
| `scatter_chiral` | chirale Kugel, Q₊/Q₋/CD | `--n 8,12 --omega 1.0 --eps1 2.25,0 --chi 0.2,0` |
| `scatter_cube` | gradierter Würfel, Kanten-/Eckblöcke | `--mesh graded --L 3,5 --eps1 -11,1.2 --R 0.125` |
| `scatter_cube_dimer` | Dimer aus gradierten Würfeln | `--L 3,5 --gap 0.5` |
| `scatter_multi` | Kugel-Dimer, auch chiral | `--n 8 --dist 3 --pol circ` |
| `scatter_mesh` | beliebige Gmsh-Geometrie | `--mesh stab.msh --media "-11,1.2" --pol circ --precond hodlr:1e-2` |
| `dipole` | Zerfallsraten und Fluoreszenzverstärkung eines Emitters (Kugel, Dimer, Gmsh; Abstandsreihe, Orientierungen x, y, z, Quantenausbeute, Anregung mit `--lambda-exc`) | `--sphere 16 --unit 20 --materials Au --nbg 1.33 --lambda 650 --dist "2,5,10,20" --q0 0.5` |
| `nearfield` | Nahfeldkarte in einer Ebene (Kugel, Dimer, Gmsh, Schichten, chirales Außenmedium); mit `--particle` Kraft und chirale Kraft auf ein kleines Teilchen | `--sphere 8 --sphere-dimer 4 --unit 20 --materials Au --nbg 1.33 --coating "1:2.25,0:0.01" --lambda 580 --pol circ --plane xz` |
| `tweezers` | Optische Pinzette: Kraft-Effizienz eines fokussierten Strahls (NA, Füllfaktor) oder Gaußstrahls auf eine Kugel, axial und seitlich | `--radius 0.25 --n-particle 1.59 --nbg 1.33 --lambda 1.064 --NA 1.2 --axial "-1:2:25"` |
| `spectrum` | Spektren, Orientierungsmittelung, Beschichtungen (auch chiral, Dünnschicht-Näherung, Zweitor mit `--twoport`) | `--mesh x.msh --unit 25 --materials Ag --nbg 1.33 --lambda 340:520:20 [--coating "2:2.89,0[:χ]" --thin 0]` |
| `scatter_coated` | beschichtete Kugel/Gmsh-Körper, gegen Aden–Kerker, neutrale Referenz, Dünnschicht-Näherung | `--n 4,8 --omega 0.5 --core -11,1.2 --coat 0.02,2.25,0 --neutral` bzw. `--thin --bare` (`--thin-model jump\|dirac1\|dirac2\|dirac2fit`, auch mit `--mesh`) |
| `bench_compression` | Asymptotik der Kompression | `--geometry sphere --n 8,16,24 --mode joint` |

Vorkonditionierung in `spectrum`, `scatter_mesh`, `scatter_multi`: `--precond point | cluster:G | hodlr:eps[:leaf]`.
Geometrien: `python3 tools/make_geometries.py sphere|rod|bornkuhn|roundcube h out.msh [--radius ρ] [--angle φ]`.
Auswertung: `tools/analyze_scattering.py`, `analyze_chiral.py`, `analyze_compression.py`, `mie_spectrum.py`,
`mie_coated.py` (Aden–Kerker, Mehrfachschichten), `mie_chiral_layered.py` (geschichtete Kugeln mit chiralen Schichten),
`analyze_coated.py`, `analyze_thin.py`, `analyze_chiral_thin.py`, `analyze_twoport.py`, `plot_coated.py`,
`mie_nearfield.py` (Nahfeld nach Mie), `plot_nearfield.py`, `mie_dipole.py` (Zerfallsraten vor der Kugel), `plot_fluorescence.py`, `mie_force.py` (Strahlungsdruck), `mie_polarizability.py` (Dipol-Polarisierbarkeiten, auch chiral), `glmt.py` (Kugel im fokussierten Strahl), `plot_tweezers.py`.

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
  2 nm Oxid so genau wie die exakte Rechnung; seit v0.16 punktweise konsistente zweite Ableitungen und mehrere Körper.
  Referenzfläche auf der Metallseite; die Krümmung muss auf dem Netz aufgelöst sein (Silberwürfel mit groben Rundungen:
  5–13 % Unterschied zur exakten Rechnung, extrapoliert ≈ 1 %). Seit v0.17 auch chirale Schichten und Kerne: Der CD
  einer dünnen chiralen Hülle auf Gold konvergiert mit Ordnung 2 gegen eine neue Mie-Lösung für geschichtete chirale
  Kugeln, während die exakte Zwei-Flächen-Rechnung ihn bei d/h ≲ 0,4 um 30–50 % verfehlt. Spektrum einer Goldkugel mit
  1 nm chiraler Molekülhülle in Wasser: plasmoninduzierter CD, an der Resonanz auf 2–7 %. Seit v0.19 (Prototyp) die
  Schicht als Zweitor nach dem Vorbild der S-Matrix von Schichtsystemen: stabil für Schichten dicker als die Elemente
  (Dünnschicht-Näherung bei d/a = 0,2: −11 %, 254 Iterationen; Zweitor: −0,3 %). Seit v0.20 Mehrfachschichten, jede
  Schicht auf ihrer eigenen Fläche und dicke Schichten automatisch unterteilt: konvergent bis d/a = 0,5; CD einer
  chiralen Schicht im Abstand 0–8 nm zu einer Goldkugel auf einen konstanten Versatz von 0,015 nm² genau. Seit v0.21 auch
  chirale Außenmedien: Teilchen in chiraler Lösung, Beiträge freier und gebundener Moleküle zum CD. Seit v0.22 mehrere
  Körper im Zweitor: Gold-Dimer mit chiraler Molekülschicht, Verstärkung des CD im Spalt. Seit v0.23 lehnen alle Verfahren
  sich berührende oder durchdringende Körper und Hüllen ab (`results_coated.md`).
- **Nahfeld (v0.24):** Feldverstärkung und optische Chiralität im Außenraum, gegen Mie auf wenige Prozent bis dicht an die
  Oberfläche. Im Spalt des Gold-Dimers bis 2 000-fache Intensität und zehnfache optische Chiralität; das Spaltfeld
  verlangt feinere Netze als die Spektren. Seit v0.25 mit H-Matrix-Auswertung: 4,5-fach schneller bei 28 800 Punkten,
  Kompressionsfehler 10⁻⁵ (`results_nearfield.md`).
- **Dipolanregung (v0.26):** Zerfallsraten und Quantenausbeute von Fluorophoren vor Nanostrukturen; das Quenching verlangt
  am Fußpunkt Elemente von etwa d/8, dafür ein konform verdichtetes Kugelnetz. Seit v0.27 die Fluoreszenzverstärkung aus
  Anregung und Emission: vor einer Goldkugel bis etwa 60-fach für schwache Emitter, in der Spaltmitte eines Gold-Dimers etwa
  260-fach selbst für helle. Seit v0.28 optional Krylov-Recycling für mehrere rechte Seiten; beim Dipolproblem nur
  10–20 % Gewinn, weil dem Operator ausgeprägte langsame Eigenrichtungen fehlen. Seit v0.29 auch GCRO-DR (harmonische
  Ritz-Vektoren): mit 20 statt 139 gespeicherten Richtungen 1,30-fach nahe einer Resonanz (`results_dipole.md`).
- **ACA+ (v0.30)** als Voreinstellung der H-Matrizen: gleiche Ränge, etwa 10 % genauer, robust gegen unbemerkt zu kleinen
  Rang (`results_compression.md`).
- **Parallelisierung (v0.31):** H-Matrix-Produkt, Nahfeld, Projektionen, Dipolraten und Zweitor mit OpenMP; mit 1, 2 und 4
  Threads auf zwölf Stellen gleich. Gemessen: Zweitor mit 5 760 Dreiecken auf 6 Kernen (8 Threads) 4,7-fach schneller
  (32,4 s → 6,9 s), etwa 79 % Effizienz je Kern, paralleler Anteil etwa 95 % (`results_parallel.md`).
- **Chirale Emitter (v0.32):** magnetischer Übergangsdipol, zirkular polarisierte Lumineszenz (g_lum) vor Nanostrukturen;
  achirale Goldstrukturen verdünnen g_lum im Orientierungsmittel, einzelne Orientierungen werden verstärkt (`results_dipole.md`).
- **Optische Kräfte (v0.33):** Maxwellscher Spannungstensor auf Körperflächen, Kugeln oder Parallelflächen; Strahlungsdruck auf
  eine Goldkugel wie Mie, optische Bindung am Dimer. Seit v0.34 die Kraft auf kleine, auch chirale Teilchen in Dipolnäherung
  einschließlich der enantioselektiven Kraft, gegen die volle BEM-Rechnung auf 1–2 %. Seit v0.35 Kraftkarten im Spaltfeld:
  Die chirale Kraft folgt dem Gradienten der optischen Chiralität mit dem aus der Theorie erwarteten Faktor. Seit v0.36 Kräfte in
  chiralen Medien: Spannungstensor und unabhängige Impulsbilanz im Fernfeld; helizitätsabhängiger Strahlungsdruck einer Goldkugel
  in chiraler Lösung (`results_force.md`).
- **Strahlanregung und optische Pinzette (v0.37):** allgemeine einfallende Felder (ebene Welle, Dipol, Gaußstrahl, fokussierter
  Strahl nach Richards–Wolf als exaktes Winkelspektrum); Kräfte bei Dipolanregung mit Impulserhaltung auf 10⁻⁴; Kennlinien einer
  optischen Pinzette (`tweezers`, `results_beams.md`). Seit v0.38 Strahlen in chiralen Medien (Winkelspektrum nach Helizitäten
  mit k±): Pinzette in chiraler Lösung und für chirale Teilchen; die Enantiomere werden sehr verschieden gefangen. Seit v0.39
  GLMT-Referenz (`tools/glmt.py`): Die Pinzette stimmt bei n = 8 auf 1,5 % mit der verallgemeinerten Lorenz-Mie-Theorie überein,
  seit v0.40 auch für chirale Kugeln und chirale Lösungen (Unterschied der Enantiomere auf 1,7 %), seit v0.41 für geschichtete
  Kugeln; die Pinzette rechnet beschichtete Teilchen über das Zweitor (`--shell`).
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
                         ThinLayerScatteringProblem (Dünnschicht-Näherung), TwoPortLayerProblem (Schicht als Zweitor)
apps/  tests/  tools/  results/  examples/
data/materials/          Johnson-Christy Au, Ag (refractiveindex.info, CC0)
prototype/               Python-Prototypen: ap1 (Theorie), ap2 (3D-Galerkin, H-Matrix), resonance (Resonanzfenster)
docs/                    Architektur, Ergebnisberichte, Arbeitspapiere AP 1–3, Zusammenfassung, Python-Anbindung
bindings/python/         pybind11-Modul cliffordbem._cbem, Paket cliffordbem, Tests; examples/python: Skripte
```

## Offene Punkte

Siehe `docs/ARCHITECTURE.md`. Wichtigste: reflexionsfreier Abschluss an 3D-Kanten im Resonanzfenster,
Streckung um Spitzen im C++-Kern, Krylov-Recycling für viele rechte Seiten, Substrate, für gekrümmte Elemente
Nahfeld und Kräfte sowie eine schnellere H-Matrix mit sieben Komponenten (`docs/results_curved.md`), parallele Tests auf Mehrkernrechnern, Typ-Stubs für die Python-Anbindung; für beschichtete Körper:
Block-/HODLR-Vorkonditionierung, ein Eindeutigkeitsbeweis für verschachtelte Gebiete; für die Dünnschicht-Näherung
Stabilität für d ≳ h (nichtlokale Formulierung) und Krümmungssprünge (Übergang von ebenen Seiten zu Rundungen).

Die Beweise und Aussagen in den Arbeitspapieren sind vorläufig und nicht unabhängig geprüft.

## Lizenz

MIT, siehe `LICENSE`. Materialdaten: refractiveindex.info-Datenbank (CC0).
