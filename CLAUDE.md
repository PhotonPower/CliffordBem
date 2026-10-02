# CLAUDE.md – Arbeitsanleitung für CliffordBem

Koordinatenfreie Galerkin-Randelementmethode (BEM) der Nano-Optik in der komplexifizierten Clifford-Algebra Cl₃(ℂ):
Streuung, Nahfeld, Dipolemission und optische Kräfte an (auch chiralen, beschichteten, mehreren) Körpern. C++17-Kern ohne
externe Abhängigkeiten, optionale Python-Anbindung (pybind11). Sprache des Projekts: **Deutsch**.

**Stand: v0.59** (siehe `CHANGELOG.md`, neueste Version oben). Drei Diskretisierungen stehen nebeneinander:

| Dichten / Geometrie | Klasse | seit | Genauigkeit an der Kugel (320 Elemente) |
|---|---|---|---|
| stückweise konstant / eben | `ScatteringProblem` (+ geschichtete Varianten) | v0.1 | Glas −5,7 %, Gold +7,3 % gegen Mie, O(h²) |
| unstetig linear / eben | `LinearScatteringProblem` | v0.45 | kein Gewinn an glatten Körpern (Geometrie dominiert) |
| unstetig linear / quadratisch gekrümmt | `CurvedScatteringProblem` | v0.47 | Glas −0,012 %, Gold −0,007 %, etwa O(h⁴) |

Die Herleitung, Messungen und Begründungen der gekrümmten Elemente stehen in `docs/results_curved.md` (Stufen 1, 2a, 1b,
2b, Gmsh, Nahquadratur, H-Matrix, Geometrie, Gauß-Regeln, Sauter-Schwab, Produkt, Aufbau, Nahfeld und Kräfte, chirales Außenmedium) – vor Arbeiten an linearen/gekrümmten Elementen lesen.

## Bauen und testen

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCBEM_BUILD_PYTHON=ON   # Python optional (pip install pybind11 numpy scipy)
cmake --build build -j
ctest --test-dir build --output-on-failure                              # 30 Tests mit Python, 29 ohne
ctest --test-dir build -R "test_curved|test_linear" --output-on-failure # Teilmenge
PYTHONPATH=build/python python3 -c "import cliffordbem as cb; print(cb.__version__)"
pip install .                                                           # Python-Paket über scikit-build-core
```

- Optionen: `CBEM_OPENMP` (an), `CBEM_NATIVE` (aus; für eigene Rechnungen lohnend, Produkt doppelt so schnell), `CBEM_FAST_COMPLEX` (an, `-fcx-fortran-rules`), `CBEM_BUILD_TESTS` (an), `CBEM_BUILD_APPS` (an), `CBEM_BUILD_PYTHON` (aus).
- **Laufzeit:** Der volle Testlauf dauert auf einem Kern etwa 45 Minuten. Die längsten Tests sind `test_optical_force`
  (~200 s), `test_hodlr`, `test_twoport`, `test_dipole` (~170 s), `test_curved` (~90 s) und `test_python` (~85 s).
  Bei Änderungen nur die betroffenen Tests laufen lassen; vor Versionen, die gemeinsam genutzten Code ändern, den vollen
  Lauf. Lange Läufe im Hintergrund mit Protokolldatei starten.
- Python-Tests: `bindings/python/tests/test_python.py` (eigener Runner, pytest-kompatibel, `test_*`-Funktionen), braucht
  `CBEM_TOOLS_DIR=tools` für die Mie-Referenzen; in ctest als `test_python`, immer mit `OMP_NUM_THREADS=1` (vergleicht Wege
  bitgenau, die H-Matrix summiert mit mehreren Threads in nicht festgelegter Reihenfolge).
- Windows: MinGW-g++ aus MSYS2 ucrt64; Python-Anbindung nur gegen das MSYS2-Python (`pacman -S
  mingw-w64-ucrt-x86_64-pybind11 mingw-w64-ucrt-x86_64-python-numpy mingw-w64-ucrt-x86_64-python-scipy`, dann
  `-DPython_EXECUTABLE=C:/msys64/ucrt64/bin/python.exe`), nicht gegen das python.org-Python (MSVC-ABI).
- Optional: `pip install gmsh` (benötigt `libglu1-mesa libxcursor1 libxinerama1 libxft2`) für `examples/python/gmsh_curved.py`.

## Aufbau

```
include/cbem/, src/   C++-Kern; Unterordner nach Modul (siehe unten)
apps/                 Kommandozeilenprogramme (spectrum, tweezers, scatter_*, nearfield, dipole, bench_compression)
tests/                C++-Tests (check.hpp: CHECK(cond, fmt, ...), REPORT()); in CMakeLists.txt in der foreach-Liste
bindings/python/      pybind11-Modul cliffordbem._cbem (bind_*.cpp), Paket cliffordbem/ (__init__.py, fields.py), Tests
tools/                Python-Referenzen: mie*.py (Mie, chiral, beschichtet, Nahfeld, Kraft, Dipol), glmt.py, Auswertungen
prototype/            Python-Prototypen (ap1, ap2, resonance) und Messungen zu gekrümmten Elementen (curved/)
docs/                 ARCHITECTURE.md (Module, Ausbauplan), results_*.md (Ergebnisberichte), python.md, README.md (Index)
results/              Rohdaten (CSV) der Ergebnisberichte
data/                 materials/ (Johnson-Christy Au, Ag), meshes/ (echte Gmsh-Datei zweiter Ordnung als Testdaten)
examples/             README.md, python/ (Skripte)
```

Module (`include/cbem/…`): `core` (Typen, Materialien), `clifford` (Multivector), `geometry` (Netze, quadratische
Elemente `quadratic_mesh`, Gmsh, Quadratur, Sauter-Schwab), `kernel` (Dirac-Kern, analytische Dreiecksintegrale auch für
lineare Gewichte), `assembly` (`KernelEntries`, `LinearKernelEntries`, `CurvedKernelEntries`), `hmatrix` (Clusterbaum,
ACA, `KernelHMatrix`, `CurvedHMatrix`), `operators` (Cauchy-, Transmissions-, chirale, Mehrkörper-, gekrümmte Operatoren),
`solvers` (GMRES, Recycling, Block-/HODLR-Vorkonditionierer), `problems` (Streuprobleme), `sources` (ebene Wellen,
Strahlen, Dipole, Nahfeld, Kräfte, chirale Anregung).

## Formulierung und Datenkonventionen

- Einheiten: Länge in frei gewählter Einheit L, ε₀ = μ₀ = 1, ω = k₀ = 2πL/λ, Zeit e^{−iωt}, Pasteur-Medien
  D = εE + iχH, B = μH − iχE. Spur F = √ε E + I √μ H. Gleichung: T₁ = ½(1 + E₂) + ½(1 − E₁)J.
- Blades als Bitmasken, Reihenfolge der 8 Koeffizienten: 1, e1, e2, e12, e3, e13, e23, e123.
- Spuren: konstant 8N Koeffizienten (Basis 1/√A); linear und gekrümmt 24N, Index `8·(3t + a) + q` mit je Element
  orthonormierter Basis ψ_a = Σ S_ak λ_k (eben: `psi_matrix(A)` analytisch; gekrümmt: `curved_psi_matrices`, Cholesky).
  Die Massenmatrix ist damit die Identität, und T₁ behält seine Form.
- Kernkomponenten: konstant/linear 4 (Skalar + Vektor von Φ_k, die Normale steht außerhalb des Integrals); gekrümmt 7
  (G(z)n = s n + v z·n + v z∧n, Blades 1, e1, e2, e3, e12, e13, e23), weil die Normale im Element variiert.
- Selbstterm: Der singuläre Anteil Φ₀ = z/(4πr³) verschwindet nur für konstante Dichten auf ebenen Elementen. Für
  lineare Dichten wird er antisymmetrisiert integriert; auf gekrümmten Elementen zusätzlich der schwach singuläre Rest
  K_s = Φ₀(n(y) − n(x))/2 (ohne ihn: Plemelj-Fehler O(h) statt O(h²) – Gegentest in `test_curved`).
- Gekrümmte Transmissionsabbildung: Galerkin-Projektion J_G (24×24 je Element) statt exakter Komposition (Prototyp: 20 %
  des Diskretisierungsfehlers, gleiche Ordnung).

## Konventionen

- **Sprache:** Code-Kommentare und Bezeichnungen in Ausgaben auf Deutsch **ohne Umlaute** (ae, oe, ue, ss);
  Markdown-Dokumentation mit Umlauten. Commit-Nachrichten ohne Umlaute.
- **Versionen:** Jedes Feature ist eine Version 0.NN. Dazu gehören:
  1. `project(CliffordBem VERSION 0.NN.0 …)` in `CMakeLists.txt` (`pyproject.toml` liest die Version von dort);
  2. eine neue Zeile **oben** in der Tabelle in `CHANGELOG.md` (Inhalt | Bericht);
  3. „Stand 0.NN“ in `README.md` und `docs/ARCHITECTURE.md` sowie die betroffene Modulzeile dort; Tabellen der Fähigkeiten
     und Tests in `README.md`; `docs/python.md` bei Python-Änderungen; der passende `docs/results_*.md`-Bericht.
- **Commits:** `v0.NN: Kurztitel` mit einem erklärenden Text (was, Ergebnisse mit Zahlen, Prüfungen). Autor wie bisher:
  `git -c user.name="Claude (Anthropic)" -c user.email="noreply@anthropic.com" commit …`. Kein Token im Repository.
- **Prüfphilosophie:** Jede Funktion bekommt einen Test gegen eine *unabhängige* Referenz (Mie, feine Quadratur, dichte
  Matrix, entarteter Grenzfall, Prototyp). Bewährt hat sich:
  - erst messen, dann bauen (Prototyp in `prototype/` mit der Python-Anbindung, z. B. Galerkin-Projektion auf feine Netze);
  - entartete Proben: neuer Pfad mit Daten des alten Pfads muss bitgleich oder auf Quadraturgenauigkeit gleich sein
    (lineare Einträge summieren sich exakt zu den konstanten; gekrümmt mit Mitten auf den Sehnen = eben);
  - Konvergenzordnung in h prüfen (ein Fehlerterm zeigt sich als falsche Ordnung oder Stagnation);
  - schlägt eine Schwelle fehl, zuerst klären, welche Seite ungenau ist (mehrfach war es die Referenz), nicht die
    Schwelle aufweichen; unerwartete Ergebnisse offen berichten.
- **H-Matrix:** Einträge seit v0.56 in einfacher Genauigkeit gespeichert (`HMatrixParams::single_precision`), gekrümmte
  Elemente mit eta = 2 und ACA-Toleranz 10⁻⁶ (`curved_hmatrix_params()`, v0.57); bitgenaue Vergleiche zwischen Bauten nur mit gleichen Schaltern.
- **Bestehende Pfade nicht verändern:** Neue Diskretisierungen kamen als zusätzliche Klassen; der konstante Pfad blieb
  Operation für Operation gleich (Regressionswerte, z. B. `test_scattering` 0,20291).
- Python-Anbindung: `keep_alive` für Objekte mit Referenzmembern, GIL freigeben (`nogil`, `call_guard`), Vektor-Member als
  Kopien zurückgeben (Referenzen werden nach Reallokation ungültig), Spuren/Punkte als NumPy-Arrays.

## Stolperfallen aus der bisherigen Arbeit

- `pkill -f <muster>` trifft auch die eigene Shell, wenn das Muster in der Befehlszeile steht; Prozesse per PID beenden.
- Lang laufende Befehle (Bau, ctest) mit `setsid nohup … &` in eine Protokolldatei starten und kurz abfragen.
- Windows/MinGW (MSYS2 ucrt64): Programme werden statisch gelinkt (`-static` in `CMakeLists.txt`). Ohne das luden sie
  über den PATH die `libstdc++-6.dll` aus Git for Windows, und `test_materials`/`test_gmsh` stürzten in den Dateistreams ab.
- `std::complex<double>`-Multiplikation in heißen Schleifen: GCC erzeugt nach jeder Multiplikation eine NaN-Prüfung mit
  Rückfall auf `__muldc3` (sichtbar mit `objdump -dr datei.obj | grep __muldc3`); das verhindert die Vektorisierung. In
  inneren Schleifen reell rechnen und in lokalen Feldern akkumulieren (v0.54: Produkt fast doppelt so schnell, bitgleich);
  seit v0.55 baut die Bibliothek mit `-fcx-fortran-rules` (`CBEM_FAST_COMPLEX`), dann gibt es keine Aufrufstellen mehr.
  Laufzeiten nie nur für den Aufbau messen: das Lösen kostete bei 1 280 Elementen ebenso viel.
- Textersetzungen in Quelldateien nur mit eindeutigem Anker und Prüfung (`assert a in s`); eine Ersetzung in
  `src/hmatrix.cpp` traf einmal zwei Stellen.
- Neue Überladungen (z. B. `translated`, `signed_volume` für `QuadraticMesh`) machen `&funktion` in den pybind11-Bindungen
  mehrdeutig: `static_cast` auf die gewünschte Signatur.
- `Medium` ist in `operators/transmission_operator.hpp` definiert, `signed_volume(TriangleMesh)` in `geometry/gmsh_io.hpp`.
- Die analytische Nahquadratur der ebenen Elemente hat selbst einen Fehlerboden (3,5·10⁻⁶ mit verschärfter Außenregel,
  1,7·10⁻⁵ in der Voreinstellung) – als Referenz für Tests verschärfen (`EntryParams::adapt_ratio = 0.1`).
- Sauter-Schwab stagniert bei gestreckten Elementen (Seitenverhältnis > 1,6) bei etwa 3·10⁻⁵; dort halbanalytisch.
- ACA erhält Symmetrien nicht exakt (Spiegelsymmetrie chiral: 4·10⁻⁶ bei ε = 10⁻⁴); Symmetrietests mit ε = 10⁻⁶.
- Gold und andere Metalle reagieren empfindlich (Resonanz, Absorption fern der Resonanz verstärkt Fehler etwa um
  |Re ε|/Im ε); kleine Quadraturfehler sind dort sichtbar.

## Offene Aufgaben (priorisiert)

1. **Aufbau gekrümmter Elemente weiter beschleunigen**: auf dem Windows-Rechner (v0.57) 23 s bei 1 280 Elementen
   (Gold, ein Kern, `prototype/curved/build_profile.cpp`): H-Matrizen ≈ 11 s (ACA auf Quadraturpunkten gemessen: −24 % auf
   den ACA-Teil, ≈ −10 % Aufbau, nicht eingebaut; Fernblöcke (t, s)/(s, t) paaren, geschätzt −8 %), getrennte Nahpaare
   6,5 s, Sauter-Schwab ≈ 4,6 s (seit v0.57 gepaart). Kein großer Einzelhebel mehr (results_curved.md).
   Sauter-Schwab auf verzerrten Elementen (längste Kante/Höhe > 1,7) konvergiert langsam (Ausreißer bis 10⁻⁴ auf
   Gmsh-Netzen) – Ordnung nach Elementform wählen? Lösen kostet so viel wie der Aufbau. Immer erst messen.
2. **Nahfeld und Kräfte auf gekrümmten Elementen (v0.58, Rest):** Kraft aus Randspuren und Fernfeld
   (`force_from_traces`, `force_from_far_field`), H-Matrix für viele Auswertepunkte (wie `NearFieldOperator`; heute
   direkte Summation), einfallende Felder aus Python (`fields(x)`) für `exterior_near_field_curved`.
3. Block-/HODLR-Vorkonditionierung und geschichtete Körper für lineare/gekrümmte Elemente.
4. Typ-Stubs (`.pyi`) für die Python-Anbindung; Wheels für weitere Plattformen.
5. Ältere offene Punkte aus `README.md`/`docs/ARCHITECTURE.md`: reflexionsfreier Abschluss an 3D-Kanten, Streckung um
   Spitzen im C++-Kern, Substrate, parallele Tests auf Mehrkernrechnern, Stabilität der Dünnschicht-Näherung für d ≳ h.

## Kurzbeispiel (Python)

```python
import cliffordbem as cb
q = cb.quadratic_icosphere(4)                       # oder cb.read_gmsh_quadratic("koerper.msh")[0][1]
P = cb.CurvedScatteringProblem(q, cb.Medium(eps=-11 + 1.2j), omega=0.5, outer=cb.Medium(eps=1.33**2))
r = P.solve_plane_wave(d=(0, 0, 1), p=(1, 0, 0))
print(r.sigma_ext / cb.pi, r.iterations)
```
