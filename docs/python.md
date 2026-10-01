# Python-Anbindung (v0.42)

Das Paket `cliffordbem` macht den C++-Kern für Skripte zugänglich: Parameterstudien, eigene Auswertungen, Abbildungen und
eigene Formulierungen aus den Bausteinen des Kerns, ohne für jede Fragestellung ein neues C++-Programm zu schreiben. Die
Rechnung bleibt vollständig in C++ (H-Matrizen, Quadratur, GMRES, OpenMP); Python übergibt Eingaben und erhält NumPy-Arrays.

## Bauen und installieren

```bash
# Installation (baut Kern und Modul mit CMake über scikit-build-core; benötigt einen C++17-Compiler)
pip install .

# Entwicklung ohne Installation
pip install pybind11 numpy
cmake -S . -B build -DCBEM_BUILD_PYTHON=ON
cmake --build build -j
PYTHONPATH=build/python python3 -c "import cliffordbem as cb; print(cb.__version__)"
ctest --test-dir build -R test_python --output-on-failure      # Python-Tests (benötigen scipy für die Mie-Referenzen)
```

`CBEM_BUILD_PYTHON` ist standardmäßig aus; der Kern bleibt ohne externe Abhängigkeiten. Mit der Option wird `cbem` mit
positionsunabhängigem Code übersetzt und das Modul `cliffordbem._cbem` nach `build/python/cliffordbem` gelegt. pybind11 wird
als CMake-Paket gesucht, sonst über `python -m pybind11 --cmakedir`. Die Version kommt aus `CMakeLists.txt` (eine Quelle für
CMake, Modul und Paket). Die Materialtabellen werden bei der Installation ins Paket kopiert; `cb.DATA_DIR` zeigt auf sie
(Vorrang: Umgebungsvariable `CBEM_DATA_DIR`, Paket, Quellbaum).

## Schnellstart

```python
import numpy as np
import cliffordbem as cb

gold = cb.make_material("Au")                        # Johnson-Christy
unit, lam = 20.0, 530.0                              # Längeneinheit 20 nm (Kugelradius 1), Wellenlänge in nm
om = cb.omega_from_wavelength(lam, unit)             # omega = k0 = 2 pi unit / lambda
water = cb.Medium(eps=1.33**2)

P = cb.ScatteringProblem(cb.make_icosphere(8), gold.medium(lam), om, outer=water)
r = P.solve_plane_wave(d=(0, 0, 1), p=(1, 0, 0))
print(r.sigma_ext * unit**2, "nm^2", r.iterations, "Iterationen")

# Nahfeld, Kraft
pts = np.array([[1.5, 0, 0], [0, 0, 2.0]])
nf = cb.exterior_near_field(P.mesh, r.h, water, om, (0, 0, 1), (1, 0, 0), pts)
F = cb.force_on_sphere(P.mesh, r.h, water, om, (0, 0, 1), (1, 0, 0), center=(0, 0, 0), R=1.5)

# Spektrum wie die App spectrum (Orientierungsmittelung, zirkular: CD)
sp = cb.spectrum(cb.make_icosphere(8), gold, np.arange(450, 651, 10), unit=20, n_bg=1.33)
```

## Konventionen

| Größe | Python |
|---|---|
| Einheiten | wie im Kern: Längen in einer frei gewählten Einheit L, ε₀ = μ₀ = 1, ω = k₀ = 2πL/λ, Zeitkonvention e^{−iωt}, σ in L² |
| `Vec3`, `CVec3` | jede Folge von drei Zahlen (Tupel, Liste, `ndarray`); Rückgabe als `ndarray` (float64 bzw. complex128) der Länge 3 |
| Spuren | flaches complex128-Array der Länge 8N wie im Kern; `cb.trace_blocks(h)` gibt (N, 8), Blade-Reihenfolge 1, e1, e2, e12, e3, e13, e23, e123 |
| Streuspur | `h - b` mit `b = cb.plane_wave_trace(mesh, outer, omega, d, p)` (dieselbe rechte Seite wie `solve_plane_wave`, auch im chiralen Medium) |
| Punkte | float64-Array (M, 3); ein einzelner Punkt als (3,) |
| Nahfeld | `dict` mit `E`, `H` (M, 3), `inside`, `too_close` (M,), `enhancement` = \|E\|²/\|E₀\|², `chirality` (M,) |
| 8×8-Blöcke | `ndarray` (8, 8), Zeile r, Spalte c (`left_matrix`, `transmission_map`, `system_entry`, `TransmissionOperator.J`: (N, 8, 8)) |
| Parameter | `HMatrixParams`, `EntryParams`, `SolveOptions`, `HodlrParams`, `NearFieldOptions`, `TwoPortOptions` mit Schlüsselwortargumenten, z. B. `cb.HMatrixParams(eps=1e-6)`; unbekannte Namen → `AttributeError` |
| Fehler | `std::runtime_error` → `RuntimeError`, `std::invalid_argument` → `ValueError`; falsche Längen und Formen werden vor der Rechnung als `ValueError` gemeldet |

Netze sind in Python unveränderlich: `points`, `triangles`, `normals`, `centroids`, `areas`, `hmax` liefern Kopien; ein
geändertes Netz wird neu erzeugt (`cb.TriangleMesh(points, triangles)`, `cb.translated`, `cb.offset_surface`). Netze,
Medien und Multivektoren lassen sich kopieren und mit `pickle` speichern.

## Umfang

| Bereich | Klassen und Funktionen |
|---|---|
| Algebra | `Multivector` (`*` geometrisches Produkt, `+`, `-`, Skalare, `reverse`, `involute`, `grade`, `inverse`, `left_matrix`, Zugriff `M["e12"]`), Konstanten `e1, e2, e3, e12, e13, e23, I`, `helicity_projector`, `transmission_map` |
| Medien | `Medium(eps, mu, chi)` mit `k(omega, helicity)`, `Material`, `ConstantMaterial`, `TabulatedMaterial` (auch `from_yaml`), `make_material("Au" \| "Ag" \| "re,im" \| Pfad.yml)` |
| Geometrie | `TriangleMesh`, `make_icosphere`, `make_icosphere_graded`, `make_cube_uniform`, `make_cube_graded`, `translated`, `offset_surface`, `distance_to_surface`, `winding_number`, `signed_volume`, `require_separated`, `read_gmsh`, `write_gmsh22`, `lebedev` |
| Probleme | `ScatteringProblem` (ein oder mehrere Körper, chiral, Recycling, GCRO-DR, Block- und HODLR-Vorkonditionierung, `solve_rhs`), `LayeredGeometry` + `LayeredScatteringProblem`, `ThinBody` + `ThinLayerScatteringProblem`, `TwoPortBody` + `TwoPortLayerProblem`, Ergebnisse `PlaneWaveResult`, `LayeredResult` |
| Anregung, Fernfeld | `plane_wave_trace`, `project_plane_wave`, `circular_polarization`, `far_field`, `far_field_E`, `extinction_cross_section`, `forward_amplitude`, `plane_wave_incidence`, `helicity_part`, `extinction_in_medium`, `forward_amplitude_in_medium`, `PlaneWaveField`, `ZeroField`, `DipoleField`, `BeamField.focused`, `BeamField.gaussian` |
| Nahfeld | `exterior_near_field` (ebene Welle oder beliebige Anregung), `scattered_field`, `NearFieldOperator`, `plane_wave_evaluator`, `near_field_evaluator` |
| Dipole | `dipole_field`, `project_dipole`, `dipole_rates` → `DipoleRates`, `fluorescence_enhancement`, `helicity_projector_sign` |
| Kräfte | `force_from_traces`, `force_on_sphere`, `force_on_offset`, `force_from_far_field`, `emitter_force`, `radiated_momentum`, `stress_dot_normal`, `fields_with_gradients` → `FieldGradient`, `DipolePolarizability`, `polarizability_from_mie`, `dipole_particle_force` |
| Bausteine | `KernelEntries`, `KernelHMatrix` (`stats` → `HStats`), `CauchyOperator`, `ChiralCauchyOperator`, `TransmissionOperator`, `gmres` (C++-Operator oder Python-Funktion, Vorkonditionierer als Python-Funktion), `group_by_clusters`, `group_by_features`, `FeatureSet` |
| Komfort | `spectrum` (wie die App, homogene Körper), `omega_from_wavelength`, `wavelength_from_omega`, `polarization_basis`, `trace_blocks`, `set_num_threads`, `omp_threads` |

Die Docstrings (`help(cb.ScatteringProblem)`) beschreiben Argumente und Rückgaben; die Bedeutung der Größen ist in den
C++-Headern unter `include/cbem` erklärt, deren Namen die Anbindung übernimmt.

## Lebensdauer und Threads

- **Problemklassen kopieren ihre Netze**; Python-Objekte dürfen nach dem Aufbau freigegeben werden.
- **Bausteine mit Referenzen** (`KernelEntries` → Netz, `KernelHMatrix` → Einträge, `CauchyOperator` → H-Matrix,
  `TransmissionOperator` → Operatoren, `NearFieldOperator` → Netz) halten ihre Eingaben über `keep_alive` am Leben; die
  Kette `CauchyOperator(m, KernelHMatrix(KernelEntries(m, k)))` ist ohne Zwischenvariablen sicher.
- **Feldauswerter** (`plane_wave_evaluator`, `near_field_evaluator`) halten eigene Kopien von Netz und Spur; im Kern
  hält `make_near_field_eval` beide nur per Referenz.
- **Vektorfelder** von `LayeredGeometry`, `ThinBody` und `TwoPortBody` werden als Kopien geliefert (Referenzen in die
  Vektoren wären nach `add_surface` ungültig).
- **Der GIL wird freigegeben**, während der Kern rechnet (Aufbau, Lösen, Nahfeld, Kräfte, Projektionen). Python-Threads
  können daher gleichzeitig rechnen und blockieren sich nicht; die OpenMP-Parallelisierung des Kerns bleibt unverändert
  (`cb.set_num_threads`). Python-Funktionen, die an `gmres` übergeben werden, holen sich den GIL bei jedem Aufruf.
- Eine laufende Rechnung lässt sich nicht mit Strg-C abbrechen; der Abbruch wirkt nach ihrem Ende.

## Beispiele

| Skript | Inhalt |
|---|---|
| `examples/python/spectrum_gold_sphere.py` | Spektrum einer Goldkugel in Wasser gegen Mie, CSV und Abbildung |
| `examples/python/nearfield_dimer.py` | Nahfeldkarte eines Gold-Dimers, Kräfte je Kugel gegen die Impulsbilanz im Fernfeld |
| `examples/python/custom_formulation.py` | T₁ in NumPy aus den Cauchy-Operatoren des Kerns, Plemelj-Spurtrennung, GMRES mit Python-Vorkonditionierer |

Mit n = 4 (Dimer: 640 Dreiecke, 41 × 21 Punkte) laufen die Beispiele in wenigen Sekunden; das Dimer ergibt Kräfte
±383,49 entlang der Achse (Anziehung) und in Ausbreitungsrichtung 29,740 gegen 29,756 aus der Impulsbilanz (0,05 %).
Die von Hand zusammengesetzte T₁ stimmt mit dem Operator des Kerns auf allen Stellen überein.

## Tests (`bindings/python/tests/test_python.py`, ctest `test_python`)

| Test | prüft |
|---|---|
| `test_algebra` | e1 e2 = e12, I² = −1, Zentralität von I, Linksmultiplikation, Inverse, Projektor als Nullteiler, pickle |
| `test_mesh` | Netzgrößen, Orientierung, Rundreise über Arrays und pickle, Parallelfläche, Windungszahl, Abstand, Gmsh-Rundreise, Fehlermeldungen |
| `test_materials_and_params` | Johnson-Christy, Tabellen, chirale Wellenzahlen, Parameter mit Schlüsselwortargumenten |
| `test_scattering_regression_and_low_level` | Glaskugel wie `test_scattering` (0,20291); derselbe Weg aus Bausteinen bitgleich; optisches Theorem, S(0), `solve_rhs` |
| `test_gold_against_mie` | Goldkugel näher an Mie als die Prototyp-Regel; Kugel aus dem Außenmedium unsichtbar |
| `test_chiral` | Unsichtbarkeit im chiralen Medium, Ablehnung linearer Polarisation, Spiegelsymmetrie σ_s(χ) = σ_{−s}(−χ) |
| `test_multibody_and_spectrum` | weit getrennte Körper additiv, `spectrum` gleich der direkten Rechnung, kein CD der Kugel |
| `test_near_field_against_mie` | Nahfeld gegen Mie (Fehler fällt mit O(h²)), Markierung innerer Punkte, Auswerter, H-Matrix = direkte Summation |
| `test_forces` | Randspuren, Kugel und Fernfeld auf 1 % gleich, gegen Mie, Gradienten und Dipolkraft |
| `test_dipole_and_beams` | Raten ohne Streuer = 1, Dipolfeld, Gaußstrahl mit Leistung 1 und Poynting-Fluss 1 |
| `test_layered_neutral` | neutrale Schicht: Dünnschicht exakt ohne Wirkung, Zweitor und exakte Rechnung bis auf die Diskretisierung |
| `test_errors_and_gmres` | Längen- und Formprüfungen, GMRES mit Python-Operator, Ausnahmen aus Python-Rückrufen |
| `test_threads_release_gil` | zwei Lösungen parallel in Python-Threads, gleiche Ergebnisse, GIL frei während der Rechnung |

Laufzeit etwa 40 s auf einem Kern.

## Grenzen

- Eigene einfallende Felder und Materialmodelle lassen sich nicht als Python-Unterklassen an den Kern übergeben (die
  Auswertung je Quadraturpunkt würde den GIL brauchen); eigene rechte Seiten gehen über `solve_rhs` mit einer selbst
  berechneten Spur, eigene Materialien über `Medium(eps=...)` je Wellenlänge.
- Typ-Stubs (`.pyi`) fehlen noch; Signaturen zeigen `help()` und die Docstrings.
- Gebaut und getestet mit g++ 13, Python 3.12, pybind11 3.1 unter Linux; andere Plattformen sind nicht geprüft.
