# Ergebnisse: Mehrkörperprobleme und Gmsh-Geometrien (C++-Kern v0.7)

## Formulierung

Für getrennte Körper B₁, …, B_m wirkt der äußere Cauchy-Operator E₂ auf der Vereinigung aller Ränder,
der innere blockdiagonal, E₁ = diag(E₁⁽ᵇ⁾), mit eigenem Medium (auch chiral) und eigener Wellenzahl je
Körper (`BlockDiagonalOperator`). J wird dreiecksweise mit dem Medium des jeweiligen Körpers gebildet.
Die physikalische Lösung erfüllt T₁h = h_inc unverändert. Einstieg für Anwendungen ist
`ScatteringProblem` (baut alle Operatoren, löst für ebene Wellen).

Für gleiche achirale Medien gibt es eine zweite, unabhängige Formulierung: ein einziger innerer
Cauchy-Operator auf der Vereinigung (Innenraum = disjunkte Vereinigung). Beide Formulierungen haben
denselben Hardy-Raum der Innenspuren, sind aber als Operatoren verschieden.

## Kontrollen (Kugel-Dimer, Radius 1, Einfall ⟂ Dimerachse, Polarisation entlang der Achse)

- **Großer Abstand:** D = 40: σ_ext = 2 σ_einzeln bis auf 2·10⁻⁴ (Wechselwirkung ∼ 1/(kD)).
- **Vertauschung verschiedener Medien** (Gold/Glas gegen Glas/Gold, spiegelsymmetrisches Netz): gleich auf 1,5·10⁻⁸.
- **Zwei Formulierungen**, Gold-Dimer, D = 3 (Spalt 1), ωa = 0,5:

| Dreiecke | blockdiagonal | vereinigt | rel. Differenz |
|---:|---:|---:|---:|
| 1 440 | 2,640471 | 2,641348 | 3,3·10⁻⁴ |
| 2 560 | 2,651207 | 2,651710 | 1,9·10⁻⁴ |
| 4 000 | 2,656496 | 2,656820 | 1,2·10⁻⁴ |

Beide Folgen konvergieren gemeinsam (Q_ext/π; die Differenz fällt mit dem Diskretisierungsfehler).
Zum Vergleich: eine einzelne Goldkugel hat Q_ext = 0,59; der Dimer mit Spalt 1 hat 2,66/2 = 1,33 je
Kugel, eine starke Verstärkung durch die Kopplung der Plasmonen.

- **Chirale Kugel neben achiraler Kugel** (ε = 2,25, χ = 0,2 bzw. 0; ωa = 1; 1 440 Dreiecke), CD in Einheiten von π:

| D | 2,5 | 3 | 4 | 6 | ∞ (einzelne chirale Kugel) |
|---|---:|---:|---:|---:|---:|
| CD | 0,275 | 0,231 | 0,188 | 0,215 | 0,225 |

## Gmsh-Import

`read_gmsh` liest ASCII-Netze im Format 2.2 und 4.1; jede physikalische Flächengruppe ist ein Körper,
die Orientierung wird über das Vorzeichen des Volumens nach außen gelegt (`test_gmsh`: Rundreise,
Format 4.1, falsch orientierte Eingabe). Geometrien erzeugt `tools/make_geometries.py` (Gmsh-Python-API):

```bash
pip install gmsh      # benötigt libGLU
python3 tools/make_geometries.py sphere 0.13 examples/sphere.msh
python3 tools/make_geometries.py bornkuhn 0.12 examples/bornkuhn.msh --angle 60
./build/scatter_mesh --mesh examples/bornkuhn.msh --omega 0.5 --media "-11,1.2" --pol circ
```

**Kugel aus Gmsh** (unstrukturiertes Netz, Gold, ωa = 0,5, Mie σ_ext = 1,853597):

| h | Dreiecke | σ_ext | Abw. |
|---:|---:|---:|---:|
| 0,25 | 540 | 1,938394 | 4,6 % |
| 0,18 | 1 014 | 1,896895 | 2,3 % |
| 0,13 | 1 948 | 1,876005 | 1,2 % |

Konvergenz zweiter Ordnung in h, wie beim Ikosaedernetz.

**Born-Kuhn-Dimer** (zwei Gold-Stäbe, Länge 3, Radius 0,3, Abstand 1,2 entlang der Einfallsrichtung, um den
Winkel φ gegeneinander gedreht; ωa = 0,5, etwa 1 800 Dreiecke): Zirkulardichroismus allein aus der
Geometrie, bei achiralem Material.

| φ | σ₊ | σ₋ | CD = σ₊ − σ₋ |
|---:|---:|---:|---:|
| +60° | 1,497574 | 1,942868 | −0,445294 |
| −60° | 1,943591 | 1,498104 | +0,445486 |
| 0° (achiral) | 2,435736 | 2,435853 | −0,000117 |

Die Enantiomere haben entgegengesetzten CD (Unterschied 0,04 %, verschiedene Netze), die spiegelsymmetrische
Anordnung hat CD ≈ 0 (Restwert auf dem Niveau der Netzasymmetrie).
