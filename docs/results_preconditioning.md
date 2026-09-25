# Ergebnisse: Blockvorkonditionierung für mehrere Körper und chirale Medien (C++-Kern v0.9)

## Umsetzung

- `BlockPreconditioner` erhält den inneren Block über eine Funktion B ↦ E₁|_{B×B}. `ScatteringProblem::inner_block`
  liefert ihn für beliebige Dreiecksgruppen: blockdiagonal über die Körper (Paare aus verschiedenen Körpern
  tragen nichts bei), für chirale Körper P₊ E_{k₊} + P₋ E_{k₋}. Der äußere Block kommt aus den Einträgen auf der
  Vereinigung. T_BB = ½(I + E₂) + ½(J − E₁J) wird dicht per LU invertiert, außerhalb der Gruppen 2(1+J)⁻¹.
- Gruppen: `group_by_features` (Kanten/Ecken, `FeatureSet::cube(a, Mittelpunkt)`, mehrere Merkmalsmengen mit
  `append`) oder `group_by_clusters` (Blätter eines Clusterbaums mit höchstens G Dreiecken, für glatte
  Geometrien).
- Apps: `--precond point|cluster:G` in `spectrum`, `scatter_mesh`, `scatter_multi`; `scatter_cube_dimer`.

## Korrektheit

Dimer aus plasmonischer und chiraler Kugel (ε = −4 + 0,3i bzw. ε = 2,25, χ = 0,15; schräger Einfall,
zirkulare Polarisation; `test_block_precond`):

| Vorkonditionierung | GMRES | σ_ext |
|---|---:|---:|
| punktweise | 55 | 22,6665530590 |
| Cluster (32 Dreiecke) | 34 | 22,6665530588 |
| ein Block über alles | **1** | 22,6665530588 |

Ein einziger Block über alle Dreiecke ist die exakte Inverse von T₁ bis auf die H-Matrix-Toleranz; GMRES
konvergiert in einer Iteration. Das prüft den inneren Block für mehrere Körper und chirale Medien vollständig.

## Würfel-Dimer mit Kanten und Ecken

Zwei zu den Kanten gradierte Goldwürfel (Kantenlänge 2, Spalt 0,5), ωa = 0,5, Einfall ⟂ Dimerachse,
Polarisation entlang der Achse; Blöcke im Radius 0,125 um die 16 Ecken und 24 Kanten beider Würfel.

| L | Unbekannte | punktweise | Kanten/Ecken | größter Block | Blockanteil |
|---:|---:|---:|---:|---:|---:|
| 2 | 6 912 | 56 | 49 | 96 | 0,33 |
| 3 | 12 288 | 62 | 52 | 192 | 0,44 |
| 4 | 19 200 | 68 | 53 | 416 | 0,64 |
| 5 | 27 648 | 72 | 54 | 608 | 0,75 |

Punktweise wachsen die Iterationen mit der Kantenauflösung (+5 je Stufe), mit Kanten-/Eckblöcken bleiben sie
nahezu konstant (+1 je Stufe). Die verbleibenden rund 50 Iterationen stammen von der plasmonischen Kopplung
über den Spalt, nicht von den Kanten. Die Lösungen stimmen auf 6–7 Stellen überein.

## Glatte plasmonische Körper (Born-Kuhn-Dimer)

Zwei Goldstäbe (Kapseln, 1 148 Dreiecke), λ = 725 nm (Resonanz), zirkulare Polarisation:

| Vorkonditionierung | GMRES | Zeit gesamt (Aufbau + 2 Lösungen) |
|---|---:|---:|
| punktweise | 68 | 13 s |
| Cluster 64 | 54 | 13 s |
| Cluster 128 | 46 | 14 s |
| Cluster 256 | 40 | 20 s |

Bei glatten Körpern gibt es keine lokale Ursache für die Iterationszahl: Sie kommt von resonanten Moden, die
den ganzen Stab umfassen. Lokale Blöcke sparen bis 40 % der Iterationen, aber Aufbau und Anwendung der Blöcke
kosten ebenso viel. Wirksam wäre hier eine Vorkonditionierung, die ganze Körper erfasst (H-LU je Körper, AP 3.4)
oder Deflation der resonanten Moden.

## Hierarchische Faktorisierung (v0.10)

`HodlrSolver` zerlegt das Gesamtsystem T₁ (alle Körper, auch chiral) hierarchisch im HODLR-Format: Auf jeder
Stufe des geometrischen Clusterbaums werden die beiden Nebendiagonalblöcke mit einer Block-ACA (8×8-Pivots,
Pseudoinverse des Pivotblocks, Nachkompression) niedrigrangig approximiert, die Blätter dicht per LU zerlegt,
und die Inverse über die Woodbury-Formel angewandt. Die Einträge liefert `ScatteringProblem::system_entry`
(T_ij als 8×8-Block aus den gespeicherten Nahfeld- bzw. Fernfeldregeln; gegen `T.apply` auf 6,5·10⁻¹³ geprüft).
Das ist eine hierarchische LU mit **schwacher Zulässigkeit**; die klassische H-LU mit starker Zulässigkeit
(H-Matrix-Arithmetik mit gekürzten Produkten) ist nicht umgesetzt.

**Korrektheit** (`test_hodlr`, Dimer aus plasmonischer und chiraler Kugel, 2 880 Unbekannte):

| Toleranz | GMRES (punktweise: 63) | max. Rang | Speicher | σ_ext |
|---:|---:|---:|---:|---:|
| 10⁻⁹ | 2 | 617 | 334 MB | 26,3047285276 |
| 10⁻² | 19 | 131 | 47 MB | 26,3047285279 |

Mit feiner Toleranz ist die Faktorisierung ein direkter Löser.

**Würfel-Dimer** (zwei gradierte Goldwürfel, L = 3, 12 288 Unbekannte, Spalt 0,5; eine rechte Seite):

| Vorkonditionierung | GMRES | Lösen | Aufbau | Rang | Speicher |
|---|---:|---:|---:|---:|---:|
| punktweise | 62 | 4,7 s | – | – | – |
| Kanten/Ecken (v0.9) | 52 | – | < 1 s | – | – |
| HODLR 10⁻¹ | 30 | 3,2 s | 10,8 s | 100 | 156 MB |
| HODLR 10⁻² | 19 | 2,6 s | 37,8 s | 182 | 352 MB |

**Born-Kuhn-Dimer an der Resonanz** (λ = 725 nm, 9 184 Unbekannte, zwei rechte Seiten):

| Vorkonditionierung | GMRES | Aufbau | Rang | Gesamtzeit |
|---|---:|---:|---:|---:|
| punktweise | 68 | – | – | 12 s |
| HODLR 10⁻¹ | 51 | 3,4 s | 38 | 15 s |
| HODLR 10⁻² | 36 | 13,7 s | 111 | 24 s |
| HODLR 10⁻³ | 26 | 60,7 s | 272 | 71 s |

**Bewertung.**
- Die Faktorisierung lohnt sich nur bei vielen rechten Seiten für denselben Operator. Beim Würfel-Dimer liegt
  der Break-even bei etwa 7 (Toleranz 10⁻¹) bzw. 18 (10⁻²) rechten Seiten; eine Orientierungsmittelung mit 26
  Richtungen und zwei Polarisationen hat 52.
- An der Plasmonresonanz ist sie schwach: T₁ hat dort sehr kleine Singulärwerte, und der Fehler der
  Näherungsinverse wird mit der Kondition verstärkt. Erst Toleranzen um 10⁻³ senken die Iterationen deutlich,
  dann dominieren Aufbau und Speicher.
- Die Ränge mit schwacher Zulässigkeit wachsen, weil benachbarte Cluster auf derselben Fläche direkt aneinander
  grenzen (bis 272 bei 10⁻³). Eine H-LU mit starker Zulässigkeit hätte kleinere Ränge, würde am
  Konditionsproblem an der Resonanz aber nichts ändern.
- Für Resonanzen mit vielen rechten Seiten sind Krylov-Recycling über die rechten Seiten (z. B. GCRO-DR) oder
  Deflation der resonanten Moden die geeigneteren Werkzeuge. Der Aufbau der Faktorisierung ist noch nicht
  parallelisiert.
