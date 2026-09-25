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
