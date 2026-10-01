# Ergebnisse: Optische Kräfte (v0.33)

Kräfte auf die BEM-Körper über den Maxwellschen Spannungstensor (`include/cbem/sources/optical_force.hpp`), zeitgemittelt,
ε₀ = μ₀ = 1, achirales Außenmedium:

    ⟨T⟩ = ½ Re[ε E ⊗ E* + μ H ⊗ H* − ½(ε|E|² + μ|H|²) I],   F = ∮_S ⟨T⟩·n dS

über eine geschlossene Fläche S im Außenraum, die genau den betrachteten Körper umschließt. Drei Wege:

- **Randspuren** (`force_from_traces`): S ist die Körperfläche selbst, Felder aus den stückweise konstanten Außenspuren; ohne
  Nahfeldauswertung, je Körper.
- **Kugel** (`force_on_sphere`): Kugel um den Körper bzw. um alle Körper, Gauß-Legendre × gleichmäßig in φ, Felder aus dem Nahfeld.
- **Parallelfläche** (`force_on_offset`): Parallelfläche des Körpers im Abstand δ (bei einem Dimer kleiner als der halbe Spalt),
  7-Punkt-Regel, Felder aus dem Nahfeld.

Normierung: Eine ebene Welle der Amplitude |E₀| trägt die Impulsstromdichte ½ε|E₀|²; der Strahlungsdruck auf einen Körper ist
F_z = σ_pr·½ε|E₀|² mit σ_pr = σ_ext − ⟨cos θ⟩ σ_sca (Referenz `tools/mie_force.py`, Bohren–Huffman; verlustfreie Kugel:
σ_ext = σ_sca auf zehn Stellen).

## Strahlungsdruck auf eine Goldkugel

ε = −11 + 1,2i, Radius 1, Wasser, ωa = 0,5, Einfall entlang z, x-polarisiert; Mie: F_z = 10,6298.

| Weg | n = 8 | n = 12 |
|---|---:|---:|
| Randspuren | 10,4498 (−1,69 %) | 10,5487 (−0,76 %) |
| Kugel R = 1,5 | 10,4680 (−1,52 %) | 10,5569 (−0,69 %) |
| Kugel R = 3 | 10,4828 (−1,38 %) | 10,5635 (−0,62 %) |
| Parallelfläche δ = 0,2 | 10,4656 (−1,55 %) | 10,5558 (−0,70 %) |
| Parallelfläche δ = 0,05 | 10,4640 (−1,56 %) | 10,5550 (−0,70 %) |

Die Wege stimmen untereinander auf etwa 0,3 % überein; die Abweichung gegen Mie ist der Diskretisierungsfehler der Streulösung
selbst (σ_ext auf denselben Netzen −1,4 bzw. −0,6 %). Die Integration des Spannungstensors fügt praktisch keinen Fehler hinzu.
Querkräfte verschwinden (10⁻⁵ bis 10⁻⁸).

## Optische Bindung am Gold-Dimer

Zwei Goldkugeln (Radius 1) mit Spalt 0,3 entlang x in Wasser, Einfall entlang z, n = 8; Parallelfläche δ = 0,06:

| | Polarisation entlang der Achse | senkrecht |
|---|---:|---:|
| Kraft entlang der Achse auf Kugel 1 / 2 | +43,26 / −43,26 | −0,50 / +0,50 |
| Kraft in Einfallsrichtung je Kugel | 17,89 | 8,02 |
| Summe beider Kugeln / Kugel um beide | 35,78 / 35,81 | 16,03 / 16,07 |

- Entlang der Achse polarisiert ziehen sich die Kugeln stark an; die Bindungskraft ist 2,4-mal größer als der Strahlungsdruck.
  Senkrecht polarisiert stoßen sie sich schwach ab (nebeneinander schwingende Dipole).
- Die Bindungskräfte sind entgegengesetzt gleich, und die Summe der Einzelkräfte trifft die Kraft auf die umschließende Kugel auf
  0,1–0,2 % (n = 6: 0,4 %). Randspuren und Parallelfläche stimmen auf 0,3–0,8 % überein.

## Grenzen und nächster Schritt

- Nur achirale Außenmedien (der Spannungstensor im Pasteur-Medium ist ein anderer); Anregung durch ebene Wellen.
- Stufe 2: Kraft auf ein kleines, auch chirales Teilchen in Dipolnäherung aus E, H und ihren Gradienten, einschließlich der
  enantioselektiven Kraft im Spaltfeld.
