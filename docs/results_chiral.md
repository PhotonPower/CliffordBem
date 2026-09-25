# Ergebnisse: chirale Kugel und Zirkulardichroismus (C++-Kern v0.6)

## Umsetzung

Pasteur-Medium mit D = ε E + iχ H, B = μ H − iχ E (AP 1, Schritt 1.6):

- **Innerer Randoperator** (`ChiralCauchyOperator`): Mit den Helizitätsprojektoren P± = (1 ± iI)/2 zerfällt
  das Feld in F± = P± F mit (∇ − ik±) F± = 0, k± = ω(√ε√μ ∓ χ). Da P± zentral sind, ist
  E₁ h = P₊ E_{k₊} h + P₋ E_{k₋} h wieder eine Involution. Pro chiralem Körper werden zwei innere
  H-Matrizen (k₊, k₋) aufgebaut.
- **Transmissionsabbildung** (`transmission_map`): Tangentialanteile wie achiral; die Normalanteile koppeln
  über (e_n, h_n)₁ = D₁ C₁⁻¹ C₂ D₂⁻¹ (e_n, h_n)₂ mit C_j = [[ε_j, iχ_j], [−iχ_j, μ_j]] (AP 1, Lemma Jchiral,
  an exakten Fresnel-Lösungen auf 10⁻¹⁶ verifiziert). Die C++-Matrix stimmt mit dem Prototyp auf
  2·10⁻¹⁶ überein (`test_chiral`).
- **Anregung:** zirkular polarisierte ebene Wellen p = u + i s v (`circular_polarization`), Extinktion
  σ_ext = 4π/k Im(p̄ · E_∞(d))/|p|² für komplexe Polarisation.
- **Referenz** (`tools/mie_chiral.py`): chirale Mie-Lösung in der Helizitätsbasis (innen
  Beltrami-Felder M ± N mit k_A = ω(√εμ + χ), k_B = ω(√εμ − χ), je Multipolordnung ein 4×4-System aus der
  Stetigkeit von E_tan und H_tan). Kontrollen: χ = 0 reproduziert `mie.py` auf alle Stellen;
  verlustfrei ist Q_ext = Q_sca auf 10⁻¹⁶; σ_s(χ) = σ_{−s}(−χ).

## Chirales Dielektrikum, ε₁ = 2,25, χ = 0,2, ωa = 1 (verlustfrei)

Mie: Q₊ = 0,355566, Q₋ = 0,130452, CD = Q₊ − Q₋ = 0,225114.

| Dreiecke | Unbekannte | Q₊ | Q₋ | CD | Abw. CD | GMRES ± |
|---:|---:|---:|---:|---:|---:|---:|
| 1 280 | 10 240 | 0,350019 | 0,128696 | 0,221323 | 1,7 % | 10/10 |
| 2 880 | 23 040 | 0,353085 | 0,129667 | 0,223418 | 0,75 % | 10/10 |
| 5 120 | 40 960 | 0,354167 | 0,130010 | 0,224157 | 0,43 % | 10/10 |
| 8 000 | 64 000 | 0,354672 | 0,130168 | 0,224504 | 0,27 % | 10/10 |

Extrapolation aus den drei feinsten Netzen (Ordnung p ≈ 2): Q₊ 6·10⁻⁵, Q₋ 5·10⁻⁵, CD 9·10⁻⁵ relativ zu Mie.

## Chirales Metall, ε₁ = −11 + 1,2i, χ = 0,05, ωa = 0,5

Mie: Q₊ = 0,584560, Q₋ = 0,595354, CD = −0,010794 (Eindeutigkeitsbedingung aus AP 1, Im k±² ≥ 0, erfüllt).

| Dreiecke | Unbekannte | Q₊ | Q₋ | CD | Abw. CD | GMRES ± |
|---:|---:|---:|---:|---:|---:|---:|
| 1 280 | 10 240 | 0,595144 | 0,606178 | −0,011034 | 2,2 % | 26/26 |
| 2 880 | 23 040 | 0,589213 | 0,600114 | −0,010901 | 0,99 % | 26/26 |
| 5 120 | 40 960 | 0,587166 | 0,598004 | −0,010838 | 0,41 % | 25/25 |
| 8 000 | 64 000 | 0,586208 | 0,597042 | −0,010834 | 0,37 % | 24/24 |

Extrapolation: Q₊ 1,7·10⁻⁴, Q₋ 6·10⁻⁵. Der CD ist hier eine Differenz zweier fast gleicher Größen
(1,8 % von Q); sein absoluter Fehler liegt bei 4·10⁻⁵, die Folge ist für eine Extrapolation nicht glatt genug.

## Bewertung

- Der chirale Löser trifft die chirale Mie-Lösung mit derselben Ordnung 2 wie der achirale.
- Die Iterationszahlen sind von N unabhängig und gleich denen des achiralen Falls (10 bzw. 24–26).
- Kosten: drei H-Matrizen statt zwei (bei 64 000 Unbekannten 2,0 GB; Aufbau 110 s; beide Polarisationen
  zusammen 37–81 s Lösen).
- Nicht umgesetzt: chirale Außenmedien, Blockvorkonditionierung mit chiralem Innenmedium (braucht die
  helizitätsweisen dichten Blöcke).
