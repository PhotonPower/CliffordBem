# Ergebnisse: Kreuzapproximation mit Multivektor-Pivots (H3, erster Teil; C++-Kern v0.11)

## Varianten

Komprimiert wird der Kern K(i,j) = ∬ Φ_k(x − y), ein Paravektor s + **v** je Dreieckspaar; der Randoperator ist
E_ij x_j = −2/√(|τ_i||τ_j|) · K_ij (**n**_j x_j) (geometrisches Produkt).

| Modus | Darstellung eines Niedrigrangblocks | Speicher je Rangschritt |
|---|---|---:|
| `comp` | vier getrennte skalare ACAs | (m + n) je Komponente |
| `joint` | K_c(i,j) ≈ Σ_k u_k(i) w_{k,c}(j): skalare Zeilenfaktoren, Paravektor-Spaltenfaktoren | m + 4n |
| `mv` (neu) | K(i,j) ≈ Σ_k u_k(i) w_k(j) über Cl₃(ℂ), Kreuz K ← K − K_{:,j} K_ij⁻¹ K_{i,:} | 8(m + n) |

Die Multivektor-ACA (`AcaMode::Multivector`) wählt als Pivot unter den acht normgrößten Einträgen der Restzeile
den ersten invertierbaren (Nullteiler werden über det L(K_ij) erkannt, `mv_inverse`), aktualisiert mit dem
geometrischen Produkt (Reihenfolge beachtet) und bricht mit dem üblichen Kriterium
‖u_k‖‖w_k‖ ≤ ε ‖S_k‖ ab (Multivektornorm ⟨A Ã*⟩₀). Das Matrix-Vektor-Produkt nutzt die Algebra direkt:
y_i += Σ_k u_k(i) t_k mit t_k = Σ_j w_k(j) z_j, z_j = **n**_j x_j/√|τ_j|. Test: `test_hmatrix` (Fehler ≤ 10 ε).

## Vergleich bei gleicher Genauigkeit

Relativer Fehler des Matrix-Vektor-Produkts gegen exakte Zeilen; „NR“ = Speicher der Niedrigrangblöcke
(die dichten Nahblöcke sind für alle Varianten gleich).

**Kugel, 2 880 Dreiecke, k = 1,5**

| Variante | ε | Fehler | mittl. Rang | NR-Speicher | Mat-Vek |
|---|---:|---:|---:|---:|---:|
| mv | 10⁻² | 1,7·10⁻⁴ | 3,9 | 164 MB | 0,31 s |
| mv | 10⁻³ | 2,8·10⁻⁵ | 5,3 | 224 MB | 0,44 s |
| mv | 10⁻⁴ | 2,9·10⁻⁶ | 7,3 | 308 MB | 0,60 s |
| joint | 10⁻³ | 9,1·10⁻⁵ | 5,1 | 67 MB | 0,12 s |
| joint | 10⁻⁴ | 1,0·10⁻⁵ | 6,9 | 91 MB | 0,14 s |
| joint | 10⁻⁵ | 9,1·10⁻⁷ | 9,6 | 127 MB | 0,17 s |

**Gleichmäßiger Würfel, 3 072 Dreiecke** (Kanten und Ecken): mv 3,4·10⁻⁶ mit 297 MB gegen joint 1,4·10⁻⁶ mit 127 MB.
**Kugel bei k = 12** (kD = 24): mv 1,7·10⁻⁴ mit 280 MB gegen joint 4,1·10⁻⁵ mit 128 MB.

## Bewertung

- Die Multivektor-ACA braucht für dieselbe Genauigkeit **2,3- bis 2,9-mal mehr Speicher** als die gemeinsame
  skalare ACA, und das Matrix-Vektor-Produkt ist 3- bis 4-mal langsamer. Das gilt für glatte und kantige
  Geometrien und bei höherer Frequenz.
- Der Rang über der Algebra ist nur wenig kleiner als der skalare Rang der gemeinsamen ACA (z. B. 5,3 gegen 6,9
  bei vergleichbarem Fehler), während jeder Schritt acht statt einer bzw. vier Zahlen je Zeile und Spalte speichert.
- Erklärung: Die Niedrigrangstruktur des Kerns ist **skalar separabel**. Eine Fernfeldentwicklung um die
  Clustermittelpunkte hat die Form Σ_α (x − c)^α ∂^α Φ(c − y)/α!, also skalare Funktionen von x mal
  paravektorwertige Funktionen von y. Genau das ist die Darstellung der gemeinsamen ACA (skalare Zeilenfaktoren,
  Paravektor-Spaltenfaktoren). Multivektorwertige Zeilenfaktoren fügen Freiheitsgrade hinzu, die der Kern nicht
  braucht.
- **H3, erster Teil:** In der Form „Multivektor-Pivots senken den Speicher“ ist die Hypothese für die Kompression
  des Cauchy-Kerns widerlegt. Die Aussage, die sich halten lässt: Die gemeinsame Behandlung der Grade mit
  skalaren Zeilen- und Paravektor-Spaltenfaktoren ist die günstigste der untersuchten Darstellungen (etwa 19 %
  weniger als komponentenweise, 2,3–2,9-mal weniger als über der Algebra). Offen bleibt, ob Multivektor-Pivots
  für Operatoren mit echt multivektorwertigen Faktoren (etwa das Gesamtsystem T₁ mit J) oder als stabilere
  Pivotstrategie Vorteile bringen.
