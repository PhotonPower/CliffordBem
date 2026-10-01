# Gekrümmte Elemente: Messung (Stufe 1, v0.44), lineare Dichten (Stufe 2a, v0.45), quadratische Geometrie (Stufe 1b, v0.46)

**Frage.** Lohnen gekrümmte Elemente? Vorab war geschätzt worden, dass der Fehler der Kugelstreuung zu rund 90 % aus der
Geometrie (eingeschriebenes Polyeder) stammt und gekrümmte Elemente ihn etwa zehnfach senken. Grundlage war ein Vergleich
mit Mie für eine Kugel gleichen Volumens. Diese Schätzung hat sich **nicht** bestätigt: Exakte Geometrie allein halbiert den
Fehler bei Glas und verdoppelt ihn bei Gold. Erst zusammen mit unstetig linearen Dichten sinkt der Fehler um zwei bis drei
Größenordnungen.

Skripte: `prototype/curved/curved_galerkin.py` (Messung), `prototype/curved/analyze.py` (Extrapolation, Tabellen);
Rohdaten: `results/curved_sphere_glass.csv`, `results/curved_sphere_gold.csv`. Der C++-Kern ist unverändert.

## Methode

Die Ikosaederkugel n entsteht durch radiale Projektion eines baryzentrischen Gitters auf den Ikosaederflächen. Das Gitter
n·m enthält alle Punkte des Gitters n, und die groben Kanten liegen auf denselben Großkreisen. Jedes gekrümmte grobe
Element J (Kugeldreieck) wird daher **exakt** von m² feinen ebenen Dreiecken i parkettiert.

- **Teilraum.** Stückweise konstante Dichten auf den gekrümmten Elementen: Basis P_iJ = √(a_i/A_J) (i ∈ J), A_J = Σ a_i,
  Spalten orthonormal. Unstetig lineare Dichten: je Element und Blade drei Funktionen (baryzentrische Koordinaten der
  Parametrisierung), auf dem feinen Netz als Mittelwerte je feinem Dreieck dargestellt (Mittelwert einer linearen Funktion
  = Wert im Schwerpunkt), je Element orthonormiert.
- **Galerkin-Projektion.** T_grob = Pᵀ T_fein P mit dem Transmissionsoperator T₁ des Kerns auf dem feinen Netz. Das erfasst
  die im Element variierende Normale, die im Cauchy-Operator und in der Transmissionsabbildung J steht. Rechte Seite
  Pᵀ b_fein, Streuspur P(h − b), Extinktion aus dem Fernfeld auf dem feinen Netz. GMRES mit dem Operator als
  Python-Funktion (`cb.gmres`), vorkonditioniert mit Pᵀ M_fein P.
- **Ebene Vergleichsgeometrie.** Die feinen Punkte werden zentral auf die Ebene ihres groben Dreiecks projiziert. Großkreise
  werden dabei zu Geraden, die Unterteilung bleibt konform; so entstehen ebene Elemente, gleich fein integriert.
- **Extrapolation m → ∞.** Die feinen Polyeder (gekrümmte Geometrie) und die Darstellung linearer Funktionen durch feine
  Konstanten nähern den Grenzfall mit einem Restfehler, der wie 1/m² fällt. Extrapoliert wird mit Q_∞ + c/m² + d/m⁴ aus den
  drei feinsten Punkten; die Abweichung zur Extrapolation mit c/m² aus zwei Punkten ist die angegebene Unsicherheit. Feine
  Netze bis 8 000 Dreiecke (Speichergrenze 3 GB), also m ≤ 5 bei 320, m ≤ 4 bei 500 und m ≤ 3 bei 720 groben Elementen.

Fälle wie in `results_scattering.md`: Glas ε = 2,25, ωa = 1 (Mie Q_ext = 0,215098); Gold ε = −11 + 1,2i, ωa = 0,5
(Mie Q_ext = 0,590018). Ebene Welle in z-Richtung, x-polarisiert; Voreinstellungen von `HMatrixParams` und `EntryParams`,
GMRES-Toleranz 10⁻⁸.

### Gegenproben

| Prüfung | Ergebnis |
|---|---|
| m = 1 gegen `ScatteringProblem` | bitgleich; bei 1 280 Elementen 0,211964 (Glas) und 0,600720 (Gold) wie in `results_scattering.md` |
| ebene Unterteilung, konstante Dichten, m = 1, 2, 3 (320 Elemente, Glas) | −5,6805 %, −5,6800 %, −5,6799 %: dieselbe Rechnung bis auf die Quadratur (5·10⁻⁶) |
| 1/m²-Gesetz (Glas, 320, gekrümmt/konstant) | aus m = 1, 2 vorhergesagt m = 3: +1,81 %, gemessen +1,85 %; Folge −5,68 / +0,64 / +1,85 / +2,28 / +2,48 % |
| jedes grobe Element enthält m² feine Dreiecke | geprüft (Zuordnung über den Kegel über dem Kugeldreieck) |

## Ergebnisse

Relativer Fehler von Q_ext gegen Mie; „eben/konstant“ ist die heutige Rechnung (m = 1), die übrigen Spalten sind nach
m → ∞ extrapoliert (Unsicherheit in Klammern; „2 P.“: nur zwei Punkte, ohne Unsicherheitsangabe).

**Glas** (ε = 2,25, ωa = 1)

| Elemente | eben/konstant (heute) | gekrümmt/konstant | eben/linear | gekrümmt/linear |
|---:|---:|---:|---:|---:|
| 320 | −5,681 % | +2,832 % (±0,002) | −6,163 % (±0,001) | **+0,013 %** (±0,002) |
| 500 | −3,680 % | +1,830 % (±0,001) | −4,000 % (±0,0001) | **+0,005 %** (±0,001) |
| 720 | −2,573 % | +1,278 % (±0,002) | −2,80 % (2 P.) | +0,001 % (2 P.) |
| 1 280 | −1,457 % | +0,72 % (2 P.) | – | – |

**Gold** (ε = −11 + 1,2i, ωa = 0,5)

| Elemente | eben/konstant (heute) | gekrümmt/konstant | eben/linear | gekrümmt/linear |
|---:|---:|---:|---:|---:|
| 320 | +7,306 % | +16,26 % (±0,02) | −5,79 % (±0,01) | **−0,073 %** (±0,006) |
| 500 | +4,684 % | +10,28 % (±0,01) | −3,87 % (±0,01) | **−0,023 %** (±0,002) |
| 720 | +3,246 % | +7,03 % (±0,01) | −2,73 % (2 P.) | −0,007 % (2 P.) |
| 1 280 | +1,814 % | +3,83 % (2 P.) | – | – |

**Modellabhängigkeit für gekrümmt/linear.** Der Grenzwert entsteht aus der Auslöschung großer Zahlen (Glas, 320 Elemente:
−1,45 / −0,65 / −0,36 / −0,23 % für m = 2…5). Mit acht Fehlermodellen (1/m²; 1/m² + 1/m⁴; 1/m + 1/m²; 1/m² + 1/m³;
jeweils mit und ohne m = 2) liegt er bei

| | 320 Elemente | 500 Elemente |
|---|---:|---:|
| Glas | +0,007 … +0,022 % | +0,003 … +0,008 % |
| Gold | −0,120 … −0,047 % | −0,037 … −0,017 % |

Die Aussage „höchstens etwa 0,02 % (Glas) bzw. 0,12 % (Gold) bei 320 Elementen“ hängt daher nicht am Modell.

## Deutung

1. **Die stückweise konstanten Dichten haben einen eigenen, großen Fehler.** Mit exakter Geometrie bleibt er allein übrig:
   +2,8 % (Glas) und +16 % (Gold) bei 320 Elementen. Bei Gold ist er besonders groß, weil die Oberflächenladung der
   plasmonischen Antwort stark variiert.
2. **Das eingeschriebene Polyeder kompensiert ihn teilweise.** Der Geometriefehler (sichtbar in eben/linear: −6,2 % bzw.
   −5,8 %) hat das entgegengesetzte Vorzeichen. Die beiden Fehler addieren sich aber nicht: −6,16 % + 2,83 % = −3,33 % gegen
   gemessen −5,68 % (Glas, 320); bei Gold +10,5 % gegen +7,3 %. Die Kompensation ist qualitativ, keine Zerlegung.
3. **Darum war die Volumenabschätzung irreführend.** Sie unterstellte, dass die BEM auf dem Polyeder fast exakt rechnet;
   tatsächlich enthält der Rest nach der Volumenkorrektur den Dichtefehler. Für Gold sagte sie die Verschlechterung durch
   exakte Geometrie sogar richtig voraus (Rest +3,6 % bei 1 280 Elementen, gemessen +3,8 %).
4. **Lineare Dichten allein helfen nicht**: Auf ebenen Elementen bleibt der Geometriefehler, der Fehler wird betragsmäßig
   bei Glas etwas größer (−6,2 % statt −5,7 %), bei Gold etwas kleiner (−5,8 % statt +7,3 %).
5. **Beides zusammen** senkt den Fehler bei 320 Elementen von 5,7 % auf höchstens 0,02 % (Glas, Faktor ≥ 250) und von
   7,3 % auf höchstens 0,12 % (Gold, Faktor ≥ 60).

**Konvergenzordnung.** Konstante Dichten konvergieren in allen drei Varianten mit O(h²) = O(N⁻¹) (gekrümmt/konstant, Glas:
2,83 % → 0,72 % von 320 auf 1 280 Elemente, Faktor 3,9). Gekrümmt/linear fällt deutlich schneller: von 320 auf 500
Elemente (N × 1,56) um den Faktor 2,6 (Glas) bzw. 3,2 (Gold), verträglich mit O(h⁴) = O(N⁻²). Das passt zur doppelten
Konvergenzordnung von Fernfeldfunktionalen in der Galerkin-BEM (h^{2(p+1)}: h² für konstante, h⁴ für lineare Dichten).
Die Ordnung ist aus zwei bis drei Netzen geschätzt und mit Vorsicht zu lesen.

**Aufwand, grob.** Gekrümmt/linear mit 320 Elementen hat 320 · 3 · 8 = 7 680 Unbekannte. Für denselben Fehler bräuchten
ebene Elemente mit konstanten Dichten (Fehler ≈ C/N, C aus 1 280 Elementen) bei Glas 90 000–150 000 Elemente (für
0,02 % bzw. 0,013 %), bei Gold je nach Modell 19 000–49 000 Elemente (0,12 % bzw. 0,047 %), also zwischen 150 000 und gut
einer Million Unbekannten. Eine Unbekannte der linearen Formulierung ist teurer
(Formfunktionen in der Nahfeldquadratur, mehr Kernkomponenten); der Vorteil bleibt dennoch groß, ist aber erst in Stufe 2
zu messen.

## Grenzen der Messung

- Nur die Kugel, je eine Frequenz, nur die Extinktion. Feldwerte im Nahfeld konvergieren punktweise langsamer als
  Fernfeldfunktionale; für Nahfeld und Kräfte ist der Gewinn eigens zu prüfen.
- Exakte Geometrie. Im Kern würde eine Näherung verwendet (z. B. quadratische Elemente aus Gmsh). Für interpolierende
  Flächen vom Grad k sind Flächen- und Volumenfehler O(h^{2k}), für k = 2 also O(h⁴); das reicht voraussichtlich für die
  volle Ordnung, ist aber nicht gemessen. Messbar wäre es mit demselben Verfahren: feine Punkte auf das quadratische
  Interpolationspolynom statt auf die Kugel legen (Stufe 1b).
- Extrapolation in m mit höchstens fünf Punkten; die Unsicherheit ist aus den Modellen geschätzt, nicht bewiesen.
- Würfel und andere Körper mit Kanten profitieren an ebenen Seiten nicht; dort dominiert die Kantensingularität.

## Konsequenz für Stufe 2

Gekrümmte Elemente lohnen sich **nur zusammen mit unstetig linearen Dichten**. Der Umbau betrifft:

- Geometrie: Elemente mit Parametrisierung (quadratisch), Jacobi-Determinante und variabler Normale.
- Eintragsauswertung: Die Normale steht im Integral, und es gibt drei Formfunktionen je Element. Dadurch entstehen mehr
  Kernkomponenten und größere H-Matrix-Blöcke.
- Nahfeldquadratur: Sauter-Schwab in Parameterkoordinaten mit Formfunktionen; die analytischen Dreiecksintegrale entfallen
  für gekrümmte Paare.
- Transmissionsabbildung als Elementmassenmatrix (24 × 24 je Element), Vorkonditionierer, rechte Seiten, Fernfeld,
  Nahfeld, Kräfte, Python-Anbindung.

Vor dem Umbau empfiehlt sich Stufe 1b (quadratische statt exakter Geometrie), weil sie entscheidet, ob quadratische Elemente
genügen.

## Stufe 2a: unstetig lineare Dichten im Kern (v0.45)

Der Kern kann jetzt mit unstetig linearen Dichten auf ebenen Dreiecken rechnen (`LinearScatteringProblem`, 24 Unbekannte je
Dreieck). Das ist die Hälfte der Kombination, die sich in Stufe 1 als lohnend erwiesen hat; die gekrümmte Geometrie folgt
in Stufe 2b. Auf ebenen Elementen allein bringt das an der Kugel erwartungsgemäß keinen Gewinn. Der Wert dieses Schritts
liegt darin, dass er sich unabhängig prüfen lässt.

**Aufbau.**
- Je Element wird eine orthonormierte Basis ψ_a = Σ_k S_ak λ_k verwendet, S = G^{−1/2}, analytisch:
  √(12/A)·1 + (√(3/A) − √(12/A))·𝟙𝟙ᵀ/3. Die Massenmatrix ist damit die Identität, und T₁ = ½(1+E₂) + ½(1−E₁)J
  behält seine Form; J wirkt auf ebenen Elementen blockweise.
- Die Einträge K^{ab}_c(i,j) haben 3×3 Formfunktionen × 4 Kernkomponenten je Elementpaar. Fern-, Nah-, Sauter-Schwab- und
  halbanalytische Zweige sind dieselben wie für konstante Dichten.
- Neu sind die analytischen Innenintegrale ∫λ_b/r und ∫λ_b(x−y)/r³ (`triangle_integrals_linear`). Sie werden über den Satz
  von Gauß in der Ebene auf Kantenintegrale von R, 1/R und ρ/R zurückgeführt.
- Im Selbstterm verschwindet der singuläre Anteil Φ₀ für lineare Dichten nicht. Er wird antisymmetrisiert integriert,
  ½[λ_a(x)λ_b(y) − λ_a(y)λ_b(x)], und ist so absolut integrierbar.
- H-Matrix: Der Clusterbaum läuft über die Elemente, mit je drei zusammenhängenden Indizes. Dichte Blöcke und ACA-Zeilen bzw.
  -Spalten werden je Elementpaar einmal ausgewertet. Der Pfad für konstante Dichten ist unverändert.

**Prüfungen** (`tests/test_linear.cpp`, Python `test_linear_densities`):

| Prüfung | Ergebnis |
|---|---|
| Innenintegrale gegen Brute-Force-Quadratur (28 672 Punkte; nahe Punkte: die Referenz konvergiert gegen die Formel) | 4,6·10⁻¹⁰ |
| Σ_ab K^{ab} = Einträge der konstanten Dichten (Summe der λ = 1), Kugel: fern, nah, Ecke, Kante, Selbstterm | ≤ 6·10⁻¹⁵ |
| dasselbe am gradierten Würfel; gestreckte Selbstterme halbanalytisch (Sauter-Schwab stagniert dort bei 3·10⁻⁵) | ≤ 4·10⁻¹⁵, Selbstterm 1,4·10⁻⁵ (fällt mit `sa_order`: 3,6·10⁻⁶ bei 20, 9,7·10⁻⁷ bei 28) |
| Mittelwert der linearen Projektion = konstante Projektion | 2,6·10⁻¹⁶ |
| Plemelj E b = b (ebene Welle, 320 Elemente, k = 0,665) | linear 4,9·10⁻⁴, konstant 4,6·10⁻³ |
| chirale Kugel: χ → 0 wie achiral; σ_s(χ) = σ_−s(−χ) | 10⁻⁹; 10⁻⁸ (ACA ε = 10⁻⁶; die Kompression erhält die Symmetrie nicht: 4·10⁻⁶ bei ε = 10⁻⁴) |

**Gegen die Vorhersage aus Stufe 1.** Stufe 1 hat eben/linear mit einem völlig anderen Verfahren berechnet: als
Galerkin-Projektion auf feinen Unterteilungen mit konstanten Dichten, extrapoliert in m. Beide Rechnungen stimmen bei
Fehlern von 4–6 % auf 0,002–0,02 Prozentpunkte überein:

| Fall | Stufe 1 (extrapoliert) | Stufe 2a (direkt) |
|---|---:|---:|
| Glas, 320 Elemente | −6,163 % (±0,001) | −6,165 % |
| Glas, 500 Elemente | −4,000 % | −4,001 % |
| Gold, 320 Elemente | −5,79 % (±0,014) | −5,81 % |
| Gold, 500 Elemente | −3,87 % (±0,012) | −3,88 % |

**Kosten** (Gold, voreingestellte Parameter): Die Zahl der Unbekannten verdreifacht sich. Die H-Matrizen werden bei 320
Dreiecken 7,2-mal und bei 1 280 Dreiecken 4,9-mal so groß (89,7 statt 12,5 MB bzw. 540 statt 110 MB). Das Verhältnis
nähert sich mit wachsendem Netz dem Faktor 3 der Niedrigrangteile, weil der Nahfeldanteil (Faktor 9) relativ abnimmt.
GMRES braucht etwa gleich viele Iterationen (26 gegen 28–29).

**Umfang von 2a.** Unterstützt sind ein oder mehrere Körper (auch chiral), ein achirales Außenmedium, ebene Wellen und
beliebige rechte Seiten, Extinktion und Vorwärtsamplitude, jeweils auch in Python. Noch nicht unterstützt sind das chirale
Außenmedium (Fernfeld je Helizität), Nahfeld und Kräfte, Block- und HODLR-Vorkonditionierung sowie die geschichteten
Probleme.

**Nächster Schritt.** Stufe 2b bringt gekrümmte Elemente im Kern. Nach Stufe 1b (unten) genügen dafür quadratische
Elemente. Dazu gehören Parametrisierung, Jacobi-Determinante und die Normale im Integral. Prüfziel ist
quadratisch/linear aus Stufe 1b: −0,020 % (Glas) und −0,12 % (Gold) bei 320 Elementen.

## Stufe 1b: quadratische statt exakter Geometrie (v0.46)

**Frage.** Erhält eine quadratisch interpolierte Geometrie die schnelle Konvergenz der exakten Kugel? Diese Geometrie bekäme
der Kern später aus Gmsh (`Mesh.ElementOrder 2`).

**Verfahren.** Das Verfahren von Stufe 1 bleibt, nur die Geometrieabbildung wird getauscht (`--geometry quadratic`):

- Parameter sind wie zuvor die Kegelkoordinaten λ über dem groben Element.
- Exakte Geometrie: X(λ) = Σλ_k V_k / |Σλ_k V_k|.
- Quadratische Geometrie (6-Knoten-Element): X(λ) = Σ λ_k(2λ_k − 1) V_k + Σ 4λ_aλ_b M_ab. Die Ecken V_k liegen auf der
  Kugel, die Kantenmitten M_ab = (V_a + V_b)/|V_a + V_b| ebenfalls, wie bei Gmsh.
- Dichtebasis, feine Unterteilung und Extrapolation sind unverändert. Der Vergleich trennt damit allein die Wirkung der
  Geometrienäherung ab.

**Gegenproben.**

| Prüfung | Ergebnis |
|---|---|
| Konformität: Punkte auf groben Kanten, von beiden Nachbarn berechnet | 10⁻¹⁵ |
| m = 1 (nur Ecken) gleich der heutigen Rechnung | 0,202879 wie zuvor |
| radialer Abstand der quadratischen Fläche von der Kugel (m = 4) | 1,07·10⁻⁴ (320 El.) → 6,9·10⁻⁶ (1 280 El.), Faktor 15,6 |
| Volumenfehler der quadratischen Fläche (in m extrapoliert, da das feine Polyeder selbst facettiert ist) | −2,2·10⁻³ (80) → −1,46·10⁻⁴ (320) → −9,5·10⁻⁶ (1 280): O(h⁴); eben: O(h²) |

Der Volumenfehler fällt wie O(h⁴), wie von der Theorie für Interpolation vom Grad k erwartet (O(h^{2k})).

**Ergebnisse** (Fehler von Q_ext gegen Mie, nach m → ∞ extrapoliert; Spanne über acht Fehlermodelle wie in Stufe 1):

| Fall | heute (eben, konstant) | exakt, linear | quadratisch, linear |
|---|---:|---:|---:|
| Glas, 320 Elemente | −5,681 % | +0,013 % (+0,007 … +0,022) | **−0,020 %** (−0,027 … −0,012) |
| Glas, 500 Elemente | −3,680 % | +0,005 % (+0,003 … +0,008) | **−0,003 %** (−0,006 … +0,005) |
| Glas, 720 Elemente | −2,573 % | +0,001 % (2 P.) | −0,003 % (2 P.) |
| Gold, 320 Elemente | +7,306 % | −0,073 % (−0,120 … −0,047) | **−0,118 %** (−0,182 … −0,080) |
| Gold, 500 Elemente | +4,684 % | −0,023 % (−0,037 … −0,017) | **−0,039 %** (−0,065 … −0,028) |
| Gold, 720 Elemente | +3,246 % | −0,007 % (2 P.) | −0,012 % (2 P.) |

**Deutung.**
- **Quadratische Elemente genügen.** Bei Gold fällt der Fehler von 320 über 500 auf 720 Elemente um die Faktoren 3,0 und
  3,2. O(h⁴) entspräche bei diesen Netzschritten den Faktoren 2,4 und 2,1; gemessen ist sogar etwas mehr, wie bei
  exakter Geometrie (3,2). Gegenüber heute ist der Fehler bei 320 Elementen 60-mal (Gold) bzw. fast 300-mal (Glas) kleiner.
- **Der Preis der Näherung** ist ein Fehler, der bei Gold etwa 1,6-mal so groß ist wie mit exakter Geometrie; bei Glas
  liegt er in derselben Größenordnung mit anderem Vorzeichen. Die Differenz zur exakten Geometrie passt zum Volumenfehler:
  Bei Glas mit 320 Elementen sagt σ ∝ V² für den Volumenfehler −1,46·10⁻⁴ etwa −0,03 % voraus, gemessen sind −0,033 %.
- **Höhere Geometrieordnung** (kubische Elemente) würde bei festem Netz höchstens diesen Faktor 1,6 bringen. Sie lohnt den
  Aufwand nicht; für Stufe 2b genügen quadratische Elemente.

**Grenzen.** Wie in Stufe 1: nur die Kugel, je eine Frequenz, nur die Extinktion; Extrapolation mit höchstens fünf Punkten.
Die Kantenmitten liegen hier exakt auf der Kugel. Bei Gmsh liegen sie auf der CAD-Fläche, bei Netzen ohne CAD-Beschreibung
muss die Fläche anders rekonstruiert werden.
