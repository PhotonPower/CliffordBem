# Gekrümmte Elemente: Messung (Stufe 1, v0.44), lineare Dichten (Stufe 2a, v0.45), quadratische Geometrie (Stufe 1b, v0.46), gekrümmte Elemente im Kern (Stufe 2b, v0.47), Gmsh-Netze zweiter Ordnung (v0.48), schnellere Nahquadratur (v0.49), schnellere H-Matrix (v0.50), schnellere Geometrie (v0.51), Gauß-Regeln in der Nahquadratur (v0.52), anisotropes Sauter-Schwab (v0.53), schnelleres H-Matrix-Produkt (v0.54), ACA-Toleranz (v0.56), gepaartes Sauter-Schwab und H-Matrix (v0.57), Nahfeld und Kräfte (v0.58), chirales Außenmedium (v0.59), Nahfeldkosten und Python-Felder (v0.60)

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

## Stufe 2b: gekrümmte Elemente im Kern (v0.47)

Der Kern rechnet jetzt auf quadratischen (gekrümmten) Elementen mit unstetig linearen Dichten (`CurvedScatteringProblem`,
24 Unbekannte je Element). Das ist die Kombination, die sich in Stufe 1 als lohnend erwiesen hat.

**Ergebnis** (Fehler von σ_ext gegen Mie, ebene Welle, Voreinstellungen):

| Elemente | Glas, heute | **Glas, gekrümmt** | Gold, heute | **Gold, gekrümmt** |
|---:|---:|---:|---:|---:|
| 320 | −5,68 % | **−0,012 %** | +7,31 % | **−0,0067 %** |
| 500 | −3,68 % | **−0,0048 %** | +4,68 % | **−0,0023 %** |
| 720 | −2,57 % | **−0,0027 %** | +3,25 % | **−0,0010 %** |
| 1 280 | −1,46 % | – | +1,81 % | −0,00006 % |

Bei 320 Elementen ist der Fehler etwa 500-mal (Glas) bzw. 1 000-mal (Gold) kleiner als heute. Gold mit 320 gekrümmten
Elementen ist genauer als die heutige Rechnung mit 11 520 Elementen (+0,19 %). Der Fehler fällt etwa wie O(h⁴): bei Gold
von 320 auf 500 auf 720 Elemente um die Faktoren 2,9 und 2,3, erwartet sind 2,4 und 2,1. Bei 1 280 Elementen liegt Gold
nahe am Vorzeichenwechsel bzw. am Rauschboden der Kompression (6·10⁻⁷ relativ); eine Ordnung ist daraus nicht abzulesen.
Strengere Parameter (ACA ε = 10⁻⁶, Sauter-Schwab-Ordnung 8, Nahquadratur 0,2/0,1) ändern Gold bei 320 Elementen nur von
−0,0067 % auf −0,0068 %.

**Aufbau.**
- `QuadraticMesh`: 6-Knoten-Elemente wie Gmsh `ElementOrder 2` (Ecken des ebenen Netzes, Kantenmitten auf der Fläche),
  Parametrisierung über dem Referenzdreieck, dS = |X_u × X_v| du dv, Normale je Punkt. Die Basis ist je Element orthonormal
  bezüglich der gekrümmten Fläche: ψ = L⁻¹λ mit der Cholesky-Zerlegung der numerisch integrierten Gram-Matrix.
- `CurvedKernelEntries`: Die Normale variiert im Element und steht deshalb im Integral. G(z)·n(y) = s·n + v·(z·n) +
  v·(z∧n) hat sieben Komponenten (Skalar, Vektor, Bivektor), nicht 4×3 = 12. Fernpaare werden mit Gauß 7×7 im
  Parameterraum gerechnet, benachbarte Paare mit Sauter-Schwab in Parameterkoordinaten. Im Selbstterm wird
  Φ₀(x−y)n(y) = K_a + K_s zerlegt: K_a = Φ₀(n(x)+n(y))/2 ist antisymmetrisch und wird als Hauptwert mit antisymmetrisierten
  Gewichten integriert, K_s = Φ₀(n(y)−n(x))/2 ist schwach singulär. Nahe, getrennte Paare werden doppelt adaptiv
  integriert; analytische Innenintegrale gibt es auf gekrümmten Elementen nicht.
- `CurvedHMatrix` für sieben Komponenten (Baum über Elementen, je drei Indizes, Joint-ACA), `CurvedCauchyOperator`.
- `CurvedTransmissionOperator`: J(n(x)) variiert im Element und bildet lineare Dichten nicht auf lineare ab. Die exakte
  Komposition E₁J wäre ein Kern mit 64 Komponenten. Verwendet wird die Galerkin-Projektion J_G (24×24 je Element). Der
  Prototyp hat das vorab gemessen: Der Zusatzfehler beträgt etwa 20 % des Diskretisierungsfehlers (Glas, 320 Elemente:
  +0,005 %; Gold: +0,023 %) und fällt selbst wie O(h⁴). Vorkonditionierung 2(1 + J_G)⁻¹ je Element.

**Prüfungen** (`tests/test_curved.cpp`, Python `test_curved_elements`):

| Prüfung | Ergebnis |
|---|---|
| ebene Elemente (Mitten auf den Sehnen): Flächen, Normalen, Basis | 10⁻¹⁵ |
| Volumen der quadratischen Kugel, 320 Elemente | −1,458·10⁻⁴ (Prototyp Stufe 1b: −1,460·10⁻⁴), O(h⁴) (Faktor 15,7) |
| ebene Gegenprobe der Einträge gegen `LinearKernelEntries`: fern, Ecke, Kante, Selbstterm | ≤ 4,4·10⁻¹⁵ |
| dasselbe für nahe, getrennte Paare | 5,8·10⁻⁶; das ist der Fehler der Referenz (Rest mit 7 Gauß-Punkten); die doppelt adaptive Quadratur konvergiert in sich auf 2·10⁻⁷ |
| ebene Gegenprobe des Streuproblems gegen `LinearScatteringProblem` | 1,9·10⁻⁶ |
| Plemelj auf der gekrümmten Fläche: E h = h (innere Lösung), E h = −h (äußere, Dipol) | O(h²) (Faktor 2,2 bei 320 → 720) |
| Gegentest: Selbstterm ohne K_s | Fehler 40-mal größer und nur O(h): der Plemelj-Test erkennt einen falschen Selbstterm |
| chirale Kugel: σ_s(χ) = σ_−s(−χ) | 10⁻⁸ |

**Abweichung von der Vorhersage bei Gold.** Für Glas trifft der Kern die Vorhersage aus Stufe 1b einschließlich der
J_G-Korrektur (vorhergesagt etwa −0,015 % bei 320 Elementen, gemessen −0,012 %). Für Gold ist er etwa zehnmal besser als
vorhergesagt (etwa −0,09 % vorhergesagt, gemessen −0,007 %). Drei Befunde sprechen für den Kern und gegen die Vorhersage:

1. Der Kern ist numerisch konvergiert; strengere Parameter ändern das Ergebnis nicht.
2. Der Fehler geht wie O(h⁴) gegen null. Ein Fehlerterm im Kern würde typischerweise stagnieren oder langsamer fallen.
3. Der Plemelj-Test bestätigt den gekrümmten Operator, und ein absichtlich verfälschter Selbstterm fällt darin deutlich auf.

Die Vorhersage für Gold beruhte auf einer Extrapolation, die große Zahlen auslöschen musste (+1,82 → +0,79 → +0,42 % für
m = 2…4, extrapoliert auf −0,07 bis −0,12 %). Ihre Unsicherheit war offenbar größer als die angegebene Spanne über die
Fehlermodelle. Endgültig bewiesen ist das nicht.

**Kosten.** Gekrümmte Elemente mit 320 Elementen brauchen auf einem Kern 27 s Aufbau, davon 25 s für die doppelt adaptive
Nahquadratur, und 155 MB. Die heutige Rechnung braucht für einen vergleichbaren Fehler (Gold +0,19 % bei 11 520
Elementen) deutlich mehr Speicher und Zeit. Die Nahquadratur ist der erste Ansatzpunkt für eine Optimierung, etwa über
Singularitätssubtraktion mit den analytischen Integralen des Sehnendreiecks.

**Umfang von 2b.** Wie in 2a: ein oder mehrere Körper (auch chiral), achirales Außenmedium, ebene Wellen und beliebige
rechte Seiten, Extinktion, Vorwärtsamplitude und Fernfeld, jeweils auch in Python. Noch nicht unterstützt: das chirale
Außenmedium, Nahfeld und Kräfte, Block- und HODLR-Vorkonditionierung sowie geschichtete Körper. Gmsh-Netze zweiter
Ordnung liest v0.48 (unten). Netze mit Kanten und Ecken (Würfel) profitieren an den ebenen Seiten
nicht.

## Gmsh-Netze zweiter Ordnung (v0.48)

`read_gmsh_quadratic` liest 6-Knoten-Dreiecke (Gmsh-Elementtyp 9) aus MSH 2.2 und 4.1 (ASCII) als `QuadraticMesh`. Die
Kantenmitten liegen dann auf der CAD-Fläche, nicht auf einer nachträglichen Projektion. Körper werden wie bisher nach
physikalischer Gruppe bzw. Entität getrennt und über das Vorzeichen des Volumens nach außen orientiert. Beim Umdrehen
(v₀, v₂, v₁) werden die Kantenmitten zu (M₂₀, M₁₂, M₀₁) mitgetauscht. Andere Elementtypen (Punkte, Linien und
Volumenelemente, auch zweiter Ordnung) werden übersprungen, parametrische Knotenkoordinaten (`Mesh.SaveParametric`)
gelesen und verworfen. `read_gmsh` liest aus denselben Dateien die Ecken; bisher brach es bei Typ 9 im Format 4.1 ab.
`write_gmsh22_quadratic` schreibt konforme Netze zweiter Ordnung (eine Mitte je Kante). Netze erster Ordnung lehnt
`read_gmsh_quadratic` mit einem Hinweis ab (mit `Mesh.ElementOrder = 2` vernetzen oder `make_quadratic` verwenden).

**Prüfungen** (`test_gmsh`, Python `test_gmsh_second_order`) an einer echten Gmsh-Datei (`data/meshes/sphere_r1_order2_gmsh41.msh`:
CAD-Kugel aus OpenCASCADE, Gmsh 4.15.2, 320 Elemente, mit Punkten, Linien zweiter Ordnung und parametrischen Koordinaten):

| Prüfung | Ergebnis |
|---|---|
| Kantenmitten auf der Kugel | \|M\| − 1 = 2·10⁻¹⁶ (CAD-genau) |
| Formate 2.2 und 4.1, mit und ohne physikalische Gruppen | gleiche Netze |
| Volumen, quadratisch / eben | −1,54·10⁻⁴ / −3,50·10⁻²; von 254 auf 540 Elemente fällt der quadratische Fehler um 4,5 = (540/254)², also O(h⁴) |
| Rundreise mit zwei Körpern, einer absichtlich innen orientiert | Mitten und Volumen exakt |
| Netz erster Ordnung | abgelehnt |

**Streuung auf echten Gmsh-Netzen** (Fehler gegen Mie; „eben/konstant“ aus den Ecken derselben Datei):

| Fall | Elemente | gekrümmt | eben/konstant |
|---|---:|---:|---:|
| Glas, ωa = 1 | 254 | −0,017 % | −7,26 % |
| Glas | 320 | −0,010 % | −5,83 % |
| Glas | 540 | −0,0044 % | −3,53 % |
| Gold, ε = −11 + 1,2i, ωa = 0,5 | 254 | −0,028 % | +10,31 % |
| Gold | 540 | −0,0023 % | +4,57 % |

**Ein empfindlicher Fall.** Eine Goldkugel mit 20 nm Durchmesser in Wasser bei 700 nm (ωa = 0,09, fern der Resonanz)
hat eine winzige Extinktion (Q = 0,018), die fast reine Absorption ∝ Im α ist. Mit |Re ε|/Im ε = 15,5 verstärken sich
relative Fehler der Polarisierbarkeit etwa um diesen Faktor:

| Elemente (Ikosaederkugel) | gekrümmt | eben/konstant |
|---:|---:|---:|
| 320 | +0,72 % | +72 % |
| 720 | +0,145 % | +32 % |
| 2 880 | – | +8,0 % |

Gekrümmt fällt der Fehler um den Faktor 5,0 bei vorhergesagten (720/320)² = 5,06, also O(h⁴); eben fällt er mit O(h²).
Fern der Resonanz absorbierende Metallteilchen sind damit ein Fall, in dem die gekrümmten Elemente den Unterschied
zwischen brauchbar und unbrauchbar ausmachen.

**Beispiel** `examples/python/gmsh_curved.py`: ein Gold-Nanostäbchen (Zylinder mit Halbkugelkappen, 50 × 20 nm) in Wasser
bei 700 nm, longitudinal angeregt, mit Gmsh vernetzt. Bei 196 bzw. 366 Elementen liegen ebene Elemente 29 % bzw. 15 % über
dem feinsten gekrümmten Wert, gekrümmte ändern sich nur um 0,3 %. Ohne das Python-Paket gmsh rechnet das Skript die
mitgelieferte Kugel.

## Schnellere Nahquadratur (v0.49)

**Ausgangslage.** Bei 320 gekrümmten Elementen gingen 25 von 27 s Aufbau in die Nahquadratur, und davon 88 % in die nahen,
nicht benachbarten Elementpaare: 12 660 Paare zu je 600 µs. Sauter-Schwab für die benachbarten Paare (Ecke 150 µs, Kante
370 µs, Selbstterm 730 µs) war mit zusammen 12 % unkritisch. Die ebenen linearen Elemente brauchen für dieselben Paare
21 µs, weil ihr Innenintegral analytisch ist.

**Verfahren: Singularitätssubtraktion am Tangentialdreieck.** Für jeden äußeren Punkt x wird der Fußpunkt u* auf dem
gekrümmten Element bestimmt (Projektion auf das Sehnendreieck, zwei Newton-Schritte). Dort wird der singuläre Kern über
dem Tangentialdreieck X_aff(u) = X(u*) + X_u(u − u*) + X_v(v − v*) mit der Normalen n(u*) und der Jacobi-Determinante
J(u*) analytisch integriert (`triangle_integrals_linear`, die Gewichte λ_b sind exakt linear) und abgezogen. Im Rest
stimmen die Positionen bei u* bis O(|u − u*|²), Normale und Jacobi-Determinante bis O(|u − u*|); er verhält sich dort wie
1/r statt 1/r³ und wird mit einem gröberen adaptiven Kriterium integriert. Ein Sehnendreieck statt des Tangentialdreiecks
genügt nicht: Die Normale variiert über ein Element um etwa h/R (15 % bei 320 Elementen).

**Genauigkeit und Zeit der Einträge** (nahe, getrennte Paare; Fehler gegen das streng gerechnete doppelt adaptive
Verfahren, 0,12/0,06):

| Verfahren | 320 Elemente | µs je Paar | 1 280 Elemente | µs je Paar |
|---|---:|---:|---:|---:|
| doppelt adaptiv 0,3/0,15 (bis v0.48) | 3,0·10⁻⁶ | 593 | 4,5·10⁻⁶ | 624 |
| **Subtraktion 0,3/0,3 (Voreinstellung)** | 5,2·10⁻⁶ | 250 | 4,8·10⁻⁶ | 263 |
| Subtraktion 0,35/0,35 | 8,7·10⁻⁶ | 156 | 1,3·10⁻⁵ | 170 |
| Subtraktion 0,5/0,5 (schnell) | 2,3·10⁻⁵ | 50 | 3,0·10⁻⁵ | 57 |

**Wirkung auf die Streurechnung:**

| Fall | Aufbau bis v0.48 | Aufbau v0.49 | Änderung von σ_ext |
|---|---:|---:|---:|
| Gold, 320 Elemente | 19,2 s (Nahfeld 17,6 s) | 9,7 s (8,3 s) | 2·10⁻⁶ relativ |
| Glas, 320 Elemente | 18,5 s (17,0 s) | 10,0 s (8,5 s) | 4·10⁻⁸ |
| Gold, 1 280 Elemente | 88,8 s (64,6 s) | 55,1 s (32,1 s) | 9·10⁻⁷ |
| Gold, 320, schnell 0,5/0,5 | | 4,9 s (3,4 s) | 2·10⁻⁵ |

Die Voreinstellung halbiert den Aufbau bei unveränderter Genauigkeit. Die schnelle Einstellung ist so genau wie die
analytische Nahquadratur der ebenen Elemente (deren Voreinstellung liegt bei 1,7·10⁻⁵). Bei Gold mit 320 Elementen
verschiebt sie σ_ext aber um etwa 30 % des Diskretisierungsfehlers (−0,0045 % statt −0,0069 %). Für gekrümmte Elemente
ist sie deshalb nicht die Voreinstellung.

**Was nicht geholfen hat.**
- Die Gewichte λ_b n J um u* linear abzuziehen (als Kombination der λ_k, ebenfalls analytisch) änderte nichts
  (2,3·10⁻⁵ wie zuvor).
- Den Rest in Polarkoordinaten um u* zu integrieren (die Jacobi-Determinante hebt die 1/r-Spitze auf), war ab 4 Punkten je
  Richtung konvergiert, aber teurer als die adaptive Korrektur.
- Begrenzend ist die äußere Integration. Mit dem Kriterium 0,5 liegt sie bei 2·10⁻⁵, wie bei den ebenen Elementen.
  Beide Varianten wurden wieder entfernt.

**Prüfung** (`test_curved`, Teil 5): Subtraktion gegen das streng gerechnete doppelt adaptive Verfahren,
Voreinstellung 4,7·10⁻⁶, schnell 2,1·10⁻⁵. Python: σ_ext mit Subtraktion und doppelt adaptiv auf 10⁻⁶ gleich. Das
doppelt adaptive Verfahren bleibt als Referenz wählbar (`CurvedNearParams(subtract=False)`).

**Nächster Engpass.** Bei 1 280 Elementen braucht die H-Matrix mit sieben Komponenten jetzt 23 von 55 s, die
Nahquadratur 32 s.

## Schnellere H-Matrix mit sieben Komponenten (v0.50)

**Messung** (Gold ε = −11 + 1,2i, ωa = 0,5, 1 280 Elemente, ein Kern, Windows/MinGW g++ 16; je Wellenzahl k = 0,5 bzw.
k = 0,09 + 1,66i). Die Zahlen gelten für diesen Rechner: Hier braucht die H-Matrix nur 14 von 56 s Aufbau, die
Nahquadratur 42 s. Die Aufteilung im vorigen Abschnitt („23 von 55 s“) stammt von einem anderen Rechner.

| Anteil je H-Matrix (k = 0,5) | Zeit |
|---|---:|
| Fernblöcke in der ACA (1,06 Mio. Elementpaare, davon 0,83 Mio. verschieden) | 4,0 s |
| dichte Blöcke (281 000 Paare, davon 66 000 nah aus dem Cache) | 0,9 s |
| Nachkompression | 0,9 s |
| übrige ACA-Algebra | 1,0 s |
| gesamt | 6,8 s |

- Die ACA wertet 61 % aller 1,36 Mio. zulässigen Elementpaare aus: Die Blöcke sind klein (im Mittel 16 × 16 Elemente),
  jede ACA-Zeile kostet eine ganze Elementzeile, und bei Rang 7 sind nach wenigen Schritten fast alle Paare berechnet.
  Zeilen- und Spaltenzwischenspeicher waren getrennt, 22 % der Paare wurden doppelt berechnet.
- Ein Fernblock (49 Punktpaare) kostete 3,0 µs: Kern 1,45 µs (komplexes `exp`), Komponenten 0,36 µs, Akkumulation in
  3 × 3 × 7 Einträge und Umrechnung in die ψ-Basis (`to_psi`) etwa 1,2 µs.

**Änderungen.**
- Zeilen- und Spaltenzwischenspeicher der ACA lesen Elementpaare, die schon in der anderen Richtung berechnet sind, von
  dort (`CurvedHMatrix`).
- `block_far` rechnet direkt in der ψ-Basis: Die Gewichte w ψ_a an den Quadraturpunkten werden je Element einmal
  berechnet, die Kontraktion läuft zweistufig (erst über y, dann über x; 21 statt 63 Produkte je Punktpaar), `to_psi`
  entfällt.
- `dirac_kernel_fast`: e^{ikr} = e^{−Im k r}(cos Re k r + i sin Re k r) in reeller Arithmetik, 22 statt 29 ns je
  Auswertung, Abweichung ≤ 7·10⁻¹⁶. Nur im gekrümmten Pfad (Fern-, Nah-, Sauter-Schwab-Quadratur); der konstante Pfad
  behält `dirac_kernel_full` und bleibt bitgleich.

**Ergebnis** (ein Kern):

| Fall | Aufbau v0.49 | Aufbau v0.50 | H-Matrizen v0.49 → v0.50 | Änderung von σ_ext |
|---|---:|---:|---:|---:|
| Glas, 320 Elemente | 11,6 s | 10,6 s | | 0 (12 Stellen) |
| Gold, 320 Elemente | 11,7 s | 10,8 s | | 0 (12 Stellen) |
| Gold, 1 280 Elemente | 56,1 s | 49,0 s | 14,4 → 10,3 s | 1,3·10⁻⁸ relativ |

Mit 12 Threads braucht der ganze Aufbau bei 1 280 Elementen 9,1 s (H-Matrizen 1,9 s, Nahquadratur 7,3 s).

**Prüfung** (`test_curved`, Teil 6): `block_far` gegen S L Sᵀ der λ-Blöcke 9,5·10⁻¹⁶; `dirac_kernel_fast` gegen
`dirac_kernel_full` 6,8·10⁻¹⁶. Alle übrigen Werte von `test_curved` und `test_gmsh` unverändert.

**Plattform.** Auf diesem Rechner liefert schon v0.49 bei Gold mit 320 Elementen −0,0065 % gegen Mie statt −0,0069 %
(Linux). Der Unterschied von 4·10⁻⁶ liegt in der Größenordnung der Nahquadratur.

**Nächste Schritte.** Weiter ließe sich die H-Matrix durch eine ACA auf den Quadraturpunkten beschleunigen (etwa 4-mal
weniger Kernauswertungen, dann dominiert die ACA-Algebra; geschätzt Faktor 1,5). Der Speicher (475 MB je H-Matrix bei
1 280 Elementen, davon 280 MB in dichten Blöcken) ist ein eigener Punkt. Auf diesem Rechner ist der größere Hebel die
Nahquadratur: 39 von 49 s; der schnellere Kern ändert dort fast nichts, die Zeit geht in die analytischen Integrale und
die Geometrie.

## Schnellere Geometrie in der Nahquadratur (v0.51)

**Messung** (Gold, 1 280 Elemente, ein Kern, je Wellenzahl 19,3 s Nahquadratur): Sauter-Schwab für benachbarte Paare
5,1 s (31 Mio. Punkte in 16 580 Paaren), nahe getrennte Paare 14,3 s. Dort laufen an 2,8 Mio. äußeren Punkten
Fußpunkt und analytische Integrale (419 ns je Aufruf), dazu 107 Mio. Korrekturpunkte (≈ 38 je äußerem Punkt) und
28 Mio. Unterteilungsknoten. Ein Korrekturpunkt kostete etwa 100 ns, davon 42 ns Geometrie: `QuadraticMesh::jacobian`
(29 ns) und `QuadraticMesh::X` (13 ns) sammeln die Knoten je Aufruf über die Konnektivität und werten die Formfunktionen
getrennt aus.

**Änderung.** `CurvedKernelEntries` speichert je Element die ausmultiplizierte Geometrie
X(u, v) = A + B u + C v + D u² + E uv + F v² (u = λ₁, v = λ₂) und wertet X, Tangenten und Jacobi-Determinante inline aus –
in Fern-, Nah-, Sauter-Schwab-Quadratur, Unterteilung und Fußpunktsuche. `QuadraticMesh` bleibt unverändert.

**Ergebnis** (ein Kern):

| Fall | Aufbau v0.50 | Aufbau v0.51 | Nahquadratur v0.50 → v0.51 | Änderung von σ_ext |
|---|---:|---:|---:|---:|
| Glas, 320 Elemente | 10,6 s | 8,2 s | 9,7 → 7,4 s | 9·10⁻¹⁰ |
| Gold, 320 Elemente | 10,8 s | 8,4 s | 10,0 → 7,6 s | 5·10⁻⁹ |
| Gold, 1 280 Elemente | 49,0 s | 39,5 s | 38,7 → 29,4 s | 6·10⁻¹⁰ |

σ_ext ändert sich um mehr als die Rundung, weil die Rundung einzelne adaptive Unterteilungsentscheidungen kippt; das liegt
drei Größenordnungen unter dem Fehler der Nahquadratur (5·10⁻⁶). Mit 12 Threads braucht der Aufbau bei 1 280 Elementen
7,0 s (Nahquadratur 5,1 s, H-Matrizen 1,8 s), zu Beginn von v0.50 waren es 9,1 s. `test_curved` und `test_gmsh`
bestehen unverändert (ebene Gegenprobe der Einträge weiterhin ≤ 2·10⁻¹⁴).

**Was die Korrektur begrenzt** (nahe getrennte Paare, 320 Elemente, k = 1,3 + 0,05i, Fehler gegen das streng gerechnete
doppelt adaptive Verfahren):

| äußere Regel / Korrektur | größter Fehler | mittlerer Fehler | µs je Paar |
|---|---:|---:|---:|
| 0,3 / 0,3 (Voreinstellung) | 5,2·10⁻⁶ | 1,7·10⁻⁶ | 278 |
| 0,3 / 0,4 | 1,1·10⁻⁵ | 5,7·10⁻⁶ | 171 |
| 0,3 / 0,5 | 1,9·10⁻⁵ | 1,0·10⁻⁵ | 138 |
| 0,4 / 0,4 | 1,3·10⁻⁵ | 6,6·10⁻⁶ | 108 |
| 0,5 / 0,5 | 2,3·10⁻⁵ | 1,2·10⁻⁵ | 51 |

Bei der Voreinstellung begrenzt die Korrektur, nicht die äußere Regel: Schon 0,3/0,4 verdoppelt den Fehler. Der Rest
nach der Subtraktion verhält sich wie 1/r, getragen von zwei Termen gleicher Ordnung (Variation von n J und Krümmung
X − X_aff). Weniger Korrekturpunkte bei gleicher Genauigkeit verlangen, beide abzuziehen; den Abzug von n J allein hatte
v0.49 ohne Gewinn versucht.

## Gauß-Regeln in der Nahquadratur (v0.52)

**Frage.** Nach v0.51 begrenzte die Korrektur die Genauigkeit der Nahquadratur. Vermutet war der Rest nach der
Subtraktion, der sich wie 1/r verhält. Er hat zwei Anteile gleicher Ordnung: die Variation von N = X_u × X_v = J n und
die Krümmung X − X_aff = q(δ) = Dδu² + Eδuδv + Fδv².

**Prototyp Stufe 2** (`prototype/curved/second_order_subtraction.cpp`, nicht im Kern): zusätzlich abgezogen
λ_b(u*) [K₀(z_a)(N_u δu + N_v δv) − DK₀(z_a)[q(δ)] N*] mit DK₀(z)[q] = (q − 3(ẑ·q)ẑ)/(4πr³). Der erste Teil folgt aus
`triangle_integrals_linear`, ∫ δ_k K₀ du = (Ig_k − u*_k ΣIg_b)/(4π J*); der Krümmungsvektor wurde für die Messung fein
numerisch integriert. Gemessen wurde das Innenintegral an 18 140 Punkten x (320 Elemente, k = 1,3 + 0,05i) gegen eine
Referenz, die gegen die ungeteilte adaptive Quadratur auf 1,2·10⁻¹⁰ stimmt:

| Blattregel (ohne Unterteilung) | ohne Subtraktion | Stufe 1 (v0.49) | Stufe 2 |
|---|---:|---:|---:|
| Dunavant 7 (Grad 5) | 9,9·10⁻⁴ | 1,9·10⁻⁴ | 1,0·10⁻⁴ |
| Gauß 4 × 4 (16 Punkte) | 9,9·10⁻⁵ | 2,4·10⁻⁵ | 1,4·10⁻⁵ |
| Gauß 5 × 5 (25 Punkte) | 6,8·10⁻⁶ | 1,7·10⁻⁶ | 1,5·10⁻⁶ |
| Gauß 6 × 6 (36 Punkte) | 5,1·10⁻⁷ | 1,4·10⁻⁷ | 1,3·10⁻⁷ |

Dunavant 7 mit der Unterteilung 0,3 braucht 30 Punkte für 8,4·10⁻⁶. Der Krümmungsterm bringt nur einen Faktor 1,4: Nicht
die Singularität am Fußpunkt begrenzt, sondern der Grad der Blattregel. Der Integrand ist auf der Skala des Abstands
glatt, dort sind Regeln hohen Grades viel wirksamer als Unterteilung. Stufe 2 wurde deshalb nicht eingebaut.

**Änderung.** `QuadRule::conical(n)`: Gauß-Legendre n × n auf dem Quadrat mit Duffy-Abbildung (exakt bis Grad 2n − 2).
`CurvedNearParams::outer_rule` und `correction_rule` wählen die Blattregel der äußeren Integration und der Korrektur
(0 = Dunavant 7). Neue Voreinstellung: Gauß 5 × 5, Kriterien 1,5/1,5. Die Adaptivität bleibt, damit sehr nahe Punkte
(gradierte, unregelmäßige Netze) weiter unterteilt werden.

**Abstimmung** (`prototype/curved/near_rules.cpp`; nahe getrennte Paare gegen Gauß 8 × 8 mit 0,3/0,3; diese Referenz stimmt mit dem streng gerechneten
doppelt adaptiven Verfahren auf 0,7–1,2·10⁻⁷):

| Einstellung | Ikosaeder 320 | Ikosaeder 1 280 | Gmsh-Kugel 320 | µs je Paar |
|---|---:|---:|---:|---:|
| Dunavant 7, 0,3/0,3 (v0.49) | 5,2·10⁻⁶ | 4,8·10⁻⁶ | 9,5·10⁻⁶ | 217–244 |
| Dunavant 7, 0,5/0,5 (schnell bis v0.51) | 2,3·10⁻⁵ | 3,9·10⁻⁵ | 3,6·10⁻⁵ | 49–53 |
| **Gauß 5, 1,5/1,5 (Voreinstellung)** | 6,2·10⁻⁷ | 1,1·10⁻⁶ | 2,6·10⁻⁶ | 69–82 |
| Gauß 4, 1,0/1,0 (schnell) | 2,5·10⁻⁵ | 2,5·10⁻⁵ | 3,7·10⁻⁵ | 50–51 |
| Gauß 6, 1,5/1,5 | 2,6·10⁻⁸ | | 1,9·10⁻⁷ | 137–174 |
| Gauß 7, 2,0/2,0 | 5,1·10⁻⁹ | | | 233 |

Größere Kriterien als 1,5 ändern auf den Kugeln fast nichts mehr (dort wird dann kaum noch unterteilt).

**Wirkung auf die Streurechnung** (ein Kern; „streng“: Gauß 7 × 7 mit 1,0/1,0):

| Fall | σ_ext v0.51 gegen streng | v0.52 gegen streng | Aufbau v0.51 → v0.52 |
|---|---:|---:|---:|
| Glas, 320 Elemente | 4,6·10⁻⁸ | 2,6·10⁻⁹ | 7,8 → 3,9 s |
| Gold, 320 Elemente | 2,7·10⁻⁶ | 4,7·10⁻⁸ | 7,8 → 3,8 s |
| Gold, 1 280 Elemente | 9,5·10⁻⁷ | 1,5·10⁻⁸ | 36,9 → 22,6 s (Nahquadratur 27,6 → 12,7 s) |

Die Voreinstellung ist 20- bis 60-mal genauer und halbiert den Aufbau. Die schnelle Einstellung (Gauß 4) liegt bei
1,3·10⁻⁶ (Gold 320) und spart gegenüber der Voreinstellung nur noch 10–15 % des Aufbaus.

**Nebenbefund.** Gold mit 320 Elementen gibt jetzt −0,0067 % gegen Mie, mit der strengen Quadratur ebenso. Die früheren
Werte −0,0069 % (Linux) und −0,0065 % (Windows) unterschieden sich um den Quadraturfehler von v0.49 (2,7·10⁻⁶), der auf
beiden Plattformen verschieden ausfiel.

**Prüfung.** `test_curved` (5): Voreinstellung 6,1·10⁻⁷ (Schwelle 3·10⁻⁶, vorher 1·10⁻⁵), schnell 1,9·10⁻⁵, Dunavant 7
mit 0,3/0,3 4,7·10⁻⁶ gegen das streng gerechnete doppelt adaptive Verfahren. Python: Felder `outer_rule`,
`correction_rule`, Voreinstellung und schnelle Einstellung gegen doppelt adaptiv.

**Nächster Engpass.** Bei 1 280 Elementen braucht die Nahquadratur 12,7 s; davon entfallen jetzt etwa 5 s auf
Sauter-Schwab (benachbarte Paare). Die H-Matrizen brauchen 10 s.

## Anisotrope Sauter-Schwab-Regeln (v0.53)

**Messung** (`prototype/curved/ss_orders.cpp`; Einträge benachbarter Paare gegen isotrope Ordnung 12, Gold innen
k = 0,09 + 1,66i). Nach v0.52 war Sauter-Schwab die größte Fehlerquelle der Quadratur: Bei isotroper Ordnung 5 (bis v0.52)
lagen die Einträge auf dem Ikosaeder (320 Elemente) bei 2,9·10⁻⁶ (Ecke), 6,0·10⁻⁶ (Kante) und 2,6·10⁻⁵ (Selbstterm), σ_ext
von Gold um 9,8·10⁻⁷ neben der Ordnung 9 – 20-mal mehr als die Nahquadratur. Die isotrope Ordnung konvergiert etwa eine
Dekade je Ordnung; auf der Gmsh-Kugel nur etwa Faktor 4, mit Ausreißern bei verzerrten Elementen (längste Kante/Höhe
1,7–3,1, gleichseitig 1,15): dort Ecke bis 3,9·10⁻⁴, Selbstterm bis 3,5·10⁻⁴ bei Medianen wie auf dem Ikosaeder.

**Welche Richtung begrenzt.** Eine Richtung der Sauter-Schwab-Abbildung (ξ, η₁, η₂, η₃) von 5 auf 8 angehoben:

| Ordnungen (ξ, η₁, η₂, η₃) | Ecke | Kante | Selbstterm |
|---|---:|---:|---:|
| 5, 5, 5, 5 | 2,9·10⁻⁶ | 6,0·10⁻⁶ | 2,6·10⁻⁵ |
| 8, 5, 5, 5 | 2,9·10⁻⁶ | 6,0·10⁻⁶ | 2,6·10⁻⁵ |
| 5, 8, 5, 5 | 1,8·10⁻⁶ | 6,0·10⁻⁶ | 2,6·10⁻⁵ |
| 5, 5, 8, 5 | 2,8·10⁻⁶ | 6,1·10⁻⁶ | 2,6·10⁻⁵ |
| 5, 5, 5, 8 | 2,9·10⁻⁶ | 6,2·10⁻⁶ | **3,9·10⁻⁸** |

Der Selbstterm hängt allein an η₃; die Kante an η₂ und η₃ gemeinsam (4, 5, 7, 7: 3,5·10⁻⁷), die Ecke an η₁ und η₂ (5, 6, 6, 5:
2,2·10⁻⁷). ξ genügt mit 4 (Kante, Selbstterm) bzw. 5 (Ecke); ξ = 3 ist überall zu wenig. Die Drehung der Ecken des
Selbstterms (größter Winkel an die Referenzecke) ändert nichts, die Regel ist über ihre sechs Teilsimplizes symmetrisch.

**Änderung.** `PairRule::sauter_schwab(Adjacency, std::array<int, 4>)` mit eigener Ordnung je Richtung;
`CurvedNearParams::ss_orders` je Nachbarschaft, Voreinstellung Ecke 5, 6, 6, 5 (900 Punkte je Teilgebiet), Kante 4, 4, 6, 6
(576), Selbstterm 4, 4, 3, 7 (336); erste Zahl 0: isotrop mit `EntryParams::ss_order`. Die Pfade konstanter und linearer
Dichten bleiben isotrop.

| Typ (Ikosaeder 320) | isotrop 5 | anisotrop | µs je Paar |
|---|---:|---:|---|
| Ecke | 2,9·10⁻⁶ | 2,2·10⁻⁷ | 102 → 152 |
| Kante | 6,0·10⁻⁶ | 6,0·10⁻⁷ | 251 → 230 |
| Selbstterm | 2,6·10⁻⁵ | 3,8·10⁻⁷ | 683 → 358 |

Auf der Gmsh-Kugel sinkt der größte Fehler des Selbstterms von 3,5·10⁻⁴ auf 1,1·10⁻⁵, der Kante von 1,0·10⁻⁴ auf 2,0·10⁻⁵,
der Ecke von 3,9·10⁻⁴ auf 9,3·10⁻⁵ (die Ausreißer bleiben, aber seltener und kleiner).

**Wirkung auf die Streurechnung** (σ_ext gegen isotrope Ordnung 9):

| Fall | isotrop 5 (bis v0.52) | anisotrop (v0.53) | isotrop 6 |
|---|---:|---:|---:|
| Ikosaeder 320, Gold | 9,8·10⁻⁷ | 2,0·10⁻⁸ | 6,2·10⁻⁸ |
| Gmsh-Kugel 320, Glas | 9,6·10⁻⁷ | 3,5·10⁻⁸ | 4,8·10⁻⁷ |
| Gmsh-Kugel 320, Gold | 2,2·10⁻⁵ | 1,4·10⁻⁷ | 5,2·10⁻⁶ |

Die anisotrope Regel ist genauer als isotrop 6 und so teuer wie isotrop 5: Aufbau Gold 1 280 auf einem Kern 22,0 s
(isotrop 5) gegen 22,1 s. Auf der Gmsh-Kugel war der Sauter-Schwab-Fehler von Gold (2,2·10⁻⁵) bisher größer als alle
übrigen Quadraturfehler zusammen.

**Prüfung.** `test_curved` (7): anisotrop gegen isotrop 12 auf 180 Elementen: Ecke 1,7·10⁻⁷, Kante 1,2·10⁻⁶,
Selbstterm 1,0·10⁻⁶ (isotrop 5: 2,5·10⁻⁶, 5,2·10⁻⁶, 2,2·10⁻⁵). Die ebene Gegenprobe (2) gegen `LinearKernelEntries`
läuft mit der isotropen Regel, damit sie bitgenau bleibt.

## Schnelleres H-Matrix-Produkt (v0.54)

**Messung** (Gold, 1 280 Elemente, ein Kern; `prototype/curved/hmatrix_params.cpp`). Das Lösen kostete so viel wie der
ganze Aufbau: bei tol = 10⁻⁶ 23 Iterationen in 20,5 s, bei 10⁻¹⁰ 45 Iterationen in 38,7 s (Aufbau 22 s). Je Iteration
werden zwei H-Matrizen angewendet, je 0,42 s. Eine H-Matrix hält 481 MB (dichte Blöcke 283 MB); das Produkt lief damit
mit 1,1 GB/s, die Speicherbandbreite eines Kerns liegt bei 7,2 GB/s.

Ursachen:
- GCC prüft nach jeder komplexen Multiplikation `std::complex<double>` auf NaN und springt dann in `__muldc3` (C99
  Anhang G). Die innere Schleife des Produkts enthielt 48 solche Aufrufstellen; die Verzweigung verhindert die
  Vektorisierung.
- Die 8 Komponenten des Ergebnisses wurden nach jeder Operation zurückgeschrieben (`y` und `z` könnten für den Compiler
  denselben Speicher bezeichnen).

**Änderung** (`CurvedHMatrix::apply`): innere Schleife y += k z in reeller Arithmetik (ac − bd, ad + bc, die Rechnung
des schnellen Zweigs von `std::complex`), Ergebnis in lokalen Akkumulatoren. Die Reihenfolge der Additionen bleibt; die
Ergebnisse sind bitgleich (σ_ext von Gold auf 17 Stellen, gleiche Iterationszahlen).

| | v0.53 | v0.54 | v0.54 mit `-DCBEM_NATIVE=ON` |
|---|---:|---:|---:|
| Produkt, dichte Blöcke | 0,27 s | 0,127 s | 0,058 s |
| Produkt, niedrigrangige Blöcke | 0,15 s | 0,098 s | 0,047 s |
| Lösen Gold 320, tol 10⁻¹⁰ (47 It.) | 7,4 s | 3,8 s | |
| Lösen Gold 1 280, tol 10⁻¹⁰ (45 It.) | 38,7 s | 22,3 s | |

Mit `-march=native` (AVX2, FMA) halbiert sich das Produkt noch einmal; die dichten Blöcke laufen dann mit 4,9 GB/s nahe
der Bandbreite, und der Aufbau wird 15 % schneller (Ergebnisse wegen FMA nicht mehr bitgleich). `CBEM_NATIVE` bleibt
in der Voreinstellung aus (portable Programme und Wheels); für eigene Rechnungen lohnt es sich.

**Partition** (gleiche Messung, Fehler des Produkts gegen ε = 10⁻⁸):

| leaf / eta / sep_factor | Aufbau | Produkt | Speicher (dicht) | Fehler |
|---|---:|---:|---:|---:|
| 32 / 1 / 3 (Voreinstellung) | 5,2 s | 0,42 s | 481 MB (283) | 1,4·10⁻⁶ |
| 64 / 1 / 3 | 4,8 s | 0,57 s | 631 MB (491) | 9,1·10⁻⁷ |
| 32 / 2 / 3 | 4,1 s | 0,39 s | 438 MB (283) | 1,6·10⁻⁶ |
| 64 / 2 / 2 | 4,2 s | 0,36 s | 405 MB (252) | 2,3·10⁻⁶ |

(Zeiten vor der Änderung des Produkts.) `eta = 2` spart 20 % des Aufbaus bei kaum größerem Fehler. `HMatrixParams` gilt
für alle Pfade; die Voreinstellung bleibt deshalb.

**Nächste Schritte.** Dieselbe Umstellung für `KernelHMatrix` (konstante und lineare Dichten; 45 Aufrufstellen von
`__muldc3`), HODLR (23) und dichte Blöcke (34) – bitgleich, also ohne Änderung der Ergebnisse. Danach Speicherung in
einfacher Genauigkeit (complex64): halbiert den Speicher und, weil die dichten Blöcke mit `-march=native` an der
Bandbreite liegen, etwa auch das Produkt.

## ACA-Toleranz 10⁻⁵ für gekrümmte Elemente (v0.56)

Seit v0.52/v0.53 liegen die Quadraturfehler bei etwa 2·10⁻⁸ auf σ_ext; der Kompressionsfehler der ACA mit ε = 10⁻⁴
(1,5·10⁻⁶) dominierte. Messung (ein Kern, H-Matrix in einfacher Genauigkeit, Fehler gegen ε = 10⁻⁷):

| Gold 1 280 | Fehler σ_ext | Aufbau | Speicher | Lösen (tol 10⁻¹⁰) |
|---|---:|---:|---:|---:|
| ε = 10⁻⁴ | 1,5·10⁻⁶ | 21,6 s | 477 MB | 45 It. × 0,499 s = 22,5 s |
| ε = 10⁻⁵ | 1,6·10⁻⁷ | 22,8 s | 549 MB | 43 It. × 0,551 s = 23,7 s |
| ε = 10⁻⁶ | 1,7·10⁻¹⁰ | 25,5 s | 626 MB | 43 It. × 0,651 s = 28,0 s |

Bei 320 Elementen: 1,7·10⁻⁶ (10⁻⁴), 4,4·10⁻⁸ (10⁻⁵), 2,7·10⁻⁹ (10⁻⁶) bei gleicher Zeit. Voreinstellung für gekrümmte
Elemente jetzt 10⁻⁵ (`curved_hmatrix_params()`, Standardwert von `CurvedHMatrix`, `CurvedScatteringProblem` und der
Python-Anbindung); 10–40-mal genauer für etwa 5 % Zeit und 15 % Speicher. Die übrigen Pfade bleiben bei 10⁻⁴. Gold mit
320 Elementen gibt jetzt −0,0068 % gegen Mie (40 statt 47 Iterationen).

## Aufbau: gepaartes Sauter-Schwab, H-Matrix eta = 2 mit ε = 10⁻⁶ (v0.57)

**Messung** (`prototype/curved/build_profile.cpp`; Gold, 1 280 Elemente, ein Kern, Stand v0.56): Aufbau 24,1 s, davon
H-Matrizen 10,9 s (45 %), getrennte Nahpaare 6,5 s (27 %), Sauter-Schwab 6,3 s (26 %; davon Eckpaare 3,5 s, Kantenpaare
1,9 s, Selbstterme 0,9 s), Rest (Transmission, ψ, Quadratur) 0,05 s.

**Gepaartes Sauter-Schwab.** Eine Regel auf τ_i × τ_j integriert auch das Paar (j, i): x′ = y auf τ_j, y′ = x auf τ_i,
z′ = −z, die Normale gehört zu x. Für Ecken- und Kantenpaare werden beide Blöcke in einem Durchgang berechnet
(`lambda_sauter_schwab_pair`); Kern und Geometrie je Punkt nur einmal. Der Block (i, j) bleibt bitgleich, der Block (j, i)
ist so genau wie direkt gerechnet (gegen isotrope Ordnung 12: Ecke 1,7·10⁻⁷, Kante 1,2·10⁻⁶). Nahfeld-Cache 13,1 → 11,4 s.

**H-Matrix.** Partition bei ε = 10⁻⁵ (`prototype/curved/hmatrix_params.cpp`, Produkt mit Zufallsvektor gegen ε = 10⁻⁸):
eta = 2 spart 14 % Aufbau, 10 % Produkt und 11 % Speicher bei gleichem Fehler des Produkts (1,7·10⁻⁷ statt 1,5·10⁻⁷);
größere Blätter oder sep_factor = 2 bringen nichts oder verschlechtern den Fehler. σ_ext reagiert aber empfindlicher
(Gold, gegen eta = 1 mit ε = 10⁻⁷):

| eta / ε | Fehler 320 | Fehler 1 280 | Aufbau 1 280 | Speicher | Lösen 1 280 |
|---|---:|---:|---:|---:|---:|
| 1 / 10⁻⁵ (v0.56) | 4,4·10⁻⁸ | 1,6·10⁻⁷ | 21,8 s | 549 MB | 23,8 s |
| 2 / 10⁻⁵ | 5,9·10⁻⁷ | 5,9·10⁻⁷ | 20,2 s | 486 MB | 21,8 s |
| 1,5 / 10⁻⁵ | 5,9·10⁻⁷ | 5,5·10⁻⁷ | 20,4 s | 499 MB | 21,8 s |
| **2 / 10⁻⁶ (v0.57)** | 3,7·10⁻⁹ | 3,1·10⁻⁸ | 22,6 s | 548 MB | 23,8 s |
| 1,5 / 10⁻⁶ | 3,7·10⁻⁹ | 2,4·10⁻⁸ | 22,6 s | 565 MB | 24,4 s |

(Zeiten mit gepaartem Sauter-Schwab.) eta = 2 mit ε = 10⁻⁶ kostet so viel wie die Voreinstellung von v0.56 und ist 5- bis
12-mal genauer; die ACA ist damit so genau wie die Quadratur (≈ 2·10⁻⁸). Neue Voreinstellung in `curved_hmatrix_params()`.

**Ergebnis.** Aufbau Gold 1 280: 24,1 s (v0.56) → 23,0 s (Nahfeld 11,4 s), σ_ext 1,6·10⁻⁷ → 3·10⁻⁸ neben der Rechnung
mit ε = 10⁻⁷. Gold 320: −0,0068 % gegen Mie (40 Iterationen).

**Prüfung.** `test_curved` (7): gepaarte Blöcke (j, i) gegen isotrop 12, Block (i, j) bitgleich zur Einzelrechnung.

## Messung: H-Matrix-Aufbau gekrümmter Elemente (nach v0.57)

**Aufteilung** (Gold, 1 280 Elemente, ein Kern, eta = 2, ε = 10⁻⁶; je Wellenzahl 5,4–5,9 s): dichte Blöcke 0,6–0,7 s
(281 000 Elementpaare), Blockauswertungen der ACA 2,2–2,3 s (850 000 von 1,36 Mio. zulässigen Paaren, 63 %), ACA-Algebra
1,0–1,1 s, Nachkompression 1,5–1,7 s. Die Blöcke haben im Mittel 22 × 22 Elemente bei Rang 15; die ACA wertet deshalb
fast jedes Paar aus.

**Nachkompression.** Sie senkt den Rang von 18,6 auf 15,0, spart 11 % Speicher (311 → 278 MB) und 12 % Zeit je Produkt
(0,308 → 0,270 s) und kostet 1,7 s Aufbau je Wellenzahl. Bei 45 Iterationen mit zwei Produkten gleicht sich das aus; sie
bleibt.

**ACA auf Quadraturpunkten** (`prototype/curved/point_aca.cpp`; 30 zufällige zulässige Clusterpaare, ε = 10⁻⁶): Statt
Elementpaaren (7 × 7 Punkte, 7 Komponenten) werden einzelne Punktzeilen und -spalten der vier Kernkomponenten (s, v z)
gezogen.

| Clustergröße | Elementebene: Rang, Kernauswertungen, Zeit | Punktebene: Rang, Kernauswertungen, Zeit |
|---|---|---|
| 22 Elemente | 13,0; 23 700; 2,18 ms | 13,2; 4 480; 1,43 ms |
| 44 Elemente | 16,0; 63 200; 5,83 ms | 16,0; 11 100; 4,07 ms |

Gleicher Rang, fünfmal weniger Kernauswertungen, aber nur ein Drittel weniger Zeit: ACA-Algebra und Nachkompression
überwiegen dann. Mit der Umrechnung der Punktfaktoren auf die Basis (≈ 0,2 ms je Block) bliebe etwa −24 % auf den
ACA-Teil, −20 % auf die H-Matrizen und etwa −10 % auf den Aufbau (23 → 21 s). Nicht eingebaut; der Aufbau ist auf
H-Matrizen (≈ 11 s), getrennte Nahpaare (6,5 s) und Sauter-Schwab (≈ 4,6 s) verteilt, und mit 12 Threads dauert er bei
1 280 Elementen etwa 4 s.

## Nahfeld und Kräfte auf gekrümmten Elementen (v0.58)

**Verfahren.** Das Streufeld im Außenraum ist das Cauchy-Integral der Streuspur h_s = h − b mit der Normalen im Integral:
F_s(x) = Σ_τ Σ_a [∫_τ ψ_a(y) Φ_k(x − y) n(y) dS_y] u_{τ,a}, mit den Koeffizienten u_{τ,a} der 24N-Spur und demselben
Vorzeichen wie bei konstanten Dichten. Die Integrale in der Klammer sind die sieben Kernkomponenten der Einträge;
`CurvedKernelEntries::point_integrals` wertet sie fern (Abstand zum Schwerpunkt > 4 Umkreisradien) mit der 7-Punkt-Regel aus,
nah mit der Singularitätssubtraktion der nahen Paare (Fußpunkt, Tangentialdreieck, Gauß-Korrektur). Die Markierungen
„innen“ und „zu nah“ beziehen sich auf das Sehnennetz; als zu nah gilt zusätzlich der Bereich der größten Wölbung, in dem
ein Punkt zwischen Sehne und gekrümmter Fläche liegen kann. Die Kraftfunktionen (`force_on_sphere`, `force_on_offset`,
`fields_with_gradients`) arbeiten unverändert über den Feldauswerter `make_plane_wave_eval_curved` bzw.
`make_near_field_eval_curved`; `project_incident_curved` projiziert beliebige einfallende Felder (Strahlen, Dipole) auf die
gekrümmte Spur.

**Nahfeld gegen Mie** (Goldkugel ε = −11 + 1,2i in Wasser, ωa = 0,5, zirkular polarisiert; Abweichung von |E|² und der
optischen Chiralität C an vier Punkten in 0,05 bis 1 Radien Abstand):

| Elemente | (1,05, 0, 0) | (0, 0, 1,2) | (0,72, 0,58, 0,77) | (2, 0, 0) |
|---|---|---|---|---|
| konstant, 1 280 | −1,7 % / −0,3 % | −5,7 % / +2,4 % | −3,1 % / −0,8 % | −1,3 % / −0,02 % |
| gekrümmt, 180 | −0,60 % / −1,13 % | −0,80 % / +0,81 % | −0,48 % / −0,09 % | −0,15 % / +0,04 % |
| gekrümmt, 320 | +0,67 % / +0,16 % | −0,17 % / −0,08 % | −0,18 % / −0,05 % | −0,05 % / +0,01 % |
| gekrümmt, 720 | +0,16 % / +0,04 % | −0,04 % / +0,02 % | −0,04 % / −0,01 % | −0,01 % / +0,002 % |
| gekrümmt, 1 280 | +0,046 % / +0,012 % | −0,015 % / +0,008 % | −0,011 % / −0,003 % | −0,004 % / +0,001 % |

Bei 1 280 Elementen 40- bis 400-mal genauer als konstante Dichten; schon 320 gekrümmte Elemente sind genauer als 1 280
ebene. Von 320 auf 1 280 Elemente fällt der Fehler um den Faktor 14 (O(h⁴): 16).

**Strahlungsdruck gegen Mie** (linear polarisiert, Mie 10,629788036; Kugel R = 1,5 bzw. Parallelfläche des Sehnennetzes
im Abstand 0,2):

| Elemente | Kugel | Parallelfläche | Querkraft |
|---|---:|---:|---:|
| konstant, 1 280 | −1,52 % | −1,54 % | |
| gekrümmt, 180 | −0,357 % | −0,368 % | 2·10⁻⁹ |
| gekrümmt, 320 | −0,120 % | −0,123 % | 4·10⁻⁹ |
| gekrümmt, 720 | −0,025 % | −0,026 % | 1·10⁻⁸ |
| gekrümmt, 1 280 | −0,0081 % | −0,0083 % | 7·10⁻⁹ |

Bei 1 280 Elementen 190-mal genauer als konstante Dichten. Die Quadratur auf der Kugel (32 bzw. 48 Ringe) ändert den
Wert nicht.

**Kosten.** Direkte Summation über Punkte × Elemente: vier Punkte bei 1 280 Elementen 0,01 s, die Kraft über die Kugel
(32 × 64 Punkte) 0,3 s mit 12 Threads.

**Prüfung** (`test_curved_near_field`, Python `test_curved_elements`): allgemeine Projektion gegen die der ebenen Welle
(6·10⁻¹⁶), Fernfeldgrenze (7,6·10⁻⁴ bei R = 400 wie im ebenen Pfad), Nahfeld und Strahlungsdruck gegen Mie bei 320 und
720 Elementen mit Konvergenzordnung, Kugel gegen Parallelfläche, Querkraft, Markierungen (innen, Wölbung).

**Offen.** Kraft aus den Randspuren und aus dem Fernfeld (`force_from_traces`, `force_from_far_field`) für gekrümmte
Elemente; H-Matrix für viele Auswertepunkte (wie `NearFieldOperator`); chirales Außenmedium; einfallende Felder aus Python.

## Chirales Außenmedium für lineare und gekrümmte Elemente (v0.59)

**Verfahren** wie im Pfad konstanter Dichten (`results_chiral.md`, v0.21): Außenoperator P₊E_{k₊} + P₋E_{k₋} aus zwei
Cauchy-Operatoren (`ChiralCauchyOperator`), einfallende Helizitätswelle mit k_σ (`plane_wave_incidence`; lineare
Polarisation ist keine Eigenmode und wird abgelehnt), Extinktion und Vorwärtsamplitude nach dem optischen Theorem im
Kanal σ aus P_σ h_s. Die Transmissionsabbildung J (`transmission_map`, an jedem Quadraturpunkt) behandelte chirale Medien
auf beiden Seiten schon. Das Nahfeld auf gekrümmten Elementen rechnet das Streufeld je Helizität mit k_± (wie
`exterior_near_field`), die Kraftfunktionen folgen über den Feldauswerter.

**Goldkugel in chiralem Wasser** (ε = 1,7689, χ = 0,05; Gold ε = −11 + 1,2i; ωa = 0,5) gegen die chirale Mie-Lösung
(σ₊ = 11,765497, σ₋ = 11,912810, `tools/mie_chiral_layered.py`):

| Elemente | konstant σ₊ / σ₋ / CD | linear σ₊ / σ₋ / CD | gekrümmt σ₊ / σ₋ / CD |
|---|---|---|---|
| 180 | | −9,69 % / −9,84 % / −21,5 % | −0,367 % / −0,370 % / −0,63 % |
| 320 | −5,08 % / −5,28 % / −21,0 % | −5,71 % / −5,80 % / −12,6 % | −0,123 % / −0,124 % / −0,21 % |
| 720 | −2,37 % / −2,46 % / −9,7 % | −2,64 % / −2,68 % / −5,8 % | −0,026 % / −0,026 % / −0,044 % |
| 1 280 | −1,36 % / −1,41 % / −5,6 % | −1,51 % / −1,53 % / −3,3 % | −0,0083 % / −0,0084 % / −0,014 % |
| 2 880 | −0,61 % / −0,64 % / −2,5 % | | |

Der Zirkulardichroismus CD = σ₊ − σ₋ ist eine Differenz von 1,2 % der Extinktion; sein relativer Fehler ist im ebenen
Pfad viermal so groß wie der von σ. Gekrümmte Elemente geben ihn bei 320 Elementen 100-mal, bei 1 280 Elementen 400-mal
genauer als konstante Dichten (etwa O(h⁴)). Lineare Dichten auf ebenen Elementen verbessern σ nicht (die Geometrie
dominiert, wie in Stufe 2a), wohl aber den CD um den Faktor 1,7.

**Gegenproben** (`test_curved` (9), `test_linear` (7), `test_curved_near_field` (6)): Kugel aus dem Außenmedium
unsichtbar (σ = 8·10⁻¹⁷ bzw. 10⁻¹⁷), lineare Polarisation abgelehnt, Grenzfall χ = 10⁻⁹ gegen achiral (gekrümmt
1,8·10⁻¹⁰, linear 3,2·10⁻⁸, Nahfeld 4,5·10⁻⁹), Spiegelsymmetrie σ₊(χ) = σ₋(−χ) (gekrümmt 7·10⁻⁹ mit ε = 10⁻⁶,
linear 1,6·10⁻⁶ mit ε = 10⁻⁴).

## Nahfeld: Kosten großer Karten und einfallende Felder aus Python (v0.60)

**Kosten** (`prototype/curved/near_field_cost.cpp`; Goldkugel in Wasser, Karte in der Ebene y = 0, 12 Threads): Die
direkte Summation kostet etwa 0,1 µs je Punkt und Element (ein Kern etwa 0,7–0,9 µs).

| Elemente | Punkte | gekrümmt, direkt | eben (ab 2 000 Punkten H-Matrix) |
|---|---:|---:|---:|
| 1 280 | 9 800 | 1,6 s | 0,43 s |
| 5 120 | 10 000 | 6,6 s | 0,88 s |
| 5 120 | 40 000 | 27,4 s | 2,7 s |

Die Kraft über die Kugel (32 × 64 Punkte) braucht 0,4 s (1 280 Elemente) bzw. 1,5 s (5 120 Elemente). Da 1 280
gekrümmte Elemente genauer sind als 5 120 ebene, kostet eine Karte mit 40 000 Punkten bei gleicher Genauigkeit etwa 7 s
statt 2,7 s.

**Aufteilung** (ein Kern, 1 280 Elemente, 2 868 Punkte): nahe Paare (Singularitätssubtraktion) 13 800 mit zusammen
0,04 s; ferne Paare 3,7 Mio. zu je 0,47 µs für die Integrale (davon die 7 Kernauswertungen etwa die Hälfte), dazu das
Produkt mit der Dichte.

**Was nicht geholfen hat.**
- Die Dichte je Quadraturpunkt vorab einzurechnen (D_{τ,p,c} = Σ_a w ψ_a(y_p) e_c u_{τ,a}) und je Punkt und Element nur
  noch Kern, Komponenten und Produkte zu rechnen: 7 × 7 × 8 = 392 komplexe Produkte je Paar statt 7 × 3 × 7 + 3 × 7 × 8 =
  315, also mehr Arbeit (es gibt 7 Quadraturpunkte, aber nur 3 Basisfunktionen); ein Kern 0,9 statt 0,7 µs je Paar.
- Kachelung (Blöcke von 64 Punkten je Element) gegen den Speicherverkehr: keine Änderung, der Speicher begrenzt nicht.

Beides wurde zurückgenommen. Wesentlich schneller wird die Karte nur mit einer H-Matrix für die Auswertepunkte (wie
`NearFieldOperator` im ebenen Pfad).

**Einfallende Felder aus Python.** `CustomField.project` wählt für ein `QuadraticMesh` die Projektion auf die gekrümmte
Spur (`projection_points_curved`, `project_samples_curved`: dieselbe Regel wie `project_incident_curved`, ein
vektorisierter Aufruf von `fields(x)`); `exterior_near_field_curved` und `near_field_evaluator_curved` nehmen Python-Felder
wie der ebene Pfad (Streufeld ohne GIL, `fields(x)` einmal je Auswertung). Prüfung (Python `test_curved_elements`): eine
ebene Welle als Python-Feld gibt dieselbe Spur (10⁻¹³), dasselbe Nahfeld (10⁻¹²) und dieselbe Kraft (10⁻¹⁰) wie die des
Kerns.

## Messung: Vorkonditionierung gekrümmter Elemente (nach v0.60)

**Iterationen** (`prototype/curved/precond_baseline.cpp`; Kugel, ebene Welle; Vorkonditionierer 2(1 + J_G)⁻¹ je Element):

| Fall | Elemente | tol 10⁻⁶ mit / ohne | tol 10⁻¹⁰ mit / ohne |
|---|---:|---:|---:|
| Glas (ε = 2,25, ωa = 1) | 320 / 1 280 | 10 / 11 | 15 / 16 |
| Gold (ε = −11 + 1,2i, ωa = 0,5) | 320 | 26 / 28 | 46 / 48 |
| Gold | 1 280 | 23 / 27 | 43 / 47 |
| Gold in Wasser | 320 / 1 280 | 27 / 28, 25 / 27 | 49 / 50, 45 / 48 |

Die Iterationszahl hängt nicht vom Netz ab. Der elementweise Vorkonditionierer ändert sie kaum. Je Iteration kostet fast
nur der Operator (1 280 Elemente, 12 Threads: 0,12 s; Vorkonditionierer 0,001 s, GMRES-Verwaltung 0,01 s).

**Spektrum** (`prototype/curved/precond_spectrum.cpp`; Arnoldi mit 100 Schritten auf T M, 320 Elemente): Bei Gold liegen
die Ritz-Werte breit zwischen etwa 0,15 − 0,52i … 0,20 − 0,65i und 1,83 + 0,51i … 1,92 + 0,55i, ohne einzelne
Ausreißer. Das sind die Eigenwerte des Hauptsymbols von 2(1 + J)⁻¹T₁, die AP 1 für ε = −11 + 1,2i angibt
(0,19 − 0,54i; 0,47 + 0,81i; 1,54 − 0,81i; 1,81 + 0,54i). Die Iterationszahl ist also durch den Materialkontrast
festgelegt, nicht durch Netz oder Geometrie. Lokale Vorkonditionierer (Blöcke, HODLR) können daran nichts ändern; das
entspricht dem Befund für glatte plasmonische Körper im ebenen Pfad (`results_preconditioning.md`). Ändern ließe es sich
nur über den Operator selbst (Vorkonditionierung vom Calderón-Typ, jede Iteration etwa doppelt so teuer) oder, bei vielen
rechten Seiten, über Recycling und Deflation.

**Toleranz** (Fehler von σ_ext gegen tol 10⁻¹², 1 280 Elemente):

| tol | Glas | Gold | Gold in Wasser |
|---|---|---|---|
| 10⁻⁴ | 7 It., 2,4·10⁻⁴ | 16 It., 7,5·10⁻⁶ | 18 It., 2,4·10⁻⁶ |
| 10⁻⁶ (Voreinstellung) | 10 It., 2,6·10⁻⁷ | 23 It., 4,0·10⁻⁷ | 25 It., 1,0·10⁻⁷ |
| 10⁻⁷ | 11 It., 7,3·10⁻⁸ | 27 It., 9,3·10⁻⁹ | 30 It., 5,0·10⁻⁹ |
| 10⁻⁸ | 12 It., 1,7·10⁻⁹ | 32 It., 7,3·10⁻¹⁰ | 35 It., 2,1·10⁻¹⁰ |

Die Voreinstellung 10⁻⁶ hält σ_ext auf 4·10⁻⁷, weit unter dem Diskretisierungsfehler (gegen Mie etwa 7·10⁻⁵ bei
1 280 Elementen). Für Vergleiche auf Quadraturgenauigkeit (≈ 3·10⁻⁸) genügt 10⁻⁷; 10⁻¹⁰ (Tests, Messungen) kostet fast
die doppelte Zahl an Iterationen.

## Messung: Krylov-Recycling gekrümmter Elemente (nach v0.60)

`prototype/curved/recycling_rhs.cpp`: zwölf rechte Seiten wie bei der Orientierungsmittelung (sechs Einfallsrichtungen, je
zwei Polarisationen), nacheinander gelöst mit GMRES, `RecyclingGmres` (bis 120 Richtungen) und `GcroDr` (k = 20, m = 80);
Kugel, ωa = 0,5. Iterationen gesamt und Zeit (12 Threads):

| Fall | Elemente, tol | GMRES | Recycling 120 | GCRO-DR 20/80 |
|---|---|---:|---:|---:|
| Gold (ε = −11 + 1,2i) | 320, 10⁻⁶ | 312 It., 6,8 s | 288 It., 7,5 s | 299 It., 7,2 s |
| Gold | 320, 10⁻¹⁰ | 552 It., 11,9 s | 575 It., 14,6 s | 554 It., 12,9 s |
| Gold | 1 280, 10⁻⁶ | 276 It., 38,5 s | 273 It., 40,6 s | 289 It., 39,9 s |
| nahe Dipolresonanz (ε = −2,3 + 0,2i) | 320, 10⁻⁶ | 396 It., 8,6 s | 382 It., 10,2 s | 348 It., 8,4 s |
| nahe Dipolresonanz | 320, 10⁻¹⁰ | 763 It., 16,1 s | 811 It., 19,6 s | 641 It., 15,0 s |
| nahe Dipolresonanz | 1 280, 10⁻⁶ | 348 It., 46,2 s | 384 It., 57,3 s | 356 It., 49,0 s |

Recycling spart höchstens 16 % der Iterationen (GCRO-DR, Resonanz, 10⁻¹⁰) und kaum Zeit, weil die Orthogonalisierung gegen
die gespeicherten Richtungen hinzukommt; oft braucht es sogar mehr Iterationen als GMRES. Das folgt aus dem Spektrum (Messung
oben): Die Eigenwerte liegen breit auf den vier Ästen des Hauptsymbols, ohne einzelne Ausreißer, die ein kleiner Unterraum
abfangen könnte; die Dipolresonanz der Kugel hebt sich bei diesem Materialkontrast nicht genug ab. Im ebenen Pfad waren es
10 % bei ε = −11 und 30–40 % nahe einer Resonanz (`results_dipole.md`), bei groberen Netzen mit stärkeren Ausreißern.
`CurvedScatteringProblem` bekommt daher kein `use_recycling`; die Lösungen stimmen auf tol überein (Abweichung ≤ 2·10⁻⁶).

## Nahfeldkarten mit H-Matrix (v0.61)

**Ansatz.** Wie `NearFieldOperator` im ebenen Pfad: Zeilen sind die Auswertepunkte, Spalten die Elemente, Clusterbäume
über Punkte und Schwerpunkte. Ein Block ist zulässig, wenn min(diam) ≤ η dist, dist > 3 h_max, |k| diam ≤ 20 und
zusätzlich dist > 4 r_max (r = größter Abstand der Ecken vom Schwerpunkt plus Wölbung): Dann rechnet `point_integrals` für
alle Paare des Blocks mit der 7-Punkt-Fernregel, und das Feld des Blocks ist eine Summe über Quadraturpunkte,
F(x) = Σ_q (s + v z) W_q mit z = x − y_q und W_q = n_q Σ_a w_q ψ_a(y_q) u_a. Nahe Blöcke rechnen direkt mit `point_integrals`
(Singularitätssubtraktion wie bisher).

**Messung** (`prototype/curved/near_field_hmatrix.cpp`; Goldkugel in Wasser, Karte in der Ebene y = 0, 12 Threads; Fehler
gegen die direkte Summation bezogen auf max |F_s|). Zwei Varianten der ACA+ in zulässigen Blöcken:
- **E (Elementebene):** Einträge sind die 21 Integrale `point_integrals` je Punkt und Element (3 Basisfunktionen × 7
  Komponenten).
- **P (Punktebene):** Einträge sind die 4 Komponenten des Dirac-Kerns (s, v z) je Punkt und Quadraturpunkt (28 Spalten je
  Element).

| Fall (η = 2, eps 10⁻⁴) | Aufbau + Anwendung | Fehler | Rang | Kernauswertungen |
|---|---:|---:|---:|---:|
| 1 280 × 9 968, direkt | 1,61 s | – | – | 100 % |
| E, Blatt 32 | 0,51 + 0,04 s | 3,0·10⁻⁵ | 8,2 | 15,5 % |
| P, Blatt 32 / 64 / 128 | 0,41 / 0,43 / 0,45 s + 0,05–0,10 s | 2,5 / 2,3 / 6,7·10⁻⁵ | 8,2–8,9 | 5,9 / 5,2 / 4,5 % |
| P, Blatt 64, ohne Nachkompression | 0,23 + 0,08 s | 2,0·10⁻⁵ | 12,5 | 5,2 % |
| 5 120 × 39 912, direkt | 26,6 s | – | – | 100 % |
| E, Blatt 32 | 3,03 + 0,14 s | 3,5·10⁻⁵ | 8,8 | 5,5 % |
| P, Blatt 64 | 2,76 + 0,27 s | 4,0·10⁻⁵ | 9,5 | 1,7 % |
| P, Blatt 64, ohne Nachkompression | 1,83 + 0,40 s | 4,2·10⁻⁵ | 13,7 | 1,7 % |

Mit eps = 10⁻⁶ (P, Blatt 64): 1 280 × 9 968 0,89 + 0,08 s mit, 0,40 + 0,11 s ohne Nachkompression; 5 120 × 39 912
7,05 + 0,38 s mit, 3,13 + 0,52 s ohne; Fehler jeweils 2,6–2,8·10⁻⁷.

Was die Zeit bestimmt:
- **Arithmetik der ACA, nicht die Kernauswertungen.** P braucht nur 1,7 % der Kernauswertungen der direkten Summation.
  Fast die ganze Zeit geht in die Arithmetik der ACA+ (Reste, Kreuzterme, Referenzen: O(r²(m + 28 n)) je Block; 5 120 × 39 912:
  19 s Rechenzeit über alle Threads) und in die Nachkompression (QR und SVD, 12 s; bei eps = 10⁻⁶ 44 s).
- **Nachkompression.** Für eine einmalige Anwendung ist sie überflüssig: Ohne sie wird der Aufbau 35–55 % schneller, das
  Produkt etwas langsamer (höherer Rang), der Fehler bleibt.
- **Betragsquadrate statt `std::abs`.** In der Pivotsuche ersetzt (`std::abs` einer komplexen Zahl ruft `hypot`), sparen sie
  nur 5–20 % der ACA-Zeit; `aca.cpp` (gemeinsamer Code) bleibt unverändert.
- **E gegen P.** E ist ähnlich schnell, braucht aber dreimal so viele Kernauswertungen (jede Zeile rechnet `point_integrals`).
- **Speicher.** Gespeichert bräuchte die H-Matrix 1,1–2,4 GB bei 40 000 Punkten, vor allem für die dichten Blöcke
  (21 komplexe Zahlen je Paar). Daher wird jeder Block nach dem Aufbau sofort angewandt und verworfen.

Bei der ersten Messung war der Prototyp ohne `-DCBEM_USE_OPENMP` gebaut; seine eigenen Schleifen liefen dann auf einem Kern,
und die H-Matrix sah fünfmal langsamer aus (Summe der Blockzeiten = Wandzeit).

**Umsetzung** (`scattered_field_curved_hmatrix`): Variante P, η = 2, Blattgröße 64, ACA+ ohne Nachkompression, jeder Block
sofort angewandt (Summen je Thread). `exterior_near_field_curved` nutzt sie ab `hmatrix_min_points` Punkten mit
`NearFieldOptions` (Voreinstellung `curved_near_field_options()`: 4 000 Punkte, eps = 10⁻⁶; eben 2 000 Punkte und 10⁻⁴).
Die Toleranz 10⁻⁶ hält den Fehler bei 3·10⁻⁷ von max |F_s|, weit unter dem Diskretisierungsfehler (|E|² 5·10⁻⁴ bei 1 280
Elementen); 10⁻⁴ gäbe 4·10⁻⁵.

**Kosten** (`prototype/curved/near_field_cost.cpp`, Bibliothek, 12 Threads):

| Punkte | 1 280: direkt / H-Matrix | eben | 5 120: direkt / H-Matrix | eben |
|---:|---:|---:|---:|---:|
| 1 000 | 0,19 / 0,22 s | 0,10 s | 0,81 / 1,38 s | 0,40 s |
| 1 976 | 0,34 / 0,31 s | 0,18 s | 1,55 / 1,82 s | 0,76 s |
| 4 000 | 0,66 / 0,40 s | 0,20 s | 2,98 / 2,26 s | 0,56 s |
| 9 968 | 1,68 / 0,56 s | 0,46 s | 7,39 / 2,79 s | 0,97 s |
| 39 912 | 6,71 / 1,36 s | 1,44 s | 30,6 / 4,36 s | 2,90 s |

Fehler des Gesamtfelds E höchstens 3,2·10⁻⁷ (bezogen auf max |E|). Die Schwelle liegt bei 4 000 Punkten, weil die H-Matrix
mit 5 120 Elementen erst ab etwa 3 000 Punkten gewinnt. Für die Kraft über die Kugel (32 × 64 = 2 048 Punkte) bleibt die
direkte Summation schneller: 1 280 Elemente 0,36 statt 0,43 s, 5 120 Elemente 1,53 statt 2,73 s. Kräfte rechnen daher wie
bisher (gleiche Werte). Da 1 280 gekrümmte Elemente genauer sind als 5 120 ebene, kostet eine Karte mit 40 000 Punkten bei
gleicher Genauigkeit jetzt 1,4 s statt 2,9 s eben (vorher gekrümmt 7 s).

**Prüfung** (`test_curved_near_field`, Teil 7; 320 Elemente, 4 140 Punkte bis 0,01 Radien vor der Kugel): eps 10⁻⁶ Fehler
1,5·10⁻⁷, eps 10⁻⁴ 2,1·10⁻⁵, Gesamtfeld im chiralen Außenmedium (zwei Helizitäten, χ = 0,05) 1,2·10⁻⁷; Python
(`test_curved_elements`): H-Matrix gegen direkte Summation und `options`.

## Geschichtete Körper auf gekrümmten Elementen (v0.62)

**Formulierung** wie `LayeredScatteringProblem` (`results_coated.md`): je Fläche die Außenspur (24 Unbekannte je Element),
je Gebiet ein `CurvedCauchyOperator` auf seinem Rand (Flächen mit ihren eigenen Außennormalen, Vorzeichen in der Dichte),
je Fläche die Galerkin-Projektion J_G der Transmissionsabbildung für ihr Mediumpaar (`curved_transmission_blocks`, aus
`CurvedTransmissionOperator` herausgezogen, bitgleich). Für einen homogenen Körper ist das System T₁ von
`CurvedScatteringProblem`.

**Messung vor dem Einbau** (`prototype/curved/layered.cpp`, System aus den Bausteinen zusammengesetzt; Goldkern ε = −11 + 1,2i
mit Radius 1, Glasschale ε = 2,25 der Dicke d, ωa = 0,5, gegen Aden-Kerker `tools/mie_coated.py`; konstante Dichten mit
`LayeredScatteringProblem` auf denselben Ikosaedernetzen; n = 8: 1 280 Elemente je Fläche; 12 Threads):

| d | n | gekrümmt: σ_ext | Schichtwirkung gegen ohne Schicht | konstant: σ_ext | Schichtwirkung gegen neutral | Aufbau gekrümmt |
|---:|---:|---:|---:|---:|---:|---:|
| 0,2 | 4 | +8·10⁻⁶ | +1·10⁻⁵ | +1,36 % | −2,6 % | 3,1 s |
| | 8 | +2·10⁻⁶ | < 10⁻⁵ | +0,32 % | −0,71 % | 17 s |
| 0,05 | 4 | +1,7·10⁻⁴ | +0,08 % | +2,23 % | −4,5 % | 24 s |
| | 8 | +2,5·10⁻⁵ | +0,012 % | +0,90 % | −1,3 % | 44 s |
| 0,02 | 4 | +4·10⁻⁵ | +0,05 % | −0,07 % | −4,4 % | 104 s |
| | 8 | +1,6·10⁻⁵ | +0,017 % | +0,48 % | −2,2 % | 135 s |
| 0,01 | 4 | −2·10⁻⁶ | −0,004 % | −1,91 % | +2,5 % | 323 s |
| | 8 | +8·10⁻⁶ | +0,017 % | −0,04 % | −1,5 % | 443 s |

(Exakt: σ_ext = 4,3428, 2,3366, 2,0368, 1,9436; ohne Schicht 1,8536.) Lösen: 2,5 s (n = 4) bzw. 15 s (n = 8), 37–44
Iterationen wie konstant.

- **Genauigkeit:** Gekrümmte Elemente sind 100- bis 1 000-mal genauer.
- **Kein neutraler Vergleich nötig:** Der Konsistenzfehler dünner Schichten im konstanten Pfad (neutrale Schale bis 6,7 %
  daneben) ist gekrümmt höchstens 0,08 %. Die Wirkung der Schicht stimmt direkt gegen die Rechnung ohne Schicht auf etwa
  2·10⁻⁴, auch bei d/h = 0,03. Die Differenz zur neutralen Rechnung ist gekrümmt sogar ungenauer (bis 1,5 % bei n = 4,
  d = 0,01): Die beschichtete Rechnung ist genauer als die neutrale, deren Fehler sich dann nicht mehr heraushebt, sondern
  hinzukommt.
- **Zufällige Treffer im konstanten Pfad:** Einzelwerte wie −0,04 % (n = 8, d = 0,01) sind Fehlerausgleich; die Schichtwirkung
  ist dort um 1,5 % falsch.

**Nahquadratur über die Schicht.** Ohne Anpassung wächst der Aufbau mit abnehmendem d stark: n = 4, d = 0,05: 76 s; d = 0,02:
mehr als 10 Minuten (abgebrochen). Die äußere Integration verfeinert das Element, bis jedes Teilstück klein gegen seinen
Abstand zum inneren Element ist, bei parallelen Flächen also flächig auf die Größe d. Zwei Änderungen
(`CurvedNearParams::adapt_to_boundary`, Voreinstellung für geschichtete Körper `curved_layered_near_params()`; sonst aus,
bestehende Rechnungen bitgleich):
1. **Randabstand** wie im ebenen Pfad (`results_coated.md`): Liegt ein Teilstück ganz auf einer Seite der Fläche des inneren
   Elements (Höhen der Ecken über ihren Fußpunkten mit gleichem Vorzeichen), ist das Innenintegral tangential glatt und nur über
   dem Rand des inneren Elements singulär; dann zählt der Abstand zum Rand.
2. **Abstände zur gekrümmten Fläche** statt zum Sehnendreieck minus 4/3 der Wölbung: über den Fußpunkt (`tangent_at`) bzw. die
   quadratischen Randkurven (Goldener Schnitt). Bei n = 4 ist die Wölbung (etwa 0,011) so groß wie die Schicht; die
   Sehnenschätzung ließ den Abstand dann auf fast null schrumpfen.

Kosten je Element der Schale gegen alle nahen Elemente des Kerns (`prototype/curved/layered_pairs.cpp`, n = 4, etwa 50 Paare,
ein Thread; in Klammern das übereinanderliegende Paar):

| d (d/h) | Voreinstellung | nur Randabstand | Randabstand und gekrümmte Abstände |
|---|---:|---:|---:|
| 0,2 (0,56) | 0,008 s (1 ms) | 0,008 s | 0,007 s (1 ms) |
| 0,05 (0,16) | 0,18 s (86 ms) | 0,16 s | 0,135 s (54 ms) |
| 0,02 (0,07) | 2,5 s (1,65 s) | 1,7 s | 0,63 s (0,36 s) |
| 0,01 (0,03) | – | 7,0 s | 2,0 s (1,15 s) |

Die Einträge ändern sich um höchstens 1,7·10⁻⁸ (bezogen auf den größten Eintrag). Lockerere Kriterien sparen weitere 30–35 %
(äußeres Kriterium 3: Abweichung 3·10⁻⁸ bis 1·10⁻⁷; Korrektur 3: 1–7·10⁻⁶), sind aber nicht eingestellt. Die verbleibenden
Kosten wachsen etwa wie (h/d)^1,5, die Hälfte trägt das übereinanderliegende Paar: Dort verfeinert die äußere Integration ein
Band der Breite d um die Kanten, über denen die Kanten des inneren Elements liegen. Weniger würde eine zu den Kanten hin
gradierte äußere Regel kosten (offen).

**Umsetzung** (`CurvedLayeredScatteringProblem`, `curved_layered_problem.hpp`): `CurvedLayeredGeometry` mit `add_body`,
`add_layered_body`, `add_coated_body` (Parallelflächen über `offset_surface(QuadraticMesh)`: Ecken und Kantenmitten entlang der
gemittelten Normalen der angrenzenden gekrümmten Elemente; Kugel n = 4, d = 0,1: Radiusfehler der Knoten 4·10⁻⁹), auch
mehrere Körper und chirale Schichten oder Außenmedien. **Prüfung** (`test_curved_layered`, n = 4):
- ohne Schicht wie `CurvedScatteringProblem` (ein und zwei Körper, gleiche Nahquadratur: σ auf 12 Stellen gleich);
- gegen Aden-Kerker: Schale 0,2: Q_ext 8·10⁻⁶; Schale 0,1: 2·10⁻⁴ (auch über `add_coated_body`); Doppelschale (Glas bis 1,1,
  ε = 4 bis 1,2): 5·10⁻⁴; S(0) jeweils 4–6·10⁻⁴ (der Imaginärteil ist schon für die homogene Kugel bei 320 Elementen auf 1·10⁻³
  genau, bei 1 280 auf 7·10⁻⁵); optisches Theorem;
- chirale Glasschale (χ = 0,1) gegen `tools/mie_chiral_layered.py`: σ₊ 2·10⁻⁶, σ₋ 2·10⁻⁵, Zirkulardichroismus auf 1 %,
  Spiegelsymmetrie σ_s(χ) = σ_{−s}(−χ) auf 2·10⁻⁹;
- Nahquadratur mit Randabstand gegen die Voreinstellung über eine Schicht d = 0,05: 1,7·10⁻⁸.
