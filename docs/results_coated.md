# Ergebnisse: beschichtete Grenzflächen (dünne Schichten, Kern-Schale, Dünnschicht-Näherung)

An Materialgrenzen liegen meist dünne Schichten (Oxide, Sulfide, Hüllen von wenigen nm), die Amplitude und Phase
des gestreuten Lichts verändern. Seit v0.13 behandelt der Kern verschachtelte Gebiete exakt: jede Schichtgrenze ist
eine eigene Fläche, jedes Gebiet hat seinen eigenen Cauchy-Operator (`LayeredScatteringProblem`, Theorie: AP 1,
Nachtrag „Verschachtelte Gebiete“). Seit v0.14 gibt es zusätzlich eine Dünnschicht-Näherung erster Ordnung auf einer
Fläche (`ThinLayerScatteringProblem`, Abschnitt „Dünnschicht-Näherung“).

## Formulierung

Geometrie als Grenzflächengraph: geschlossene, nach außen orientierte Flächen Γ_s, jede zwischen einem Innengebiet
in(s) und einem Außengebiet out(s); Gebiet 0 ist der Außenraum. Unbekannt ist je Fläche die Außenspur h_s, die
Innenspur ist J_s h_s. Jedes Gebiet R liefert auf seinem Rand die Cauchy-Bedingung u − E_R(σu) = 0 (σ = +1, wenn R
innen an der Fläche liegt, sonst −1), und Zeile s des Systems ist die halbe Summe der Bedingungen beider
Nachbargebiete:

    ½ (h_s − [E_out(σu)]_s) + ½ (J_s h_s − [E_in(σu)]_s) = h_inc,s      (h_inc nur an Flächen zum Außenraum)

Für einen homogenen Körper ist das genau T₁ = E₂⁺ + E₁⁻ J. Kern-Schale-Teilchen, Mehrfachschichten, mehrere
beschichtete Körper, Einschlüsse und chirale Schichten sind derselbe Code. Die lokale Struktur an jeder Fläche ist
die von T₁; die punktweise Vorkonditionierung 2(1 + J_s)⁻¹ bleibt. Fredholm-Eigenschaft und Konsistenz sind für glatte,
getrennte Flächen bewiesen, die Eindeutigkeit ist eine (numerisch gestützte) Vermutung.

Schichtflächen erzeugt `offset_surface(m, d)`: Knotenversatz mit n_f · δ = d für alle angrenzenden Flächennormalen
(kleinste Quadrate, Pseudoinverse), also auf glatten Stücken δ = d n, an Kanten und Ecken auf Gehrung. Ein Würfel
wird zum Würfel mit Kante 2 + 2d (exakt), eine Ikosaederkugel zur Kugel mit Radius 1 + d (Fehler 4·10⁻⁴ bei n = 6),
ein abgerundeter Würfel (ρ = 0,4) mit d = 0,08 behält seine flachen Seiten exakt bei 1,08; die Volumenzunahme weicht
wegen des Polyederfehlers um 0,8 % vom exakten Wert ab. Umklappende Dreiecke und Umstülpen (Versatz durch das
Innere) werden erkannt.

```cpp
LayeredGeometry g(Medium{1.33 * 1.33});                           // Außenraum Wasser
add_coated_body(g, mesh, Medium{eps_ag}, {Coating{0.08, Medium{2.89}}});   // Schicht wächst nach außen
// add_coated_body(..., /*outward=*/false): Schicht innerhalb der Netzfläche (Oxidation verbraucht Metall)
LayeredScatteringProblem P(g, omega);
auto r = P.solve_plane_wave(d, p);                                  // r.sigma_ext, r.forward = S(0)
```

Neu ausgegeben wird die Vorwärtsamplitude S(0) = −ik p̄·E_∞(d)/|p|² in der Normierung von Bohren–Huffman
(σ_ext = 4π/k² Re S(0); Kugel: S(0) = Σ (2n+1)/2 (a_n + b_n)). arg S(0) ist die Phase des vorwärts gestreuten Lichts
relativ zur einfallenden Welle; sie bestimmt z. B. den effektiven Brechungsindex einer Teilchensuspension.

Referenz: `tools/mie_coated.py` (Aden–Kerker, beliebig viele Schichten über Riccati-Bessel-Transfer). Geprüft:
Schale = Außenmedium ergibt die Kernkugel, Schale = Kernmaterial die größere Kugel (beides auf 10⁻¹⁵),
quasistatische Polarisierbarkeit der beschichteten Kugel mit Fehler ∝ ω².

## Validierung: dicke Schale

Goldkern (ε = −11 + 1,2i, Radius 1) mit Glasschale (ε = 2,25) bis Radius 1,2, ωa = 0,5, konzentrische
Ikosaedernetze gleicher Feinheit n, H-Toleranz 10⁻⁶, GMRES 10⁻⁸.

| n | Dreiecke | Q_ext | S(0) | GMRES |
|---:|---:|---:|---|---:|
| 3 | 360 | 0,98207 | 0,08839 − 0,27516i | 45 |
| 4 | 640 | 0,97305 | 0,08757 − 0,28631i | 48 |
| 6 | 1 440 | 0,96563 | 0,08691 − 0,29468i | 40 |
| 8 | 2 560 | 0,96301 | 0,08667 − 0,29770i | 40 |
| Aden–Kerker | | 0,95996 | 0,08640 − 0,30163i | |

Q_ext und die komplexe Amplitude (also auch ihre Phase) konvergieren mit Ordnung 2 (Fehler von Q: 1,31·10⁻²,
5,67·10⁻³, 3,05·10⁻³ für n = 4, 6, 8). Ohne Schicht reproduziert `LayeredScatteringProblem` `ScatteringProblem`
bitgenau (ein und zwei Körper, `test_layered`). Eine chirale Schale erfüllt σ_s(χ) = σ_{−s}(−χ) auf 6·10⁻⁸.

## Dünne Schichten: Differenz zur neutralen Rechnung

Oxidschichten sind dünner als die Elemente (d ≪ h). Dann trägt die Rechnung einen systematischen Fehler von der
Größe des Diskretisierungsfehlers: Für d → 0 entartet das diskrete System zu den einzelnen Projektorbedingungen
E₂⁺h = h_inc und E₁⁻Jh = 0, deren Konsistenzfehler (diskret ist E² ≠ 1) die Lösung verschiebt (AP 1, Nachtrag).
Deutlich wird das an einer **neutralen Schale** aus Außenmedium, die physikalisch nichts ändert:

| n = 4 (h ≈ 0,26) | d = 0,05 | d = 0,02 | d = 0,01 |
|---|---:|---:|---:|
| σ ohne Schicht | 1,98902 | 1,98902 | 1,98902 |
| σ neutrale Schale | 1,92729 | 1,86028 | 1,81434 |

Der Fehler hängt aber kaum vom Schichtmaterial ab und hebt sich in der Differenz zur neutralen Rechnung auf
denselben Netzen heraus. Die Wirkung der Schicht ist deshalb als

    Δ = (beschichtet) − (neutral, gleiche Netze)

zu bestimmen, nicht als Differenz zur Rechnung ohne Schicht. Goldkern, Glasschale, ωa = 0,5
(`results/coated_thin.csv`, `tools/analyze_coated.py`); relativer Fehler von Δσ_ext bzw. Δ arg S(0) gegen
Aden–Kerker:

| d | n | d/h | Δσ gegen ohne Schicht | Δσ gegen neutral | Δ arg S gegen neutral |
|---:|---:|---:|---:|---:|---:|
| 0,05 | 4 | 0,19 | −17 % | −4,5 % | −21 % |
| | 6 | 0,29 | −5,7 % | −2,2 % | −10 % |
| | 8 | 0,38 | −2,6 % | −1,3 % | −6,1 % |
| | 12 | 0,57 | −0,9 % | −0,6 % | −2,8 % |
| 0,02 | 4 | 0,08 | −75 % | −4,4 % | −33 % |
| | 8 | 0,15 | −13 % | −2,2 % | −12 % |
| | 12 | 0,23 | −4,3 % | −1,2 % | −6,2 % |
| 0,01 | 4 | 0,04 | −192 % | +2,5 % | −36 % |
| | 8 | 0,08 | −38 % | −1,5 % | −16 % |
| | 12 | 0,11 | −14 % | −1,3 % | −8,9 % |

(Exakt: Δσ = 0,483, 0,183, 0,090; Δ arg S = 0,0209, 0,0082, 0,0041 rad.)

![Fehler der Schichtwirkung](fig_coated_thin.png)

- **Extinktion:** Gegen die neutrale Rechnung ist Δσ schon bei d/h ≈ 0,04 auf wenige Prozent genau, unabhängig von
  d/h. Die Differenz zur Rechnung ohne Schicht hat dort sogar das falsche Vorzeichen und wird erst für h ≲ 2d
  brauchbar.
- **Phase:** Die Phasenänderung ist die empfindlichere Größe. Gegen neutral konvergiert sie für d = 0,05 mit Ordnung
  etwa 2, für d = 0,01 langsamer (Ordnung etwa 1,3); auf dem feinsten Netz 3–9 %.
- **Arbeitsweise:** Absolutwert aus einer konvergierten Rechnung ohne Schicht (günstig, eine Fläche; auch mit HODLR
  oder Blockvorkonditionierung) plus Schichtwirkung aus beschichteter und neutraler Rechnung auf gröberen Netzen:
  σ ≈ σ₀ + (σ_beschichtet − σ_neutral), ebenso für arg S. `scatter_coated --neutral`, in `spectrum` über eine
  zweite Rechnung mit dem Außenmedium als Schichtmaterial.
- **Iterationen:** 42–44 beschichtet, 37–38 neutral, 36 ohne Schicht (punktweise Vorkonditionierung); kein Anstieg für
  dünne Schichten.

## Nahfeld paralleler Flächen

Liegen zwei Flächen im Abstand d ≪ h, sind die Paare über die Schicht hinweg nahe Paare mit analytischem
Innenintegral und adaptiver Außenregel. Die Außenregel verfeinerte bisher, bis jedes Teilstück klein gegen seinen
Abstand zum inneren Dreieck ist, also flächig auf die Größe d, mit (h/d)² Teilstücken. Liegt das Teilstück ganz auf
einer Seite der Ebene des inneren Dreiecks, sind die Innenintegrale (Φ₀ und 1/r) dort aber reell-analytisch mit
Singularitäten nur auf dem Rand des inneren Dreiecks; dann genügt der Abstand zum Rand
(`EntryParams::adapt_to_boundary`, Voreinstellung für geschichtete Körper: `layered_entry_params()`). Verfeinert
wird nur noch in einem Band der Breite d um die projizierten Kanten.

Zwei parallele Dreiecke (h ≈ 1) im Abstand d, Fehler gegen eine Referenz mit feinem Restglied (`k = 2 + 0,1i`):

| d | bisher | mit Randabstand |
|---:|---|---|
| 0,1 | 2,0·10⁻⁴ | 2,0·10⁻⁴ |
| 0,03 | 1,3·10⁻⁴ | 1,3·10⁻⁴ |
| 0,01 | 1,2·10⁻⁴ | 1,5·10⁻⁴ |

Der verbleibende Fehler kommt in beiden Fällen von der 7-Punkt-Regel für das glatte Restglied; die Rechenzeit sinkt
je Paar um den Faktor 3–4 (d = 0,01: 0,069 s → 0,019 s). Im Gesamtproblem ist der Gewinn kleiner, weil die meisten
Nahpaare keine eng benachbarten Paare über die Schicht hinweg sind (Goldkern mit Glasschale, n = 8, 2 560 Dreiecke,
`results/coated_nearrule.csv`):

| d | Nahfeld bisher | mit Randabstand | Ergebnis |
|---:|---:|---:|---|
| 0,05 (d/h = 0,38) | 5,9 s | 6,1 s | identisch auf 6 Stellen |
| 0,01 (d/h = 0,08) | 18,4 s | 14,9 s | identisch auf 6 Stellen |

Die Nahfeldzeit wächst mit abnehmendem d (d = 0,05 → 0,01: etwa Faktor 2,5); gegen die Gesamtzeit (Aufbau der
H-Matrizen, GMRES) bleibt sie aber untergeordnet. Für bestehende Rechnungen mit einer Fläche bleibt die
Voreinstellung unverändert (`adapt_to_boundary = false`).

Parallelflächen statt konzentrischer Kugeln (`scatter_coated --offset`, d = 0,05, n = 8): Q_ext = 0,681306 statt
0,680712 (0,09 %), entsprechend dem Radiusfehler der Parallelfläche von 4·10⁻⁴.

## Anwendung: Silberteilchen mit 2 nm Oxidschicht

### Silberkugel, Radius 20 nm, in Wasser

Silber nach Johnson–Christy, Hintergrund n = 1,33, Schicht 2 nm mit n = 1,7 (Modellwert für eine dünne Oxid-/
Sulfidschicht; eigene Daten über `--coating "2:Datei.yml"`), nach außen gewachsen. `spectrum --sphere 8` (1 280
Dreiecke je Fläche, h ≈ 2,6 nm, d/h ≈ 0,8), H-Toleranz 10⁻⁴; neutral mit Wasser als Schichtmaterial. Referenz:
Aden–Kerker mit denselben interpolierten Daten (`results/agsphere20_*.csv`, `tools/plot_coated.py`).

![Silberkugel mit Oxid](fig_coated_agsphere.png)

| λ (nm) | σ Mie ohne | σ Aden–Kerker | σ BEM korrigiert | Δσ exakt | Δσ BEM | Δ arg S exakt | Δ arg S BEM |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 360 | 1 617 | 1 096 | 1 100 | −520 | −524 | 0,056 | 0,058 |
| 370 | 3 228 | 1 739 | 1 751 | −1 489 | −1 400 | 0,172 | 0,157 |
| 380 | 5 364 | 3 747 | 3 591 | −1 617 | −1 835 | 0,040 | 0,068 |
| 390 | 14 458 | 5 928 | 5 961 | −8 530 | −8 471 | 0,342 | 0,350 |
| 400 | 25 689 | 15 740 | 15 528 | −9 949 | −9 251 | 0,880 | 0,875 |
| 410 | 10 442 | 27 135 | 26 309 | +16 693 | +15 991 | 0,712 | 0,709 |
| 420 | 4 390 | 12 203 | 12 197 | +7 813 | +7 776 | 0,292 | 0,297 |
| 430 | 2 272 | 5 145 | 5 232 | +2 873 | +2 916 | 0,137 | 0,140 |

(σ in nm², Phase in rad; „korrigiert“ = BEM ohne Schicht + (beschichtet − neutral).)

- Die 2-nm-Schicht verschiebt die Dipolresonanz um etwa 10 nm nach Rot (Maximum von knapp 400 auf etwa 408 nm) und
  hebt sie leicht an; die Phase des vorwärts gestreuten Lichts ändert sich um bis zu 0,9 rad (bei 400 nm). Die
  Phasenänderung ist auf der kurzwelligen Flanke der Resonanz am größten, weil dort die Resonanz über die
  Wellenlänge „hinwegwandert“.
- Die BEM trifft die Phasenänderung an der Resonanz auf 0,5 % (400, 410 nm) bis 2 % (390, 420 nm) und Δσ auf 1–7 %.
  Im Bereich des Quadrupols (370–380 nm) sind die Abweichungen größer (Δσ 6–13 %, Phase 0,015–0,028 rad absolut):
  höhere Multipole brauchen feinere Netze.
- Hier ist d/h ≈ 0,8; die Differenz zur Rechnung ohne Schicht stimmt deshalb mit der zur neutralen Rechnung auf
  10⁻³ überein. Bei gröberen Netzen oder dünneren Schichten gilt das nicht mehr (siehe oben).
- Kosten je Wellenlänge (ein Kern): ohne Schicht 11 s, mit Schicht 30 s (2 560 Dreiecke, 22–61 Iterationen).

### Silberwürfel, Kante 50 nm, Rundungsradius 10 nm, in Wasser

Derselbe Würfel wie in `results_roundcube.md` (Gmsh, ρ = 0,4 Einheiten, c = 12, 1 208 Dreiecke, Einheit 25 nm),
2 nm Oxid (n = 1,7) nach außen über `offset_surface` (auf den flachen Seiten h ≈ 7,5 nm, an den Rundungen ≈ 5 nm,
d/h ≈ 0,3–0,4). `spectrum --coating`, H-Toleranz 10⁻³; ohne Schicht die Werte aus v0.12
(`results/agcube_rho0.4.csv`), Schichtwirkung als Differenz zur neutralen Rechnung
(`results/agcube_rho0.4_ox2*.csv`). Eine Referenz gibt es hier nicht.

![Silberwürfel mit Oxid](fig_coated_agcube.png)

| λ (nm) | 340 | 360 | 380 | 400 | 420 | 440 | 460 | 480 | 500 | 520 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| σ ohne Schicht | 3 704 | 5 343 | 5 165 | 9 028 | 20 254 | **40 030** | 31 609 | 14 889 | 7 842 | 4 646 |
| σ mit 2 nm Oxid | 2 895 | 5 674 | 7 815 | 6 650 | 13 900 | 26 360 | **43 039** | 26 535 | 12 940 | 7 007 |
| Δ arg S(0) (rad) | +0,07 | −0,27 | −0,13 | −0,11 | +0,23 | +0,40 | **+0,50** | +0,25 | +0,12 | +0,07 |

- Die Hauptresonanz verschiebt sich um etwa 20 nm nach Rot (Maximum von etwa 440 auf etwa 460 nm, Raster 20 nm),
  doppelt so stark wie bei der Kugel mit derselben Schicht: Das Feld ist an Kanten und Ecken konzentriert, wo die
  Schicht den größten Anteil des Nahfelds füllt. Die Nebenstruktur bei 380 nm wird stärker.
- Die Phase des vorwärts gestreuten Lichts ändert sich um bis zu 0,5 rad; unterhalb von 410 nm wechselt das Vorzeichen.
  Die Verschiebung ist so groß wie die zwischen den Rundungsradien 10 und 7,5 nm (`results_roundcube.md`): Schicht
  und Kantenform sind bei Nanowürfeln gleich wichtige Parameter und aus der Streuung allein kaum zu trennen.
- Die neutrale Rechnung weicht hier nur um 0,5–3,5 % von der Rechnung ohne Schicht ab (d/h ≈ 0,3); die Korrektur
  ist klein gegen die Schichtwirkung, aber an der Resonanzflanke nicht vernachlässigbar.
- Kosten je Wellenlänge: 26–80 s mit Schicht (2 416 Dreiecke, 56–115 Iterationen), ohne Schicht 8–12 s.


## Dünnschicht-Näherung erster Ordnung (eine Fläche)

Idee: eine Greensche Funktion, die die Schichtfolge enthält. Exakt gibt es sie nur für ebene Schichtungen
(Sommerfeld-Integrale); für eine dünne Schicht genügt die lokal ebene Näherung, Krümmung geht erst in O(d²) ein. In
erster Ordnung in d ist die Wirkung der Schicht ein Sprung der Felder an der Referenzfläche Γ (X = verdrängtes
Medium):

    [E_t] = d ( ∇_Γ(D_n (1/ε_c − 1/ε_X)) − iω (μ_c − μ_X) n × H_t ),   [D_n] = −d (ε_c − ε_X) ∇_Γ·E_t
    [H_t] = d ( ∇_Γ(B_n (1/μ_c − 1/μ_X)) + iω (ε_c − ε_X) n × E_t ),   [B_n] = −d (μ_c − μ_X) ∇_Γ·H_t

In Dirac-Form ist das der Term erster Ordnung des Transfers exp(d n(ik_c − ∇_Γ)) über die Schicht. Die Innenspur wird
J_eff h = J h + L_d h, und T_eff = E₂⁺ + E₁⁻ J_eff lebt auf **einer** Fläche (`ThinLayerScatteringProblem`). E₁(L_d h)
ist das Cauchy-Integral der Schicht-Greenschen Funktion erster Ordnung, nach partieller Integration der Ableitungen
vom Kern auf die Dichte.

Diese Reihenfolge ist wesentlich: Als Kern angewandt wäre die Schichtkorrektur hypersingulär (∼ d/ρ³). Auf stückweise
konstanten Dichten ergäbe das an jeder Elementkante Beiträge ∼ d log(h/d), also einen inkonsistenten Operator. Auf der
Dichte genügen dagegen lokale Flächenableitungen über die drei Kantennachbarn (`SurfaceFV`), und alle H-Matrizen
bleiben die der unbeschichteten Rechnung:

- **Gradient:** Kleinste-Quadrate-Fit in der Tangentialebene (exakt für lineare Funktionen). Der naheliegende
  Green-Gauß-Gradient mit Kantenmitteln ist auf Ikosaedernetzen inkonsistent: Sein maximaler Fehler bleibt bei
  0,18 → 0,15 → 0,14 → 0,13 (n = 6 … 48) stehen, der Kleinste-Quadrate-Gradient fällt mit Ordnung 1
  (0,023 → 0,011 → 0,006 → 0,003).
- **Divergenz:** Flussform mit Kantenmitteln; konservativ (keine Nettoladung der Schicht), Ordnung 1
  (max. Fehler 0,12 → 0,066 → 0,034 → 0,017 für div_Γ der tangentialen Projektion eines konstanten Feldes).

Die Lage der Referenzfläche ist wählbar (Anteil f der Schicht innerhalb, `--thin-ref f`); beide Anteile addieren sich in
erster Ordnung. Nutzung:

```bash
scatter_coated --n 8 --omega 0.5 --core -11,1.2 --coat 0.01,2.25,0 --thin-ref 0.5 --bare     # Kugel, Wirkung = Differenz zu bare
spectrum --sphere 8 --unit 20 --materials Ag --nbg 1.33 --coating "2:2.89,0" --thin 0.5 ...
```

### Genauigkeit

Goldkern, Glasschale, ωa = 0,5 (`results/thin_layer.csv`, `tools/analyze_thin.py`); Schichtwirkung als Differenz zur
Rechnung ohne Schicht auf demselben Netz (keine neutrale Rechnung nötig), relative Fehler gegen Aden–Kerker, Δσ / Δ arg S:

| d | n = 4 | n = 8 | n = 16 | Referenz Mitte, n = 12 | zwei Flächen, n = 12 |
|---:|---:|---:|---:|---:|---:|
| 0,005 | +4,0 / +4,4 % | +1,6 / +2,0 % | +1,0 / +1,0 % | – | – |
| 0,01 | +4,9 / +5,3 % | +2,5 / +2,6 % | +1,9 / +1,9 % | – | −1,3 / −8,9 % |
| 0,02 | +6,8 / +7,5 % | +4,4 / +4,7 % | +3,8 / +4,1 % | +2,6 / +2,6 % | −1,2 / −6,2 % |
| 0,05 | +13,3 / +14,5 % | +10,9 / +11,5 % | +10,4 / +10,8 % | +7,0 / +6,9 % | −0,6 / −2,8 % |
| 0,1 | +28,0 / +28,8 % | +25,8 / +25,6 % | – | +17,1 / +16,0 % | – |

(Referenzfläche an der Kernoberfläche, außer „Referenz Mitte“; exakt: Δσ = 0,0446, 0,0900, 0,1833, 0,4830, 1,0535.)

![Dünnschicht-Näherung](fig_coated_thinlayer.png)

- **Diskretisierung:** Der Diskretisierungsfehler ist klein und fällt mit dem Netz; bei n = 8 liegt er unter 1 %,
  mit Referenz in der Schichtmitte ist die Rechnung schon bei n = 4–6 konvergiert. Die Iterationen bleiben für d < h bei
  35–60 (wie ohne Schicht) und steigen für d ≳ h (d = 0,1, n = 12: 150), weil die Norm von L_d wie d/h wächst.
- **Modellfehler:** Übrig bleibt der Abbruchfehler erster Ordnung, etwa 2 d/a relativ zur Schichtwirkung
  (Referenz an der Kernoberfläche) bzw. 1,3–1,7 d/a (Schichtmitte). Die Rechnung ist dabei nicht streng linear in d
  (das System wird selbstkonsistent gelöst) und liegt näher an der exakten als an der linearisierten Lösung.
- **Vergleich:** Für die Extinktionsänderung ist die exakte Zwei-Flächen-Rechnung mit neutraler Referenz auf feinen
  Netzen genauer. Für die Phasenänderung ist die Näherung unterhalb von d/a ≈ 0,03 genauer (d = 0,01: 1,9 % gegen
  8,9 %), und sie kostet etwa ein Fünftel (n = 12: rund 50 s gegen rund 300 s für beschichtete und neutrale Rechnung).
- **Kontrollen** (n = 8, d = 0,02): Schicht aus Kernmaterial (entspricht der größeren Kugel) −1,6 / −1,8 %, Schicht
  nach innen (verdrängt Gold) +2,7 / −0,3 %, verlustbehaftete Schicht (ε = 4 + i) +4,9 / +5,6 %.

### Silberkugel mit 2 nm Oxid (d/a = 0,1)

Derselbe Fall wie oben mit `spectrum --thin 0.5` (`results/agsphere20_ox2_thin.csv`, grüne Rauten in der Abbildung der
Silberkugel). Das Maximum bei 410 nm trifft Aden–Kerker (27 116 gegen 27 135 nm²), an den Flanken überschätzt die
Näherung die Verschiebung (400 nm: 13 453 gegen 15 740; 420 nm: 14 229 gegen 12 203), die Phasenänderung liegt an der
Resonanz 8–20 % zu hoch. Das entspricht dem Modellfehler bei d/a = 0,1. Für diesen Fall ist die Zwei-Flächen-Rechnung
(0,5–2 % in der Phase) die richtige Wahl.

### Einordnung

| | exakt (zwei Flächen) | Dünnschicht-Näherung |
|---|---|---|
| Gültigkeit | beliebige Dicke, Mehrfachschichten, chiral | d ≪ a (Fehler ≈ 1,5 d/a), d ≲ h, achiral, ein Körper |
| dünne Schichten (d ≪ h) | Differenz zur neutralen Rechnung nötig | direkt |
| Kosten | zwei Flächen, dazu die neutrale Rechnung | wie ohne Schicht |
| stärkste Seite | Extinktion, dicke Schichten | Phase bei d/a ≲ 0,03 |

Nächster Schritt wäre die zweite Ordnung in d: Terme d²∇_Γ² und die Formoperatoren der Fläche (Krümmung). Damit sollte
der Modellfehler auf O((d/a)²) fallen und die Näherung auch für 2 nm auf 20-nm-Teilchen reichen.

## Bewertung und Grenzen

- Die exakte Schichtformulierung ist allgemein (beliebige Dicke, Mehrfachschichten, chiral, mehrere Körper) und
  validiert. Für dicke Schalen (Kern-Schale-Teilchen) genügen die Absolutwerte.
- Für Schichten dünner als die Elemente ist die Differenz zur neutralen Rechnung nötig. Ohne sie wären Netze mit
  h ≲ d nötig, bei 1–2 nm Oxid auf 50-nm-Teilchen also 10⁴ bis 10⁵ Dreiecke je Fläche.
- Die Kosten verdoppeln sich je Schicht etwa (eine Fläche mehr, zwei Gebietsoperatoren je Fläche). Die Block- und
  HODLR-Vorkonditionierer sind für geschichtete Körper noch nicht umgesetzt; punktweise Vorkonditionierung genügt
  bisher (40–60 Iterationen).
- `offset_surface` kann an konkaven Stellen oder bei Versatz größer als der Krümmungsradius Selbstdurchdringungen
  erzeugen, die nur teilweise erkannt werden (umklappende Dreiecke, Umstülpen).
- Eindeutigkeit für verschachtelte Gebiete ist nicht bewiesen (AP 1, Vermutung).
- Für d ≪ h und d ≪ a gibt es die Dünnschicht-Näherung erster Ordnung (Abschnitt „Dünnschicht-Näherung“).
