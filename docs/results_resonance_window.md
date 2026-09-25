# Resonanzfenster an Ecken: komplexe Streckung in log r (2D-Prototyp)

Prototypen: `prototype/resonance/pml2d_exact.py` (Kollokation) und `pml2d_galerkin.py` (Galerkin); Quadrat,
quasistatisch, T₁ mit stückweise konstanten Multivektor-Dichten, geometrisch zu den Ecken gradiertes Netz.
Beobachtungsgröße: Streufeld F₁ an x = (5, 3) bei homogener Anregung. Aufruf:
`python3 pml2d_galerkin.py κ θ₁,θ₂ q₁,q₂ Tiefe [r0]`, `python3 depth_study.py κ q θ₁,θ₂ Tiefe₁,Tiefe₂`.

## Ansatz

Nahe jeder Ecke wird der Randparameter r komplex gestreckt,
r ↦ r̃ = r (r/r₀)^{iθ} für r < r₀, d. h. log r̃ = log r + iθ (log r − log r₀) („PML in log r“). Der Kern wird
analytisch fortgesetzt (Φ₂(z) = z/(2π z·z) mit bilinearem z·z), dl ↦ dr̃, die Normalen bleiben reell. Die
Randstücke werden exakt komplex parametrisiert (Gauß im reellen Parameter), der Hauptwert am eigenen Element
analytisch. Zwischen Punkten verschiedener Kanten entstehen keine Singularitäten (z·z = 0 verlangte gleiche
Beträge bei verschiedener Phase).

Eine Lösung ∼ r^a geht über in r^{a(1+iθ)}. Die physikalische (verlustbehaftete) Lösung enthält den
Singulärexponenten a = c − iτ mit c > 0, τ > 0; gestreckt hat er den Realteil c + τθ.

## Befunde mit Kollokation

1. **Vorzeichen von θ.** Mit θ > 0 wird die physikalische Komponente zur Ecke hin gedämpft; mit θ < 0
   wird sie verstärkt, und das gestreckte Problem liefert eine andere, unphysikalische Lösung, die in der
   Zahl der Schichten trotzdem konvergiert. Kontrolle bei κ = −2,8 + 1,5i, wo das ungestreckte Verfahren
   konvergiert (q = 0,8): Imaginärteil von F₁ für θ = 0 / 0,5 / 1: 0,02321 / 0,02329 / 0,02333; für θ = −1:
   0,02206. Eine Konvergenz in der Schichtzahl allein ist also kein Nachweis.
2. **L²-Bedingung.** Die gestreckte Dichte ∼ r^{c+τθ−1} liegt in L², wenn c + τθ > ½, also
   θ > (½ − c)/τ. Rechter Winkel: κ = −2,8 + 1,5i: a = 0,448 − 0,261i (θ > 0,2); κ = −2,8 + 0,3i:
   a = 0,151 − 0,238i (θ > 1,5); κ = −2,8 + 0,05i: a = 0,029 − 0,212i (θ > 2,2); verlustfrei κ = −2:
   a = ±0,613i (θ > 0,8).
3. **Abschneiden.** Mit zulässigem θ spielt die Tiefe der Gradierung keine Rolle mehr (κ = −2,8 + 0,3i,
   θ = 2, q = 0,8: 0,0412209 + 0,0372974i bei Tiefe 10⁻⁴ gegen 0,0412176 + 0,0373007i bei 10⁻⁶). Der
   ursprüngliche Mechanismus des Versagens, die Reflexion der Eckwellen am inneren Netzende, ist damit
   beseitigt.
4. **Auflösung.** Im Fenster konvergiert die Diskretisierung aber zu langsam: κ = −2,8 + 0,3i, θ = 2,
   r₀ = 1: F₁ = 0,0325 + 0,0392i (q = 0,7), 0,0382 + 0,0408i (0,8), 0,0465 + 0,0409i (0,9). Der Imaginärteil
   stabilisiert sich, der Realteil nicht; bei erreichbaren Auflösungen hängt das Ergebnis noch um etwa 10 %
   von θ ab. Ursache ist vermutlich, dass stückweise konstante Kollokation das stark singuläre, in log r
   oszillierende Feld (Dichte ∼ r^{−0,85} außerhalb der Streckung) nur mit Ordnung < 1 auflöst; zudem liegt
   κ im Bereich der Plasmonresonanzen des Quadrats, wo die Beobachtungsgröße empfindlich ist.

## Galerkin statt Kollokation (`pml2d_galerkin.py`)

Galerkin mit stückweise konstanten Ansatz- und Testfunktionen, getestet mit dem Konturintegral ∫ (·) dl̃
(bilinear, ohne Konjugation, damit die analytische Fortsetzung erhalten bleibt): G = ½(M + E_g) + ½(M − E_g)J,
G h = M h_inc. Selbstterm des Doppelintegrals = 0 (antisymmetrischer Kern); benachbarte Elemente auf derselben
Kante mit analytischem Innenintegral und zum gemeinsamen Endpunkt gradierter Außenregel; Paare an derselben
Ecke über Versätze zur Ecke (vermeidet Auslöschung); der Beginn der Streckung r₀ liegt auf einem Netzknoten
(sonst Knick der Abbildung innerhalb eines Elements, sichtbar als Fehler von 10⁻³).

**Kontrollen außerhalb des Fensters** (F₁ an x = (5, 3)):

| κ | θ | q = 0,5 | q = 0,7 | q = 0,8 |
|---|---:|---:|---:|---:|
| Gold −11 + 1,2i | 0 | 0,0249364 + 0,0007079i | 0,0249605 + 0,0006831i | 0,0249666 + 0,0006766i |
| | 1 | | 0,0249598 + 0,0006789i | 0,0249661 + 0,0006747i |
| | 2 | | 0,0249637 + 0,0006730i | 0,0249674 + 0,0006722i |
| −2,8 + 1,5i | 0 | 0,0319240 + 0,0213543i | 0,0320747 + 0,0213189i | 0,0321121 + 0,0213107i |
| | 1 | | | 0,0321145 + 0,0213031i |
| | 2 | | | 0,0321242 + 0,0212952i |

Galerkin konvergiert deutlich schneller als Kollokation (Gold: Differenzen 2,4·10⁻⁵, 6·10⁻⁶), und die
Streckung ändert das Ergebnis nur im Rahmen des Diskretisierungsfehlers.

**Im Fenster, κ = −2,8 + 0,3i** (a = 0,151 − 0,238i, zulässig θ > 1,5):

| θ | q = 0,5 | 0,7 | 0,8 | 0,87 |
|---:|---:|---:|---:|---:|
| 2 | 0,059774 + 0,034205i | 0,060023 + 0,033138i | 0,060116 + 0,033132i | 0,060156 + 0,033163i |
| 3 | 0,055680 + 0,038961i | 0,059966 + 0,033176i | 0,060137 + 0,033051i | 0,060173 + 0,033124i |

Tiefe der Gradierung (q = 0,7):

| θ | 10⁻³ | 10⁻⁴ | 10⁻⁶ | 10⁻⁸ |
|---:|---:|---:|---:|---:|
| 0 | 0,062402 + 0,034503i | 0,060162 + 0,034656i | 0,059836 + 0,033143i | 0,060163 + 0,033323i |
| 2 | 0,059998 + 0,033132i | 0,060023 + 0,033138i | 0,060029 + 0,033137i | 0,060029 + 0,033137i |

Ohne Streckung schwankt das Ergebnis mit der Tiefe (Reflexion am Netzende); mit Streckung ist es ab Tiefe
10⁻⁶ auf 10⁻⁷ fest, und θ = 2 und 3 stimmen bei q = 0,87 auf 6·10⁻⁴ überein.

**Kleiner Verlust und verlustfrei:**

| κ | θ | q = 0,7 | q = 0,8 |
|---|---:|---:|---:|
| −2,8 + 0,05i (θ > 2,2) | 3 | 0,078163 + 0,040740i | 0,078658 + 0,040521i |
| | 4 | 0,077275 + 0,041471i | 0,078592 + 0,040215i |
| −2 verlustfrei (θ > 0,8) | 2 | 0,030312 + 0,071090i | 0,030522 + 0,071200i |
| | 3 | 0,030504 + 0,070817i | 0,030519 + 0,071262i |
| −2 + 0,02i | 3 | | 0,030499 + 0,070056i |

Im verlustfreien Fall hat das Streufeld einen Imaginärteil: Energie fließt in die Ecke ab („Black-Hole“-Wellen,
kontinuierliches Spektrum). Der Wert bei κ = −2 + 0,02i liegt nahe am verlustfreien, wie es das
Grenzabsorptionsprinzip verlangt. Ohne Streckung ist der verlustfreie Fall unbrauchbar (F₁ = 0,0189 + 0,0444i,
0,0078 + 0,0538i, 0,0003 + 0,0681i bei Tiefe 10⁻⁴, 10⁻⁶, 10⁻⁸).

## Stand

Mit Galerkin und komplexer Streckung in log r (θ > 0, θ > (½ − c)/τ) ist das Resonanzfenster im 2D-Prototyp
gelöst: Die Ergebnisse sind unabhängig von der Tiefe der Gradierung, konvergieren mit der Auflösung und hängen
nur noch im Rahmen des Diskretisierungsfehlers (≤ 0,3 %) von θ ab, bis hin zum verlustfreien Fall. Die
Kollokation war die zweite Ursache der schlechten Konvergenz.

Offen: (1) die Übertragung auf 3D. Der C++-Kern ist bereits Galerkin; nötig ist die Streckung im Abstand zu
Kanten und Ecken auf dem Rand. Dabei tritt eine neue Schwierigkeit auf: Mit einer reellen Koordinate entlang der
Kante kann z·z = t² + (Δr̃)² auf dem komplexen Rand verschwinden, ohne dass die Punkte zusammenfallen; das ist in
2D ausgeschlossen und muss in 3D durch die Wahl der Streckung vermieden werden. (2) Eine unabhängige
Referenz für das Quadrat im Fenster. (3) Die Wahl von θ setzt den Singulärexponenten voraus; an der Würfelecke
ist er noch nicht bekannt. Nächste Schritte wären: Ansatzfunktionen, die an r^{a(1+iθ)}
angepasst sind (Anreicherung mit dem gestreckten Exponenten) oder polynomiell höherer Ordnung in log r;
Galerkin statt Kollokation; eine Mellin-Analyse des diskreten gestreckten Eckproblems, um θ und die nötige
Auflösung vorherzusagen; und ein Referenzwert für das Quadrat im Fenster aus einer unabhängigen Methode.

## Übertragung auf 3D: Befunde zu Kanten und Ecken (`edge3d_zeros.py`, `enrich2d.py`)

**Kanten: die Spiralstreckung erzeugt Singularitäten.** An einer Kante bleibt die Koordinate t entlang der Kante
reell, gestreckt wird nur der Kantenabstand. Dann ist z·z = (t₁ − t₂)² + w_T mit dem transversalen Anteil w_T, und
z·z verschwindet für ein reelles t₁ − t₂, sobald w_T auf der negativen reellen Achse liegt. Das passiert für die
Spirale r̃ = r (r/r₀)^{iθ} bei jedem θ > 0, auf derselben Fläche wie zwischen Nachbarflächen (numerisch: 254 bis
532 Paare mit z·z = 0 auf einem Gitter von 1 500 × 1 500 Abständen). Der fortgesetzte Kern ist dort singulär;
die Konturverschiebung ist nicht zulässig. In 2D fehlt die reelle Koordinate t, deshalb tritt das dort nicht auf.

**Begrenzte Phase.** Mit r̃ = r e^{iψ(r)} und −π/2 < ψ ≤ 0 (ψ monoton) gibt es keine Nullstellen: auf derselben
Fläche, weil r cos ψ(r) streng monoton ist, zwischen Flächen im 90°-Winkel, weil r̃₁² + r̃₂² dann nicht auf der
negativen reellen Achse liegen kann (numerisch bestätigt, max |arg w_T| = 2ψ_max). Eine begrenzte Phase dämpft die
Eckwellen aber nur um den Faktor e^{τψ_max} ≤ e^{τπ/2} (τ = 0,24: 0,69) und verbessert tief in der Schicht nichts
mehr; damit ist sie als transparente Bedingung unbrauchbar.

**Ecken und Kegelspitzen: sicher.** Streckt man alle drei Koordinaten isotrop um die Spitze (ρ ↦ ρ e^{iψ(ρ)}), dann
verlangt z·z = 0 gleiche Beträge bei verschiedener Phase, was bei ψ = ψ(ρ) ausgeschlossen ist (numerisch
min |z·z|/|z|² = 0,38 über 200 000 Zufallspaare). Die 2D-Methode überträgt sich also auf Spitzen, nicht auf Kanten.

**Anreicherung als Alternative für Kanten.** Eine transparente Eckbedingung ohne komplexe Geometrie wäre auf
3D-Kanten übertragbar: Netz nur bis r_min, Eckelement mit Ansatzfunktion r^{a−1} (physikalischer Exponent), getestet
mit 1 (Petrov-Galerkin). Im 2D-Galerkin-Prototyp (Referenz aus der Streckung: 0,060029 + 0,033137i):

| Tiefe r_min | 10⁻² | 10⁻³ | 10⁻⁴ | 10⁻⁶ |
|---|---:|---:|---:|---:|
| Gold (Kontrolle; Referenz 0,02496 + 0,00068i) | 0,0249595 + 0,0006769i | 0,0249602 + 0,0006815i | 0,0249604 + 0,0006826i | |
| −2,8 + 0,3i (Fenster) | 0,05772 + 0,03997i | 0,05482 + 0,03560i | 0,05598 + 0,03233i | 0,05994 + 0,03137i |

Außerhalb des Fensters wirkt das Eckelement wie erwartet (schon ab r_min = 10⁻² richtig), im Fenster hängt das
Ergebnis weiter von der Tiefe ab: Eine Ansatzfunktion r^{a−1} mit freien Koeffizienten in allen acht Komponenten
ist keine transparente Bedingung. Nötig wäre die exakte Modenstruktur (der Nullvektor des Mellin-Symbols zu a,
also das feste Verhältnis der Komponenten), damit nur die ausfallende Welle dargestellt wird.

## Stand für 3D

- Kegelspitzen und Würfelecken: isotrope Streckung um die Spitze ist zulässig; Umsetzung im C++-Kern steht aus
  (komplexe Quadraturpunkte, analytisch fortgesetzter Kern mit √(z·z) für k ≠ 0).
- Kanten: die Streckung in log r ist nicht übertragbar. Ein aussichtsreicher Weg ist ein transparentes Kantenelement
  mit der exakten Mellin-Mode; es lässt sich im 2D-Prototyp gegen die Streckungsreferenz prüfen, bevor es in 3D
  umgesetzt wird.
