# Resonanzfenster an Ecken: komplexe Streckung in log r (2D-Prototyp, Zwischenstand)

Prototyp: `prototype/resonance/pml2d_exact.py` (Quadrat, quasistatisch, Kollokation von T₁ mit stückweise
konstanten Multivektor-Dichten, geometrisch zu den Ecken gradiertes Netz). Beobachtungsgröße: Streufeld
F₁ an x = (5, 3) bei homogener Anregung. Aufruf: `python3 pml2d_exact.py κ θ₁,θ₂ q₁,q₂ Tiefe [r0]`
(benötigt `cl3.py`, `corner.py` aus `prototype/ap1`).

## Ansatz

Nahe jeder Ecke wird der Randparameter r komplex gestreckt,
r ↦ r̃ = r (r/r₀)^{iθ} für r < r₀, d. h. log r̃ = log r + iθ (log r − log r₀) („PML in log r“). Der Kern wird
analytisch fortgesetzt (Φ₂(z) = z/(2π z·z) mit bilinearem z·z), dl ↦ dr̃, die Normalen bleiben reell. Die
Randstücke werden exakt komplex parametrisiert (Gauß im reellen Parameter), der Hauptwert am eigenen Element
analytisch. Zwischen Punkten verschiedener Kanten entstehen keine Singularitäten (z·z = 0 verlangte gleiche
Beträge bei verschiedener Phase).

Eine Lösung ∼ r^a geht über in r^{a(1+iθ)}. Die physikalische (verlustbehaftete) Lösung enthält den
Singulärexponenten a = c − iτ mit c > 0, τ > 0; gestreckt hat er den Realteil c + τθ.

## Befunde

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

## Stand

Die Streckung in log r löst das Reflexionsproblem am Netzende, aber nicht das Genauigkeitsproblem. Das
Resonanzfenster ist damit nicht gelöst. Nächste Schritte wären: Ansatzfunktionen, die an r^{a(1+iθ)}
angepasst sind (Anreicherung mit dem gestreckten Exponenten) oder polynomiell höherer Ordnung in log r;
Galerkin statt Kollokation; eine Mellin-Analyse des diskreten gestreckten Eckproblems, um θ und die nötige
Auflösung vorherzusagen; und ein Referenzwert für das Quadrat im Fenster aus einer unabhängigen Methode.
