# Ergebnisse: Streulöser T₁ auf der Kugel (Validierung gegen Mie)

`scatter_sphere`: ebene Welle (d = e₃, p = e₁) auf die Ikosaeder-Kugel (Radius 1), resonanzfreie
Gleichung T₁ h = h_inc mit Standardwahl von J, beide Randoperatoren E_{k₁}, E_{k₂} als H-Matrizen
(ε = 10⁻⁴), GMRES (Neustart 200, Toleranz 10⁻⁶) mit punktweiser Vorkonditionierung 2(1+J)⁻¹,
Extinktion über das optische Theorem. Mie-Referenz: `tools/mie.py`. Auswertung:
`cd tools && python3 analyze_scattering.py ../results/scatter_gold.csv ../results/scatter_glass.csv`.

Ab v0.4 werden alle Dreieckspaare mit gemeinsamer Fläche, Kante oder Ecke mit
**Sauter-Schwab-Quadratur** (Ordnung 5) integriert. Vorher wurde das innere Integral analytisch
und das äußere mit einer festen Regel (112 Punkte) berechnet; diese Außenintegration konvergiert
an gemeinsamen Kanten wegen einer logarithmischen Singularität nur langsam. Einzeltest
(`test_sauter_schwab`): Für ein kantenbenachbartes Paar hatte die alte Regel 1,7 % Fehler (und mit
48² × 7 Außenpunkten noch 1,5·10⁻³), Sauter-Schwab erreicht 4·10⁻⁵ (Ordnung 4), 3·10⁻⁷ (6) und
2·10⁻⁹ (8), mit exponentieller Konvergenz.

Kontrolle: Mit der alten Nahfeldmethode (`EntryParams::sauter_schwab = false`) reproduziert der Kern
auf dem Netz mit 320 Dreiecken die dichte Lösung des Python-Prototyps auf alle Stellen (Test
`test_scattering`).

## Gold, ε₁ = −11 + 1,2i, ωa = 0,5 (Mie: Q_ext = 0,590018)

| Dreiecke | Unbekannte | Q_ext (Sauter-Schwab) | Abw. zu Mie | Q_ext (vorher) | Abw. vorher | GMRES | Aufbau |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 280 | 10 240 | 0.600720 | 1.81 % | 0.619868 | 5.1 % | 26 | 12 s |
| 2 880 | 23 040 | 0.594725 | 0.80 % | 0.607641 | 3.0 % | 26 | 31 s |
| 5 120 | 40 960 | 0.592644 | 0.45 % | 0.602371 | 2.1 % | 25 | 62 s |
| 8 000 | 64 000 | 0.591686 | 0.28 % | 0.599483 | 1.6 % | 24 | 105 s |
| 11 520 | 92 160 | 0.591164 | 0.19 % | 0.597669 | 1.3 % | 24 | 162 s |

Extrapolation Q(n) = Q_∞ − C n⁻ᵖ aus drei aufeinanderfolgenden Netzen:

| Netze | Sauter-Schwab: p | Q_∞ | Abw. zu Mie | vorher: p | Abw. zu Mie |
|---|---:|---:|---:|---:|---:|
| 1 280 / 2 880 / 5 120 | 2,02 | 0,590008 | 1,7·10⁻⁵ | 1,42 | 3,2·10⁻³ |
| 2 880 / 5 120 / 8 000 | 2,02 | 0,590005 | 2,3·10⁻⁵ | 1,35 | 2,1·10⁻³ |
| 5 120 / 8 000 / 11 520 | 1,99 | 0,589967 | 8,7·10⁻⁵ | 1,29 | 1,4·10⁻³ |

Die langsamere Konvergenz bei Gold (Ordnung ≈ 1,3–1,5, AP 2/AP 3) lag also an der
Nahfeldquadratur, nicht an der Formulierung oder den stückweise konstanten Ansätzen. Mit
Sauter-Schwab konvergiert Gold wie Glas mit Ordnung 2, und der Fehler auf dem feinsten Netz sinkt
von 1,3 % auf 0,19 %. Der Aufbau wird zugleich schneller (162 s statt 197 s bei 92 160 Unbekannten).

## Glas, ε₁ = 2,25, ωa = 1 (Mie: Q_ext = 0,215098)

| Dreiecke | Unbekannte | Q_ext | Abw. zu Mie | GMRES | Aufbau |
|---:|---:|---:|---:|---:|---:|
| 1 280 | 10 240 | 0.211964 | 1.46 % | 10 | 13 s |
| 2 880 | 23 040 | 0.213697 | 0.65 % | 10 | 33 s |
| 5 120 | 40 960 | 0.214308 | 0.37 % | 10 | 65 s |
| 8 000 | 64 000 | 0.214590 | 0.24 % | 10 | 107 s |
| 11 520 | 92 160 | 0.214745 | 0.16 % | 10 | 166 s |

Ordnung p = 1,94–2,01, Extrapolation auf 3·10⁻⁵ bis 6·10⁻⁵. Bei Glas ändert Sauter-Schwab die Werte
nur in der fünften Stelle (vorher z. B. 0,214752 bei 11 520 Dreiecken): Dort dominiert der
Geometriefehler des eingeschriebenen Polyeders.

## Bewertung

- **Stabilität (H1, glatte Ränder):** GMRES-Iterationen unabhängig von N (Glas 10, Gold 24–26)
  bis 92 160 Unbekannte.
- **Genauigkeit:** Glas und Gold konvergieren mit Ordnung 2 gegen Mie; nach Extrapolation
  Abweichungen von 10⁻⁵ bis 10⁻⁴.
- **Kosten:** bei 92 160 Unbekannten 1,9 GB für beide H-Matrizen, Aufbau 162–166 s, Lösen 19 s
  (Glas) bzw. 40 s (Gold) auf einem Rechenkern.

Die Werte der Vorversion (analytisches Innenintegral) liegen in
`results/scatter_gold_analytic.csv` und `results/scatter_glass_analytic.csv`.
