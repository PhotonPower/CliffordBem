# Ergebnisse: Streulöser T₁ auf der Kugel (Validierung gegen Mie)

`scatter_sphere`: ebene Welle (d = e₃, p = e₁) auf die Ikosaeder-Kugel (Radius 1), resonanzfreie
Gleichung T₁ h = h_inc mit Standardwahl von J, beide Randoperatoren E_{k₁}, E_{k₂} als H-Matrizen
(ε = 10⁻⁴), GMRES (Neustart 200, Toleranz 10⁻⁶) mit punktweiser Vorkonditionierung 2(1+J)⁻¹,
Extinktion über das optische Theorem. Mie-Referenz: `tools/mie.py`. Auswertung:
`cd tools && python3 analyze_scattering.py ../results/scatter_glass.csv ../results/scatter_gold.csv`.

Kontrolle: Auf dem Netz mit 320 Dreiecken reproduziert der C++-Kern die dichte Lösung des
Python-Prototyps auf alle angegebenen Stellen (Glas 0,202913, Gold 0,669150; Test `test_scattering`).

## Glas, ε₁ = 2,25, ωa = 1 (Mie: Q_ext = 0,215098)

| Dreiecke | Unbekannte | Q_ext | Abw. zu Mie | GMRES | H-Speicher | Aufbau | Lösen |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 280 | 10 240 | 0,211982 | 1,45 % | 10 | 109 MB | 15 s | 1 s |
| 2 880 | 23 040 | 0,213709 | 0,65 % | 10 | 338 MB | 39 s | 4 s |
| 5 120 | 40 960 | 0,214317 | 0,36 % | 10 | 679 MB | 77 s | 9 s |
| 8 000 | 64 000 | 0,214598 | 0,23 % | 10 | 1 323 MB | 129 s | 18 s |
| 11 520 | 92 160 | 0,214752 | 0,16 % | 10 | 1 918 MB | 201 s | 26 s |

Konvergenzordnung p = 1,96–2,01; Extrapolation aus den drei feinsten Netzen:
Q_∞ = 0,215111 (Abweichung 6·10⁻⁵), aus n = 12, 16, 20: 0,215096 (8·10⁻⁶).

## Gold, ε₁ = −11 + 1,2i, ωa = 0,5 (Mie: Q_ext = 0,590018)

| Dreiecke | Unbekannte | Q_ext | Abw. zu Mie | GMRES | H-Speicher | Aufbau | Lösen |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 280 | 10 240 | 0,619868 | 5,1 % | 26 | 110 MB | 14 s | 3 s |
| 2 880 | 23 040 | 0,607641 | 3,0 % | 26 | 341 MB | 38 s | 10 s |
| 5 120 | 40 960 | 0,602371 | 2,1 % | 25 | 683 MB | 76 s | 21 s |
| 8 000 | 64 000 | 0,599483 | 1,6 % | 24 | 1 330 MB | 125 s | 38 s |
| 11 520 | 92 160 | 0,597669 | 1,3 % | 24 | 1 922 MB | 197 s | 57 s |

Die beobachtete Ordnung fällt von 1,42 auf 1,29, und die Extrapolation nähert sich Mie
(3,2·10⁻³ → 2,1·10⁻³ → 1,4·10⁻³). Plausibel ist ein Übergang zur Ordnung 1 der stückweise
konstanten L²-Dichten, die bei Metallen stärker sichtbar ist als bei Glas; der Geometriefehler
(eingeschriebenes Polyeder) ist von zweiter Ordnung. Nahfeldquadratur an gemeinsamen Kanten
(Sauter-Schwab) und Ansätze höherer Ordnung sind die naheliegenden Verbesserungen.

## Bewertung

- **Stabilität (H1, glatte Ränder):** GMRES-Iterationen unabhängig von N (Glas 10, Gold 24–26)
  bis 92 160 Unbekannte.
- **Genauigkeit:** Glas konvergiert mit Ordnung 2 gegen Mie (Extrapolation auf 10⁻⁵),
  Gold mit Ordnung zwischen 1 und 1,5.
- **Kosten:** Der Speicher wird von den beiden H-Matrizen bestimmt (O(N log² N)); bei 92 160
  Unbekannten 1,9 GB, Aufbau 200 s, Lösung 26–57 s auf einem Rechenkern.
