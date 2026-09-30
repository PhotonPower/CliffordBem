# Ergebnisse: Asymptotik der Kompression (Hypothese H3)

Alle Rechnungen: `bench_compression`, ε = 10⁻⁴, gemeinsame ACA (`joint`), η = 1, Blattgröße 32,
Zulässigkeit zusätzlich dist > 3 h_max, ein Rechenkern (g++ 13, -O3). „Fehler“ ist der relative
Fehler des Matrix-Vektor-Produkts gegenüber exakt assemblierten Zeilen. Rohdaten in `results/`.

## 1. Glatte Geometrie: Kugel, k = 1,5

| Dreiecke | Unbekannte | gesamt | Niedrigrang | 10³·MB/(N log²N) | Aufbau | Fehler |
|---:|---:|---:|---:|---:|---:|---:|
| 1 280 | 10 240 | 55 MB | 25 MB | 0,83 | 6 s | 1,1·10⁻⁵ |
| 2 880 | 23 040 | 169 MB | 91 MB | 0,93 | 16 s | 9,0·10⁻⁶ |
| 5 120 | 40 960 | 340 MB | 221 MB | 0,91 | 33 s | 9,8·10⁻⁶ |
| 8 000 | 64 000 | 662 MB | 342 MB | 1,02 | 52 s | 7,0·10⁻⁶ |
| 11 520 | 92 160 | 959 MB | 644 MB | 0,95 | 80 s | 6,9·10⁻⁶ |
| 15 680 | 125 440 | 1 474 MB | 866 MB | 1,01 | 113 s | 5,9·10⁻⁶ |
| 20 480 | 163 840 | 1 881 MB | 1 388 MB | 0,93 | 156 s | 8,2·10⁻⁶ |

## 2. Geometrie mit Kanten und Ecken: gleichmäßiger Würfel, k = 1,5

| Dreiecke | Unbekannte | gesamt | Niedrigrang | 10³·MB/(N log²N) | Aufbau | Fehler |
|---:|---:|---:|---:|---:|---:|---:|
| 1 200 | 9 600 | 51 MB | 22 MB | 0,84 | 7 s | 1,6·10⁻⁵ |
| 3 072 | 24 576 | 170 MB | 93 MB | 0,86 | 21 s | 1,5·10⁻⁵ |
| 6 912 | 55 296 | 488 MB | 290 MB | 0,90 | 53 s | 1,4·10⁻⁵ |
| 12 288 | 98 304 | 960 MB | 661 MB | 0,88 | 100 s | 1,4·10⁻⁵ |
| 19 200 | 153 600 | 1 651 MB | 1 187 MB | 0,88 | 171 s | 9,6·10⁻⁶ |

Kanten und Ecken ändern das Wachstum nicht: Der Quotient ist innerhalb von ±4 % konstant,
die lokalen Exponenten liegen bei 1,18–1,30, praktisch wie bei der Kugel.

## 3. Frequenzabhängigkeit

Feste Kugel mit 11 520 Dreiecken (Durchmesser D = 2):

| k | kD | Niedrigrang | gesamt | mittl. / max. Rang | Aufbau | Fehler |
|---:|---:|---:|---:|---:|---:|---:|
| 1,5 | 3 | 644 MB | 959 MB | 7,0 / – | 80 s | 6,9·10⁻⁶ |
| 3 | 6 | 661 MB | 976 MB | 7,0 / 10 | 88 s | 1,3·10⁻⁵ |
| 6 | 12 | 716 MB | 1 031 MB | 7,2 / 11 | 95 s | 2,2·10⁻⁵ |
| 12 | 24 | 812 MB | 1 128 MB | 7,8 / 15 | 108 s | 3,7·10⁻⁵ |
| 24 | 48 | 1 042 MB | 1 357 MB | 9,0 / 16 | 135 s | 5,8·10⁻⁵ |

Feste Elementzahl je Wellenlänge (k = n/4, etwa 10 Elemente je Wellenlänge):

| Dreiecke | k | gesamt | Niedrigrang | mittl. Rang | Aufbau | Fehler |
|---:|---:|---:|---:|---:|---:|---:|
| 1 280 | 2 | 55 MB | 25 MB | 7,1 | 6 s | 1,3·10⁻⁵ |
| 5 120 | 4 | 354 MB | 236 MB | 7,2 | 34 s | 2,1·10⁻⁵ |
| 11 520 | 6 | 1 031 MB | 716 MB | 7,2 | 95 s | 2,2·10⁻⁵ |
| 20 480 | 8 | 2 075 MB | 1 582 MB | 7,1 | 192 s | 2,8·10⁻⁵ |

Im Bereich kD ≤ 16 bei fester Auflösung je Wellenlänge kostet die Frequenz nur etwa 10 %
gegenüber k = 1,5 (2 075 gegenüber 1 881 MB bei 20 480 Dreiecken); Exponent 1,31 gegenüber 1,28.
Bei festem Netz wächst der Niedrigrang-Speicher bis kD = 48 um 62 %. Der ACA-Fehler steigt mit
k von 7·10⁻⁶ auf 6·10⁻⁵, bleibt also unter der Schranke 10 ε, zeigt aber, dass das heuristische
Abbruchkriterium bei oszillierenden Kernen weniger scharf ist (Kandidat für ACA+).

## 4. Zu den Kanten gradierter Würfel (Tensornetz)

| Stufe L | Dreiecke | gesamt | dicht | Niedrigrang | Aufbau |
|---:|---:|---:|---:|---:|---:|
| 3 | 768 | 35 MB | 34 MB | 0,6 MB | 8 s |
| 5 | 1 728 | 129 MB | 114 MB | 14 MB | 30 s |
| 7 | 3 072 | 293 MB | 246 MB | 48 MB | 81 s |

Bei gleicher Dreieckszahl (3 072) braucht der gradierte Würfel 293 MB, der gleichmäßige 170 MB;
84 % davon sind dichte Nahblöcke. Ursache: Die gestreckten Elemente an den Kanten haben eine
lange Kante (bis 0,5), und sowohl das Nahfeldkriterium als auch die Zulässigkeit messen mit
h_max. Das ist ein Problem der Quadratur- und Zulässigkeitskriterien, nicht der
Niedrigrangstruktur. Abhilfe: Fernfeldquadratur mit Unterteilung gestreckter Elemente entlang
der langen Achse, dann Kriterien mit dem tatsächlichen Abstand zur Elementausdehnung.

## Bewertung

- **H3, zweiter Teil:** Für glatte und kantige Geometrien bis etwa 160 000 Unbekannte wächst der
  Speicher wie O(N log² N), die Aufbauzeit mit Exponent 1,1–1,2. Bei konstanter Auflösung je
  Wellenlänge bleibt das bis kD ≈ 16 erhalten. Damit ist die Hypothese in diesem Bereich
  numerisch gestützt.
- **Grenzen:** Für stark gradierte Netze ist die Kompression derzeit durch die Nahfeld- und
  Zulässigkeitskriterien begrenzt. Für elektrisch große Streuer (kD ≫ 50) ist das Wachstum der
  Ränge noch nicht untersucht; dort sind H²- bzw. richtungsabhängige Verfahren zu erwarten
  (vgl. Rückfallebene in AP 3 des Antrags).

## ACA+ (v0.30)

Die teilpivotisierte ACA wählt das nächste Pivot nur aus der zuletzt berechneten Zeile bzw. Spalte und bricht ab, sobald das
neue Kreuz klein ist. Sie kann Teile eines Blocks nie sehen, etwa Zeilen und Spalten, auf denen die bisherigen Kreuze
verschwinden, und bricht dann unbemerkt mit zu kleinem Rang ab. ACA+ (Grasedyck 2005; `aca_plus` in `aca.hpp`) führt zusätzlich
eine Referenzzeile und eine Referenzspalte mit, deren Reste nach jedem Schritt aktualisiert werden: Das Pivot kommt aus dem
größeren Rest der beiden Referenzen, eine Referenz wird ersetzt, sobald sie Pivot war oder ihr Rest verschwindet (relativ zum
größten bisher gesehenen Eintrag -- nach exakter Erfassung ist der Rest Rundungsrauschen, nicht null), und abgebrochen wird erst,
wenn das neue Kreuz und die auf den Block hochgerechneten Referenzreste unter eps·‖S‖ liegen. Seit v0.30 Voreinstellung
(`HMatrixParams::aca_plus`, Modi Joint und Separate sowie das Nahfeld); die teilpivotisierte ACA bleibt mit `aca_plus = false`.

| Fall | teilpivotisierte ACA | ACA+ |
|---|---|---|
| A = u₁v₁ᵀ + 0,3 u₂v₂ᵀ, disjunkte Träger (`test_aca`) | Rang 2, Fehler 1,6·10⁻¹ | Rang 2, Fehler 10⁻¹⁶ |
| Kugel n = 8, eps = 10⁻⁴, k = 0,5 | Fehler 5,5·10⁻⁶, Rang 6,9, 2,3 s | 4,8·10⁻⁶, 6,9, 2,4 s |
| Würfel n = 10 (ebene Flächen), eps = 10⁻⁴ | 7,7·10⁻⁶, 8,2, 2,0 s | 6,6·10⁻⁶, 8,3, 2,1 s |
| Kugel n = 12, eps = 10⁻⁶ | 5,6·10⁻⁸, 12,1, 11,9 s | 5,0·10⁻⁸, 12,1, 12,3 s |
| Würfel n = 16, eps = 10⁻⁶ | 8,2·10⁻⁸, 14,9, 12,9 s | 7,3·10⁻⁸, 14,9, 13,1 s |
| Kugel n = 12, eps = 10⁻⁴, k = 3 | 1,7·10⁻⁵, 7,1, 9,2 s | 1,6·10⁻⁵, 7,1, 9,8 s |
| Nahfeld, 7 200 Punkte × 2 560 Dreiecke, eps = 10⁻⁴ | max. Fehler 8,2·10⁻⁶, 2,8 s | 6,0·10⁻⁶, 2,9 s |

(Fehler der Anwendung gegen die dichte Matrix, η = 0; Speicher in allen Fällen gleich.)

- Bei den bisherigen BEM-Matrizen versagt die teilpivotisierte ACA nicht, auch nicht auf ebenen Flächen, wo einzelne der vier
  gestapelten Kernkomponenten verschwinden. ACA+ liefert dieselben Ränge und denselben Speicher, ist etwa 10–30 % genauer und
  braucht 2–6 % mehr Aufbauzeit.
- ACA+ ist damit eine Versicherung gegen unbemerkt falsch approximierte Blöcke bei Geometrien, die hier noch nicht vorkamen (dünne
  Strukturen, sehr unterschiedliche Clustergrößen, stark variierende Kernkomponenten), für wenige Prozent Aufbauzeit.
