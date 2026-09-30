# Parallelisierung mit OpenMP (v0.31)

**Vorbemerkung:** Die Entwicklungsumgebung hat einen einzigen Rechenkern. Eine Beschleunigung konnte daher nicht gemessen
werden. Geprüft wurden die Korrektheit unter mehreren Threads (auf einem Kern wechseln sich die Threads ab, das deckt
Wettlaufsituationen ebenso auf) und der Anteil der Rechenzeit in parallelisierten Abschnitten, aus dem sich die erreichbare
Beschleunigung abschätzen lässt.

## Was parallel läuft

| Bereich | seit | Verfahren |
|---|---|---|
| H-Matrix-Aufbau (dichte Blöcke, ACA-Blöcke), Kerneinträge, Blockvorkonditionierer | vorher | je Block bzw. Eintrag, eigene Speicherbereiche |
| **H-Matrix-Produkt** (`KernelHMatrix::apply`) | v0.31 | je Block; mehrere Blöcke schreiben dieselben Zeilen, daher je Thread ein eigener Ausgabepuffer, am Ende aufsummiert; bei einem Thread direkt (unverändert) |
| **Nahfeld-H-Matrix** (Aufbau und Anwendung) | v0.31 | wie oben |
| **direkte Nahfeldsummation**, Nahfeld je Punkt (Windungszahl, Felder) | v0.31 | je Auswertepunkt |
| **Abstandssuche** (`distance_to_surface`) | v0.31 | Gitter seriell, Abfrage je Punkt |
| **Projektion** ebene Welle und Dipolfeld | v0.31 | je Dreieck; beim Dipol werden die adaptiv benötigten Quadraturregeln vorab seriell angelegt (die Tabelle wurde vorher in der Schleife ergänzt) |
| **Dipolraten** (Fernfeldintegration) | v0.31 | je Polarwinkel, Summen über eine OpenMP-Reduktion |
| **Zweitor** (Schleifen je Dreieck, Dirac-Operator der Flächen, dünnbesetzte Resolventen) | v0.31 | je Dreieck |

`CBEM_OMP(...)` in `core/types.hpp` wird ohne OpenMP zu nichts; `omp_threads()` gibt die Zahl der Threads.

## Korrektheit unter mehreren Threads

Streuung an einer Goldkugel (H-Matrix-Produkt), Zweitor mit Glasschale (Dreiecksschleifen, Resolventen), Nahfeldkarte mit
H-Matrix und direkt (3 200 Punkte), Dipolraten: Mit 1, 2 und 4 Threads stimmen Streuquerschnitte, Iterationszahlen,
Nahfeldsummen und Raten auf zwölf Stellen überein. Die volle Testsuite (23 Tests) besteht mit `OMP_NUM_THREADS=4`.

## Anteil paralleler Rechenzeit

Goldkugel, n = 12 (2 880 Dreiecke), ein Thread:

| Schritt | Zeit | parallel |
|---|---:|---|
| Aufbau (Kerneinträge, H-Matrix) | 30,0 s | ja |
| Lösen, 43 Iterationen (0,44 s je Iteration) | 18,9 s | fast vollständig: je Iteration ein H-Matrix-Produkt; Orthogonalisierung und punktweiser Vorkonditionierer sind vernachlässigbar |
| Nahfeld, 7 200 Punkte | 3,5 s | ja |
| Dipolraten | 11,2 s | ja |

Damit liegen deutlich über 90 % der Rechenzeit in parallelen Abschnitten; seriell bleiben Clusterbäume, Blockaufteilung,
GMRES-Vektoroperationen und die kleinen dichten Rechnungen (Nachkompression, Recycling-Unterräume). Nach Amdahl ergäbe ein
paralleler Anteil von 95–98 % mit 4 Kernen höchstens 3,5–3,8-fach und mit 8 Kernen 5,9–7,0-fach.

**Erwartung in der Praxis** (nicht gemessen): Der Aufbau ist rechenintensiv und sollte nahe an diese Grenzen kommen. Das
H-Matrix-Produkt ist speicherbandbreitenbegrenzt (jeder Block wird einmal gelesen) und dürfte auf 4 Kernen eher 2,5–3,5-fach
erreichen. Die Summation der Thread-Puffer kostet je Produkt O(p·N). Messung siehe unten.

## Messung auf einem Mehrkernrechner

Windows, `scatter_coated --n 12 --omega 0.5 --core "-11,1.2" --coat "0.05,2.25,0" --twoport` (Zweitor, Goldkern mit
Glasschale, 5 760 Dreiecke in zwei Flächen), Gesamtlaufzeit mit `Measure-Command`:

| Threads | Laufzeit | Beschleunigung |
|---:|---:|---:|
| 1 | 32,4 s | 1 |
| 8 | 6,86 s | **4,72** |

- Die Ergebnisse (σ, S(0), 30 Iterationen) stimmen auf alle ausgegebenen Stellen überein.
- Aus Amdahl, 1/((1 − f) + f/8) = 4,72, folgt ein paralleler Anteil f ≈ 0,90 und eine Effizienz von knapp 60 %. Der serielle Rest
  (etwa 3,2 s) macht bei 8 Threads fast die Hälfte der Laufzeit aus; mehr Kerne brächten ohne seine Verringerung wenig
  (16 Kerne: etwa 6,4-fach).
- Ob der Rechner 8 Kerne oder 4 Kerne mit Hyperthreading hat, ist offen; im zweiten Fall ist die Effizienz besser, als sie
  aussieht.
- Vermutete serielle Anteile beim Zweitor: Vorkonditionierer-Blöcke (16×16-Inversion je Dreieck), Aufbau der
  Flächenoperatoren (`SurfaceFV`), Transmissionsabbildungen, Clusterbäume und Blockaufteilung, Vektoroperationen von GMRES
  und der inneren Resolventen. Eine Zerlegung der Laufzeit in Phasen steht aus.

## Hinweise

- Threads über `OMP_NUM_THREADS`; Voreinstellung ist die Zahl der Kerne.
- Summationsreihenfolgen hängen von der Thread-Zahl ab; Abweichungen liegen im Bereich der Rundung (in den Proben unter 10⁻¹²).
- Verschachtelte Parallelität wird nicht genutzt (etwa `far_field` innerhalb der parallelen Richtungsschleife läuft seriell).
