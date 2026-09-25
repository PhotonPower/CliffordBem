# Architektur und Ausbauplan

## Leitlinien

1. **Eintragsauswertung als zentrale Schnittstelle.** Alle Operatoren werden über
   `KernelEntries` ausgewertet: Fernfeld (Gauß), Nahfeld (analytisch plus Gauß), später
   Sauter-Schwab. H-Matrix, dichte Referenz und Vorkonditionierer rufen nur diese Schnittstelle
   auf. Neue Kerne (Helmholtz, Elastodynamik) kommen als weitere Implementierungen hinzu.
2. **Kernkomponenten statt 8×8-Blöcke.** Für den Cauchy-Operator werden je Dreieckspaar die vier
   Komponenten K_c = ∬ Φ_k (Skalar + Vektor) komprimiert. Die Multivektorstruktur
   E_ij = −2/√(|τ_i||τ_j|) L((s + v) n_j) steckt im Operator (`CauchyOperator`). Das spart Faktor 16
   und trennt Geometrie/Kern von der Algebra.
3. **Algebra ohne Abhängigkeiten.** `Multivector` mit Bitmasken-Blades; Linksmultiplikation als
   8×8-Matrix. Für Leistung später spezialisierte Grad-Darstellungen (Grad 1+2 für
   Maxwell-Anteile, 0+3 für Hilfsgrade) entsprechend der Blockstruktur aus AP 1.
4. **Keine externen Bibliotheken im Kern.** QR und SVD für die Nachkompression sind eigene
   Implementierungen (kleine Matrizen). LAPACK/BLAS kann später optional angebunden werden.
5. **Reproduzierbarkeit.** Jede Funktionalität hat einen Test gegen eine unabhängige Referenz
   (feine Quadratur, dichte Matrix, Python-Prototyp, Mie).

## Module und nächste Schritte

| Modul | Stand 0.1 | Nächste Schritte |
|---|---|---|
| geometry | Kugel, Würfel (gleichmäßig/gradiert), Mehrkörpernetze, Gmsh-Import (2.2/4.1), Dunavant, Sauter-Schwab | gekrümmte Elemente, nichtkonforme Kantennetze |
| kernel | Dirac-Kern, Wilton-Integrale | – |
| assembly | Fernfeld (Gauß 7×7); benachbarte Paare: Sauter-Schwab (gleichseitig) bzw. halbanalytisch gradiert (gestreckt); nahe Paare: adaptive Außenregel; Nahfeld-Cache | nichtkonforme Nachbarschaften, schnellere Nahpaare auf gestreckten Elementen |
| hmatrix | Clusterbaum, ACA (joint/comp), Nachkompression, ACA mit exakten Einträgen (gestreckte Elemente) | Block-ACA mit 8×8-Pivots, ACA+, H-LU, complex64-Speicher, parallele Mat-Vek |
| operators | Cauchy-Operator E_k, chiraler Cauchy-Operator, blockdiagonaler Mehrkörper-Innenoperator, Transmissionsoperator T₁ (Medium je Dreieck) | chirale Außenmedien, Substrate (geschichtete Außenmedien) |
| problems | ScatteringProblem (ein/mehrere Körper, chirale Medien, ebene Wellen) | Frequenzscans, Orientierungsmittelung |
| solvers | GMRES, punktweise 2(1+J)⁻¹, Kanten-/Eck-Blockvorkonditionierung (dichte LU) | H-LU für große Blöcke, überlappende Blöcke |
| sources | ebene Welle, Fernfeld, Extinktion | Dipolquellen, Nahfeldauswertung, Streumatrix |
| apps | Kompressionsbenchmark, Kugel- und Würfelstreuung | Parameterstudien, AP-4-Geometrien |
| bindings | – | pybind11-Modul für Skripting und Vergleich mit dem Prototyp |

## Bekannte Grenzen

- Auf gestreckten Elementen ist die Zahl der Nahpaare groß (D < 2,5 h_max ist für die Fernfeldregel
  nötig); die Vorberechnung dominiert dort die Aufbauzeit (L = 9: 93 s je Wellenzahl auf einem Kern).
- Sauter-Schwab konvergiert auf gestreckten Dreiecken langsam und orientierungsabhängig; es wird
  daher nur bis Seitenverhältnis 1,6 verwendet.
- Sauter-Schwab setzt konforme Netze voraus (Nachbarschaft über gemeinsame Knotenindizes). Bei
  hängenden Knoten fällt der Kern auf das analytische Innenintegral zurück.
- Resonanzfenster an Ecken (AP 1, Teil IV): kein reflexionsfreier Eckabschluss vorhanden.
