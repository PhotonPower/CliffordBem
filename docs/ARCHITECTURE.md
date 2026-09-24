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
| geometry | Kugel, gradierter Würfel, Dunavant-Regeln | Gmsh-Import, anisotrope/nichtkonforme Kantennetze, gekrümmte Elemente |
| kernel | Dirac-Kern, Wilton-Integrale | chirale Kerne (k±), Helizitätsprojektoren |
| assembly | Fern-/Nahfeld, Selbstterm | Sauter-Schwab, adaptive Nahfeldkriterien für gestreckte Elemente |
| hmatrix | Clusterbaum, ACA (joint/comp), Nachkompression | Block-ACA mit 8×8-Pivots, ACA+, H-LU, complex64-Speicher, parallele Mat-Vek |
| operators | Cauchy-Operator E_k | Transmissionsoperator T₁ mit Standardwahl von J, Projektion ebener Wellen, Fernfeld/Q_ext |
| solvers | – | GMRES, punktweise 2(1+J)⁻¹, Kanten-/Eck-Blockvorkonditionierung (AP 2.5), H-LU-Blöcke |
| apps | Kompressionsbenchmark | Mie-Validierung, Würfel-Studien, Parameterstudien |
| bindings | – | pybind11-Modul für Skripting und Vergleich mit dem Prototyp |

## Bekannte Grenzen

- Die Zulässigkeitsbedingung dist > 3 h_max ist konservativ und bestraft gestreckte Elemente
  (vgl. AP 3, Abschnitt 5).
- Die Nahfeldquadratur an gemeinsamen Kanten konvergiert nur langsam (logarithmische
  Singularität der Außenintegration); für hohe Genauigkeit ist Sauter-Schwab vorzusehen.
- Resonanzfenster an Ecken (AP 1, Teil IV): kein reflexionsfreier Eckabschluss vorhanden.
