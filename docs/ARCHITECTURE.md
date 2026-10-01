# Architektur und Ausbauplan (Stand 0.35)

## Leitlinien

1. **Eintragsauswertung als zentrale Schnittstelle.** Alle Operatoren werden über
   `KernelEntries` ausgewertet: Fernfeld (Gauß), benachbarte Paare (Sauter-Schwab bzw. halbanalytisch),
   nahe Paare (analytisches Innenintegral, adaptive Außenregel), einmal vorberechnet im Nahfeld-Cache.
   H-Matrix, dichte Blöcke, Systemeinträge und Vorkonditionierer rufen nur diese Schnittstelle auf. Neue Kerne (Helmholtz, Elastodynamik) kommen als weitere Implementierungen hinzu.
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

| Modul | Stand 0.12 | Nächste Schritte |
|---|---|---|
| core | Grundtypen, Materialmodelle (konstant, n/k-Tabellen) | Drude-Lorentz-Fits, Größenkorrektur der Dämpfung |
| clifford | Multivektoren, geometrisches Produkt, Inverse (Nullteiler-Erkennung), Linksmultiplikation | spezialisierte Grad-Darstellungen |
| geometry | Kugel, Würfel (gleichmäßig/gradiert), Mehrkörpernetze, Parallelflächen (`offset_surface`, auf Gehrung, Prüfung auf Faltung/Durchdringung), Abstand Punkt–Fläche mit Gittersuche, Windungszahl, konform verdichtete Ikosaederkugel, Prüfung auf sich berührende oder durchdringende Körper, Gmsh-Import (2.2/4.1), Dunavant, Sauter-Schwab | gekrümmte Elemente, nichtkonforme Kantennetze, Selbstdurchdringung einzelner Gmsh-Netze |
| kernel | Dirac-Kern, Wilton-Integrale (asinh-Form) | analytisch fortgesetzter Kern für komplexe Punkte (Streckung um Spitzen) |
| assembly | Fernfeld (Gauß 7×7); benachbarte Paare: Sauter-Schwab (gleichseitig) bzw. halbanalytisch gradiert (gestreckt); nahe Paare: adaptive Außenregel, optional mit Randabstand (parallele Flächen); Nahfeld-Cache | nichtkonforme Nachbarschaften, schnellere Nahpaare auf gestreckten Elementen |
| hmatrix | Clusterbaum (Netze und beliebige Punkte), ACA (gemeinsam, komponentenweise, Multivektor-Pivots), ACA+ (Voreinstellung), Nachkompression, ACA mit exakten Einträgen | complex64-Speicher |
| operators | Cauchy-Operator, chiraler Innenoperator, blockdiagonaler Mehrkörper-Innenoperator, T₁ mit Medium je Dreieck, geschichtete Transmission (Gebietsoperatoren je Gebiet), dichte Blöcke | Substrate (geschichtete Außenmedien, Green-Funktion der Schichtung) |
| solvers | GMRES, punktweise 2(1+J)⁻¹, Blockvorkonditionierung (Kanten/Ecken, Cluster; mehrere Körper, chiral), HODLR-Faktorisierung, Krylov-Recycling (rohe Krylov-Vektoren und GCRO-DR mit harmonischen Ritz-Vektoren), Eigenlöser für kleine komplexe Matrizen | H-LU mit starker Zulässigkeit, Deflation resonanter Moden, parallele GMRES-Vektoroperationen |
| sources | ebene Welle (linear/zirkular), Fernfeld, Extinktion, Vorwärtsamplitude S(0) (Betrag, Phase), Lebedev-Richtungen; chirales Außenmedium: Helizitätswelle, Extinktion je Kanal (`chiral_incidence`); Nahfeld im Außenraum mit halbanalytischer Nahquadratur, Feldverstärkung, optische Chiralität, rechteckige H-Matrix `NearFieldOperator` (`near_field`); elektrischer und magnetischer Dipol als Quelle, Zerfallsraten, Fluoreszenzverstärkung, zirkular polarisierte Lumineszenz (`dipole`); optische Kräfte über den Spannungstensor und in Dipolnäherung, auch enantioselektiv (`optical_force`) | Dipole in Schichten und in chiralen Medien, Felder in Kernen und Schichten, Streumatrix, Spannungstensor im chiralen Medium, Kraftkarten |
| problems | ScatteringProblem (ein/mehrere Körper, chirale Medien, Hintergrundmedium, Systemeinträge); LayeredScatteringProblem (Grenzflächengraph: Kern-Schale, Mehrfachschichten, Einschlüsse, chirale Schichten); ThinLayerScatteringProblem (eine Fläche je Körper, mehrere Körper; Dirac-Form 2. Ordnung mit Formoperator, glatten Normalen und quadratischer Anpassung der zweiten Ableitungen über die Knotennachbarschaft; Sprungform 1. Ordnung mit Kleinste-Quadrate-Gradient und Divergenz in Flussform; chirale Schichten und Kerne über den zentralen Multivektor K = k₊P₊ + k₋P₋); TwoPortLayerProblem (Schichten als Zweitore nach dem Vorbild der S-Matrix, E₂ auf der Außenfläche, E₁ auf dem Kern, g(s) über Partialbrüche mit dünnbesetzten Resolventen; Mehrfachschichten blockbidiagonal, jede Schicht auf ihrer eigenen Fläche, automatische Unterteilung dicker Schichten; mehrere Körper mit E₂ auf der Vereinigung und E₁ blockdiagonal) | mehrere rechte Seiten gleichzeitig, Wiederverwendung über Wellenlängen; Block-/HODLR-Vorkonditionierung für geschichtete Körper; neutrale Referenz mit gemeinsamem Aufbau; Dünnschicht: Krümmungssprünge an Rundungsübergängen, Nahfelder und Kräfte in chiralen Medien; Zweitor: Vorkonditionierer mit Kopplung benachbarter Schichten (Iterationen wachsen mit der Zahl der Teilschichten), Schichten nach innen |
| apps | Kugel, chirale Kugel, Würfel, Würfel-Dimer, Kugel-Dimer, Gmsh-Geometrien, Spektren (auch beschichtet), beschichtete Körper, Kompressionsbenchmark | Parameterstudien für AP 4 |
| prototype/resonance | 2D-Galerkin mit Streckung in log r, Absorber, angereicherte Eckelemente, 3D-Nullstellenanalyse | transparenter Kantenabschluss (diskrete DtN, Hardy-Raum-Ansatz) |
| bindings | – | pybind11-Modul für Skripting |

## Bekannte Grenzen

- **Resonanzfenster an scharfen 3D-Kanten:** Für Kontraste im kritischen Intervall (rechter Winkel: [−3, −1/3])
  laufen Eckwellen in die Kante, die am innersten Element reflektiert werden. In 2D löst Galerkin mit komplexer
  Streckung in log r das Problem; an 3D-Kanten erzeugt dieselbe Streckung Kernsingularitäten (z·z = 0 auf dem
  komplexen Rand), an Spitzen ist sie zulässig, ist im C++-Kern aber noch nicht umgesetzt
  (`results_resonance_window.md`). Praktischer Weg: abgerundete Kanten mit aufgelöstem Radius (`results_roundcube.md`).
- **Nahfeld auf gestreckten Elementen:** Die Zahl der Nahpaare ist groß (D < 2,5 h_max ist für die Fernfeldregel
  nötig); die Vorberechnung dominiert dort die Aufbauzeit. Sauter-Schwab wird nur bis Seitenverhältnis 1,6 verwendet.
- **Sauter-Schwab** setzt konforme Netze voraus (Nachbarschaft über gemeinsame Knoten).
- **Vorkonditionierung an Plasmonresonanzen:** lokale Blöcke und HODLR senken die Iterationen nur begrenzt; nahe
  Resonanzen ist T₁ schlecht konditioniert.
- **Dünne Schichten (d ≪ h):** Die Absolutwerte tragen einen systematischen Fehler von der Größe des
  Diskretisierungsfehlers (diskret E² ≠ 1); die Schichtwirkung ist als Differenz zu einer neutralen Rechnung auf
  denselben Netzen zu bestimmen (`results_coated.md`). Eindeutigkeit für verschachtelte Gebiete ist nicht bewiesen.
- **Dünnschicht-Näherung:** zweite Ordnung mit Restfehler ≈ 4 (d/a)² bei dielektrischen Schichten, grob um |ε_Metall| verstärkt,
  wenn die Kette durch Metall propagiert (Referenzfläche daher auf die Metallseite). Die Terme wachsen wie (d/h)², für d ≳ h
  steigen die Iterationen stark (d/h = 1,5: 476). Die Krümmung muss aufgelöst sein: am Silberwürfel mit 3 Elementen je
  Viertelrundung 5–13 % Unterschied zur exakten Rechnung (Diskretisierung, fällt mit dem Netz). Chirale Schichten und
  Kerne nur in der Dirac-Form; das Außenmedium muss achiral sein.
- **Zweitor:** Krümmung je (Teil-)Schicht bis O(d²); durch die automatische Unterteilung konvergent bis d/a = 0,5, dafür
  wachsen die Iterationen (d/a = 0,5: 156 bei n = 12). Schichten nach außen; mehrere Körper; Kugel, Kugel-Dimer und Gmsh.
- **CD dünner chiraler Schalen:** Die exakte Zwei-Flächen-Rechnung verfehlt ihn bei d/h ≲ 0,4 um 30–50 %; die
  Dünnschicht-Näherung ist dort vorzuziehen.
- **Rechenumgebung der Studien:** ein Kern, 3 GB Speicher; die Netzgrößen der Berichte sind dadurch begrenzt.
