# Python-Prototypen

Referenzimplementierungen aus der Ausarbeitung der Arbeitspakete (NumPy/SciPy/SymPy/mpmath).
- `ap1/`: Theorie-Prüfskripte (Clifford-Algebra `cl3.py`, Symbolanalyse, Kugel/Torus-Studien,
  Kanten- und Kegelanalyse, chirale Medien, 2D-Quadrat-Prototyp). Übersicht der Skripte und
  Ergebnisse im Arbeitspapier `docs/papers/ap1/AP1_Ausarbeitung.pdf`.
- `ap2/`: 3D-Galerkin-Prototyp (Kugel, Würfel, Symmetriesektoren, Vorkonditionierung) und
  H-Matrix-Prototyp (`hmat.py`, AP 3).

- `curved/`: gekrümmte Elemente, Stufe 1 (v0.44) und 1b (v0.46, `--geometry quadratic`): Messung über Galerkin-Projektion auf feinen Unterteilungen mit der
  Python-Anbindung (`PYTHONPATH=build/python`, aus dem Wurzelverzeichnis starten); Ergebnisse in `docs/results_curved.md`.

Die Skripte in `ap1/` und `ap2/` erwarten, aus ihrem eigenen Verzeichnis gestartet zu werden.
