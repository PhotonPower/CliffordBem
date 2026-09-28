# Beispiele

Die Netze werden mit `tools/make_geometries.py` erzeugt (Gmsh-Python-API, `pip install gmsh`) und nicht
eingecheckt:

```bash
python3 tools/make_geometries.py sphere 0.13 examples/sphere.msh
python3 tools/make_geometries.py rod 0.12 examples/rod.msh
python3 tools/make_geometries.py bornkuhn 0.12 examples/bornkuhn_60.msh --angle 60
python3 tools/make_geometries.py bornkuhn 0.12 examples/bornkuhn_-60.msh --angle -60
./build/scatter_mesh --mesh examples/bornkuhn_60.msh --omega 0.5 --media "-11,1.2" --pol circ
```

Beschichtete Körper (Oxidschicht 2 nm auf einem abgerundeten Silberwürfel, Kante 50 nm, in Wasser; die Schicht
wirkt als Differenz zur neutralen Rechnung mit Wasser als Schichtmaterial, siehe `docs/results_coated.md`):

```bash
python3 tools/make_geometries.py roundcube 0.3 examples/rc04.msh --radius 0.4 --curv 12
./build/spectrum --mesh examples/rc04.msh --unit 25 --materials Ag --nbg 1.33 --lambda 400:480:20 --coating "2:2.89,0"
./build/spectrum --mesh examples/rc04.msh --unit 25 --materials Ag --nbg 1.33 --lambda 400:480:20 --coating "2:1.7689,0"
./build/scatter_coated --n 4,8 --omega 0.5 --core -11,1.2 --coat 0.02,2.25,0 --neutral   # Kugel gegen tools/mie_coated.py
./build/scatter_coated --n 4,8 --omega 0.5 --core -11,1.2 --coat 0.05,2.25,0 --thin --bare   # Dünnschicht-Näherung 2. Ordnung
./build/spectrum --mesh examples/rc04.msh --unit 25 --materials Ag --nbg 1.33 --lambda 400:480:20 --coating "2:2.89,0" --thin 0
```
