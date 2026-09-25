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
