"""Gekruemmte Elemente aus Gmsh (v0.48): ein Gold-Nanostaebchen (Zylinder mit Halbkugelkappen) in Wasser, mit Gmsh in
zweiter Ordnung vernetzt (Mesh.ElementOrder 2, Kantenmitten auf der CAD-Flaeche), longitudinal angeregt. Fuer mehrere
Netzweiten wird die Extinktion mit ebenen Elementen (konstante Dichten, aus denselben Ecken) und mit gekruemmten Elementen
(lineare Dichten) gerechnet. Eine exakte Referenz gibt es fuer das Staebchen nicht; die gekruemmten Elemente konvergieren
bei viel groeberen Netzen.

Benoetigt das Python-Paket gmsh (pip install gmsh). Ohne gmsh rechnet das Skript die mitgelieferte Kugel
data/meshes/sphere_r1_order2_gmsh41.msh gegen Mie.

    PYTHONPATH=build/python python3 examples/python/gmsh_curved.py --h 0.6,0.45,0.35
"""
import argparse
import os
import sys
import tempfile
import time

import numpy as np

import cliffordbem as cb

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--h", default="0.6,0.45,0.35", help="Netzweiten (Einheiten von --unit)")
ap.add_argument("--unit", type=float, default=10.0, help="Laengeneinheit in nm (Staebchenradius)")
ap.add_argument("--length", type=float, default=3.0, help="Laenge des Zylinderteils (Einheiten)")
ap.add_argument("--lambda", dest="lam", type=float, default=700.0, help="Wellenlaenge in nm")
a = ap.parse_args()

water = cb.Medium(eps=1.33 ** 2)
gold = cb.make_material("Au").medium(a.lam)
om = cb.omega_from_wavelength(a.lam, a.unit)
so = cb.SolveOptions(tol=1e-8)


def solve(flat, quad, d, p):
    t0 = time.time()
    s_flat = cb.ScatteringProblem(flat, gold, om, outer=water).solve_plane_wave(d, p, so).sigma_ext
    t1 = time.time()
    s_curved = cb.CurvedScatteringProblem(quad, gold, om, outer=water).solve_plane_wave(d, p, so).sigma_ext
    return s_flat * a.unit ** 2, s_curved * a.unit ** 2, t1 - t0, time.time() - t1


try:
    import gmsh
except Exception:                                                       # ohne gmsh: mitgelieferte Kugel gegen Mie
    gmsh = None

if gmsh is None:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
    import mie
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "meshes", "sphere_r1_order2_gmsh41.msh")
    quad = cb.read_gmsh_quadratic(path)[0][1]
    flat = cb.read_gmsh(path)[0][1]
    ref = float(np.real(mie.qext(om, gold.eps, water.eps))) * np.pi * a.unit ** 2
    sf, sc, tf, tc = solve(flat, quad, (0, 0, 1), (1, 0, 0))
    print(f"gmsh nicht installiert: mitgelieferte Gmsh-Kugel ({len(quad)} Elemente), Radius {a.unit:g} nm, {a.lam:g} nm")
    print(f"  Mie {ref:.3f} nm^2;  eben/konstant {sf:.3f} ({sf / ref - 1:+.3%});  gekruemmt {sc:.3f} ({sc / ref - 1:+.4%})")
    print("  (fern der Resonanz ist die Extinktion fast reine Absorption ~ Im alpha; |Re eps|/Im eps = "
          f"{abs(gold.eps.real) / gold.eps.imag:.0f} verstaerkt Fehler der Polarisierbarkeit entsprechend)")
    sys.exit(0)


def rod_mesh(h, path):
    """Zylinder (Radius 1, Laenge L entlang z) mit Halbkugelkappen, zweite Ordnung."""
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("rod")
    L = a.length
    c = gmsh.model.occ.addCylinder(0, 0, -L / 2, 0, 0, L, 1.0)
    s1 = gmsh.model.occ.addSphere(0, 0, L / 2, 1.0)
    s2 = gmsh.model.occ.addSphere(0, 0, -L / 2, 1.0)
    gmsh.model.occ.fuse([(3, c)], [(3, s1), (3, s2)])
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(2, [s[1] for s in gmsh.model.getEntities(2)], 1)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h); gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.model.mesh.generate(2); gmsh.model.mesh.setOrder(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1); gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(path)
    gmsh.finalize()


L_nm = (a.length + 2) * a.unit
print(f"Gold-Staebchen {L_nm:g} x {2 * a.unit:g} nm in Wasser, {a.lam:g} nm, longitudinale Anregung (E parallel zur Achse)")
print(f"{'h':>6} {'Elemente':>9} | {'eben/konstant nm^2':>19} {'s':>6} | {'gekruemmt nm^2':>15} {'s':>6}")
rows = []
with tempfile.TemporaryDirectory() as tmp:
    for h in (float(v) for v in a.h.split(",")):
        f = os.path.join(tmp, f"rod_{h}.msh")
        rod_mesh(h, f)
        quad = cb.read_gmsh_quadratic(f)[0][1]
        flat = cb.read_gmsh(f)[0][1]
        sf, sc, tf, tc = solve(flat, quad, (1, 0, 0), (0, 0, 1))
        rows.append((len(quad), sf, sc))
        print(f"{h:6.3f} {len(quad):9d} | {sf:19.3f} {tf:6.1f} | {sc:15.4f} {tc:6.1f}", flush=True)
if len(rows) >= 2:
    best = rows[-1][2]
    print("Abweichung vom feinsten gekruemmten Wert:")
    for n, sf, sc in rows:
        print(f"  {n:5d} Elemente: eben/konstant {sf / best - 1:+.3%}, gekruemmt {sc / best - 1:+.4%}")
