"""Beispielgeometrien mit Gmsh (Python-API, OpenCASCADE) als ASCII-.msh 4.1 fuer scatter_mesh.
Jeder Koerper erhaelt eine physikalische Flaechengruppe (Tag 1, 2, ...) aus allen seinen Randflaechen.
  python3 tools/make_geometries.py sphere 0.15 examples/sphere_h015.msh
  python3 tools/make_geometries.py rod 0.12 examples/rod.msh            (Stab: Laenge 4, Radius 0.5, Achse x)
  python3 tools/make_geometries.py bornkuhn 0.1 examples/bornkuhn.msh --angle 60   (zwei gekreuzte Staebe)"""
import sys, math, argparse
import gmsh
def capsule(L, r, angle_deg=0.0, z=0.0):
    occ = gmsh.model.occ
    c = occ.addCylinder(-L/2 + r, 0, 0, L - 2*r, 0, 0, r)
    s1 = occ.addSphere(-L/2 + r, 0, 0, r); s2 = occ.addSphere(L/2 - r, 0, 0, r)
    v, _ = occ.fuse([(3, c)], [(3, s1), (3, s2)])
    if angle_deg: occ.rotate(v, 0, 0, 0, 0, 0, 1, math.radians(angle_deg))
    if z: occ.translate(v, 0, 0, z)
    return v
def main():
    ap = argparse.ArgumentParser(); ap.add_argument('kind'); ap.add_argument('h', type=float); ap.add_argument('out')
    ap.add_argument('--angle', type=float, default=60.0); ap.add_argument('--gap', type=float, default=1.2)
    ap.add_argument('--length', type=float, default=3.0); ap.add_argument('--radius', type=float, default=0.3)
    a = ap.parse_args()
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add(a.kind)
    occ = gmsh.model.occ
    if a.kind == 'sphere': vols = [[(3, occ.addSphere(0, 0, 0, 1.0))]]
    elif a.kind == 'rod': vols = [capsule(4.0, 0.5)]
    elif a.kind == 'bornkuhn': vols = [capsule(a.length, a.radius, 0.0, -a.gap/2), capsule(a.length, a.radius, a.angle, +a.gap/2)]
    else: raise SystemExit('unbekannte Geometrie')
    occ.synchronize()
    for i, v in enumerate(vols):
        surf = [s[1] for s in gmsh.model.getBoundary(v, oriented=False)]
        gmsh.model.addPhysicalGroup(2, surf, i + 1)
    gmsh.option.setNumber("Mesh.MeshSizeMax", a.h); gmsh.option.setNumber("Mesh.MeshSizeMin", 0.3*a.h)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1); gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(a.out); gmsh.finalize()
if __name__ == '__main__': main()
