"""Stufe 1 der gekruemmten Elemente: Wie genau waere die BEM mit exakter Geometrie, aber weiterhin stueckweise konstanten
Dichten? Messung an der Kugel gegen Mie, ohne den Kern zu aendern.

Konstruktion. Die Ikosaederkugel n entsteht durch radiale Projektion eines baryzentrischen Gitters auf den Ikosaederflaechen;
das Gitter n m enthaelt alle Punkte des Gitters n, und grobe Kanten liegen auf denselben Grosskreisen. Das gekruemmte grobe
Element J (Kugeldreieck) wird daher exakt von m^2 feinen Dreiecken i parkettiert. Stueckweise konstante Dichten auf den
gekruemmten Elementen bilden den Teilraum mit der Basis

    P_iJ = sqrt(a_i / A_J)   (i in J),   A_J = sum_{i in J} a_i,   Spalten orthonormal (P^T P = 1).

Die Galerkin-Matrix auf diesem Teilraum ist T_grob = P^T T_fein P -- einschliesslich der im Element variierenden Normalen,
die im Cauchy-Operator und in der Transmissionsabbildung J steht (fein: je Dreieck konstant). Rechte Seite P^T b_fein,
Streuspur P (h - b). Fuer m -> oo geht das in gekruemmte Elemente mit exakter Geometrie ueber; der Restfehler der feinen
Polyeder faellt wie 1/m^2 und wird extrapoliert. m = 1 ist die heutige Rechnung mit ebenen Dreiecken.

Erweiterung (Frage nach Stufe 2: lohnen gekruemmte Elemente erst mit Dichten hoeherer Ordnung?):
  --geometry flat   die feinen Punkte werden zentral auf die Ebene ihres groben Dreiecks projiziert (Grosskreise -> Geraden,
                    die Unterteilung bleibt konform): ebene Elemente, gleich fein integriert
  --space p1        unstetig lineare Dichten je grobem Element (3 Funktionen je Blade, baryzentrische Koordinaten der
                    Parametrisierung), auf dem feinen Netz als Mittelwerte je feinem Dreieck dargestellt (Mittelwert einer
                    linearen Funktion = Wert im Schwerpunkt), je Element orthonormiert. Die Darstellung naehert die echten
                    linearen Funktionen mit O(h/m); auch dieser Anteil wird in m extrapoliert.

    PYTHONPATH=build/python python3 prototype/curved/curved_galerkin.py --material glass --out results.json
"""
import argparse
import json
import os
import sys
import time

import numpy as np

import cliffordbem as cb

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import mie  # noqa: E402

CASES = {  # wie docs/results_scattering.md
    "glass": dict(eps=2.25, omega=1.0),
    "gold": dict(eps=-11 + 1.2j, omega=0.5),
}


def coarse_map(fine, coarse):
    """Grobes Element je feinem Dreieck (Schwerpunkt im Kegel ueber dem groben Kugeldreieck) und baryzentrische
    Koordinaten des Schwerpunkts in der Parametrisierung des groben Elements."""
    cf = fine.centroids / np.linalg.norm(fine.centroids, axis=1, keepdims=True)
    V = coarse.points[coarse.triangles]                                   # (Nc, 3, 3) Ecken auf der Kugel
    cc = coarse.centroids / np.linalg.norm(coarse.centroids, axis=1, keepdims=True)
    cand = np.argsort(-(cf @ cc.T), axis=1)[:, :8]                        # naechste grobe Elemente
    owner = np.full(len(cf), -1)
    bary = np.zeros((len(cf), 3))
    for r in range(cand.shape[1]):
        todo = owner < 0
        if not todo.any():
            break
        J = cand[todo, r]
        A = np.transpose(V[J], (0, 2, 1))                                 # Spalten = Ecken
        lam = np.linalg.solve(A, cf[todo][:, :, None])[:, :, 0]           # u = l1 A + l2 B + l3 C
        ok = np.all(lam > -1e-9, axis=1)
        idx = np.where(todo)[0][ok]
        owner[idx] = J[ok]
        bary[idx] = lam[ok] / lam[ok].sum(axis=1, keepdims=True)
    if (owner < 0).any():
        raise RuntimeError("Zuordnung fein -> grob unvollstaendig")
    counts = np.bincount(owner, minlength=len(coarse))
    return owner, counts, bary


def flatten(fine, coarse, owner):
    """Feine Punkte zentral auf die Ebene ihres groben Dreiecks projizieren (ebene Elemente, konforme Unterteilung)."""
    P = fine.points.copy()
    T = fine.triangles
    n = coarse.normals[owner]                                             # (Nf, 3)
    A = coarse.points[coarse.triangles[owner, 0]]
    for c in range(3):
        p = fine.points[T[:, c]]
        P[T[:, c]] = p * (np.sum(n * A, axis=1) / np.sum(n * p, axis=1))[:, None]
    return cb.TriangleMesh(P, T)


class Projection:
    """P (8 Nf x 8 Nc s) mit s = 1 (konstant) oder 3 (linear) Funktionen je Element und Blade; Spalten orthonormal
    bezueglich des L2-Produkts der normierten feinen Basis (Koeffizient i = sqrt(a_i) * Mittelwert auf Dreieck i)."""

    def __init__(self, fine, owner, Nc, bary=None, space="pc"):
        a = fine.areas
        self.owner, self.Nc, self.Nf = owner, Nc, len(a)
        self.A = np.bincount(owner, weights=a, minlength=Nc)             # Elementflaechen (fein summiert)
        self.s = 1 if space == "pc" else 3
        self.Q = np.zeros((self.Nf, self.s))
        for J, idx in enumerate(np.split(np.argsort(owner, kind="stable"), np.cumsum(np.bincount(owner, minlength=Nc))[:-1])):
            B = np.sqrt(a[idx])[:, None] * (np.ones((len(idx), 1)) if self.s == 1 else bary[idx])
            q, _ = np.linalg.qr(B)
            if q.shape[1] < self.s:
                raise RuntimeError("Element zu grob unterteilt fuer lineare Dichten (m >= 2 noetig)")
            self.Q[idx] = q

    @property
    def size(self):
        return 8 * self.Nc * self.s

    def up(self, xc):        # P xc
        x = xc.reshape(self.Nc, self.s, 8)[self.owner]                    # (Nf, s, 8)
        return np.einsum("is,isb->ib", self.Q, x).ravel()

    def down(self, xf):      # P^T xf
        out = np.zeros((self.Nc, self.s, 8), complex)
        np.add.at(out, self.owner, self.Q[:, :, None] * xf.reshape(self.Nf, 1, 8))
        return out.ravel()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--material", choices=CASES, default="glass")
    ap.add_argument("--coarse", default="4,5,6,8", help="grobe Ikosaederkugeln n")
    ap.add_argument("--max-fine", type=int, default=20, help="feinstes Gitter n m (20: 8000 Dreiecke)")
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--geometry", choices=["curved", "flat"], default="curved")
    ap.add_argument("--space", choices=["pc", "p1"], default="pc")
    ap.add_argument("--min-m", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    case = CASES[a.material]
    eps, om = case["eps"], case["omega"]
    q_mie = mie.qext(om, eps)
    coarse_n = [int(v) for v in a.coarse.split(",")]
    fine_n = sorted({n * m for n in coarse_n for m in range(a.min_m, a.max_fine // n + 1)})
    coarse_mesh = {n: cb.make_icosphere(n) for n in coarse_n}
    rows = json.load(open(a.out)) if os.path.exists(a.out) else []
    done = {(r["n"], r["m"]) for r in rows if r.get("geometry", "curved") == a.geometry and r.get("space", "pc") == a.space}
    print(f"{a.material}: eps = {eps}, omega a = {om}, Mie Q_ext = {q_mie:.6f}; feine Netze {fine_n}", flush=True)
    for nf in fine_n:
        todo = [n for n in coarse_n if nf % n == 0 and (n, nf // n) not in done and nf // n >= a.min_m]
        if a.geometry == "flat":                                          # ein feines Netz je grobem Netz
            groups = [(n, [n]) for n in todo]
        else:
            groups = [(None, todo)] if todo else []
        sphere = cb.make_icosphere(nf)
        for flat_of, members in groups:
            maps = {n: coarse_map(sphere, coarse_mesh[n]) for n in members}
            fine = sphere if flat_of is None else flatten(sphere, coarse_mesh[flat_of], maps[flat_of][0])
            run_fine(fine, nf, members, maps, coarse_mesh, a, eps, om, q_mie, rows)   # gibt das feine Problem wieder frei


def run_fine(fine, nf, members, maps, coarse_mesh, a, eps, om, q_mie, rows):
    """Ein feines Problem aufbauen und fuer alle groben Netze, die es verfeinert, das projizierte System loesen."""
    t0 = time.time()
    Pf = cb.ScatteringProblem(fine, cb.Medium(eps=eps), om)
    t_build = time.time() - t0
    bf = cb.plane_wave_trace(fine, cb.Medium(), om, (0, 0, 1), (1, 0, 0))
    T = Pf.T
    for n in members:
        m = nf // n
        t1 = time.time()
        owner, counts, bary = maps[n]
        assert np.all(counts == m * m), "jedes grobe Element muss genau m^2 feine Dreiecke enthalten"
        Pr = Projection(fine, owner, len(coarse_mesh[n]), bary, a.space)
        bc = Pr.down(bf)
        x, info = cb.gmres(lambda v: Pr.down(T.apply(Pr.up(v))), bc, M=lambda v: Pr.down(T.precondition(Pr.up(v))),
                           tol=a.tol, restart=200, max_iter=2000)
        hs = Pr.up(x - bc)
        sig = cb.extinction_cross_section(fine, hs, om, 1.0, (0, 0, 1), (1, 0, 0))   # Aussenraum Vakuum: k = omega
        q = sig / np.pi
        row = dict(n=n, m=m, geometry=a.geometry, space=a.space, N=len(coarse_mesh[n]), N_fine=len(fine), Q=q,
                   err=q / q_mie - 1, area=float(Pr.A.sum()), its=info.iterations, conv=bool(info.converged),
                   t_build=t_build, t_solve=time.time() - t1)
        if m == 1 and a.space == "pc":                                    # Gegenprobe: gleich ScatteringProblem
            row["Q_direct"] = Pf.solve_plane_wave((0, 0, 1), (1, 0, 0), cb.SolveOptions(tol=a.tol)).sigma_ext / np.pi
        rows.append(row)
        json.dump(rows, open(a.out, "w"), indent=1)
        print(f"  {a.geometry}/{a.space} n = {n:2d} ({row['N']:5d} Elemente), m = {m} ({row['N_fine']:5d} fein): "
              f"Q = {q:.6f}, Fehler {row['err']:+.4%}, {info.iterations} It.  [Aufbau {t_build:.0f} s]", flush=True)

if __name__ == "__main__":
    main()
