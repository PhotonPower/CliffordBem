"""Eigene Formulierungen aus den Bausteinen des Kerns: Cauchy-Operatoren als H-Matrizen aus C++, die Algebra und
die Gleichung selbst in NumPy. Hier wird T_1 = 1/2 (1 + E_2) + 1/2 (1 - E_1) J von Hand zusammengesetzt, mit dem
Operator des Kerns verglichen und mit GMRES und einer in Python geschriebenen Vorkonditionierung geloest.

    PYTHONPATH=build/python python3 examples/python/custom_formulation.py --n 6
"""
import argparse

import numpy as np

import cliffordbem as cb

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--n", type=int, default=6)
ap.add_argument("--omega", type=float, default=0.5)
a = ap.parse_args()

inner, outer = cb.Medium(eps=-11 + 1.2j), cb.Medium(eps=1.7689)
m = cb.make_icosphere(a.n)
N = len(m)

# Cauchy-Operatoren E_k (Eintraege -> H-Matrix -> Operator); keep_alive haelt die Kette am Leben
E1 = cb.CauchyOperator(m, cb.KernelHMatrix(cb.KernelEntries(m, inner.k(a.omega)), cb.HMatrixParams(eps=1e-6)))
E2 = cb.CauchyOperator(m, cb.KernelHMatrix(cb.KernelEntries(m, outer.k(a.omega)), cb.HMatrixParams(eps=1e-6)))

# 1) Plemelj: E h = +h fuer Spuren innerer Loesungen (ebene Welle), E h = -h fuer Spuren aeusserer (Streuspur; unten)
d, p = (0, 0, 1), (1, 0, 0)
b = cb.plane_wave_trace(m, outer, a.omega, d, p)
print(f"ebene Welle (innere Loesung): |E2 b - b| / |b| = {np.linalg.norm(E2 @ b - b) / np.linalg.norm(b):.2e}")

# 2) J je Dreieck aus der Algebra (8x8 je Normale), T_1 von Hand
J = np.array([cb.transmission_map(nrm, inner, outer) for nrm in m.normals])        # (N, 8, 8)


def apply_J(x, Js=J):
    return np.einsum("tij,tj->ti", Js, x.reshape(N, 8)).ravel()


def T1(x):
    Jx = apply_J(x)
    return 0.5 * (x + E2 @ x) + 0.5 * (Jx - E1 @ Jx)


T = cb.TransmissionOperator(m, E1, E2, inner, outer)
h = np.random.default_rng(0).normal(size=8 * N) + 0j
print(f"|T1_python - T1_kern| / |T1| = {np.linalg.norm(T1(h) - T.apply(h)) / np.linalg.norm(T.apply(h)):.1e}")

# 3) punktweise Vorkonditionierung 2 (1 + J)^-1, in NumPy
Pinv = 2 * np.linalg.inv(np.eye(8)[None] + J)
x, info = cb.gmres(T1, b, M=lambda v: apply_J(v, Pinv), tol=1e-8)
print(f"GMRES (Operator und Vorkonditionierer in Python): {info}")
k = outer.k(a.omega)
s_py = cb.extinction_cross_section(m, x - b, k, outer.eps, d, p)
s_ref = cb.ScatteringProblem(m, inner, a.omega, outer=outer, hmatrix=cb.HMatrixParams(eps=1e-6)).solve_plane_wave(d, p, cb.SolveOptions(tol=1e-8)).sigma_ext
print(f"sigma_ext: Python {s_py:.8f}, ScatteringProblem {s_ref:.8f}")
hs = x - b
print(f"Streuspur (aeussere Loesung): |E2 hs + hs| / |hs| = {np.linalg.norm(E2 @ hs + hs) / np.linalg.norm(hs):.2e}")

# 4) Fernfeld als Multivektor: F_inf = sqrt(eps) E_inf + I sqrt(mu) H_inf, H_inf = d x E_inf in Vorwaertsrichtung
F = cb.far_field(m, x - b, k, d)
E_inf = F.vector_part() / np.sqrt(outer.eps)
print(f"F_inf(z) = {F}\nE_inf(z) = {np.round(E_inf, 6)}")
