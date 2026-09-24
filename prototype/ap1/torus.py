"""Torus-Geometrie fuer die Subraumwinkel-Tests (Schritt 1.4: nicht einfach zusammenhaengendes Gebiet)."""
import numpy as np
from mfs import Sphere

def torus_points(R, r, nu, nv, shift=0.0):
    u = 2*np.pi*(np.arange(nu)+shift)/nu; v = 2*np.pi*(np.arange(nv)+shift)/nv
    U, V = np.meshgrid(u, v, indexing='ij'); U, V = U.ravel(), V.ravel()
    X = np.stack([(R+r*np.cos(V))*np.cos(U), (R+r*np.cos(V))*np.sin(U), r*np.sin(V)], -1)
    N = np.stack([np.cos(V)*np.cos(U), np.cos(V)*np.sin(U), np.sin(V)], -1)
    W = r*(R+r*np.cos(V))*(2*np.pi/nu)*(2*np.pi/nv)
    return X, N, W

class Torus(Sphere):
    def __init__(self, nu=64, nv=32, su=24, sv=10, R=1.0, r=0.4, fin=0.5, fout=1.75):
        self.X, self.N, self.W = torus_points(R, r, nu, nv)
        self.src_in = torus_points(R, fin*r, su, sv, 0.5)[0]
        self.src_out = torus_points(R, fout*r, su, sv, 0.5)[0]

if __name__ == '__main__':
    from mfs import min_angle
    X, N, W = torus_points(1, 0.4, 64, 32)
    print("Oberflaeche num/exakt:", W.sum(), 4*np.pi**2*0.4)
    for args in [dict(nu=48,nv=24,su=18,sv=8), dict(nu=64,nv=32,su=24,sv=10), dict(nu=80,nv=40,su=30,sv=12)]:
        t = Torus(**args)
        print(args, ["%.5f" % min_angle(t, kd, om, e1, 1, 1.0, 1) for kd in ['maxwell','dirac'] for om, e1 in [(1.0,2.25),(1e-3,2.25)]], flush=True)
