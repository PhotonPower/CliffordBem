"""Pruefung der Flussidentitaet div J = -2 Im(k) |G|^2,  J_j = <G^dagger e_j G>_0,
fuer eine allgemeine Multivektor-Dirac-Loesung (alle Grade) mit komplexem k."""
import numpy as np
from fields import gp, dirac_kernel, blade
def dagger(G):   # Reversion + komplexe Konjugation
    s = np.array([1,1,1,-1,1,-1,-1,-1])   # Grade: 0,1,1,2,1,2,2,3 -> (+,+,+,-,+,-,-,-)
    return np.conj(G)*s
rng = np.random.default_rng(3)
c = rng.normal(size=8)+1j*rng.normal(size=8)
for k in [1.7, 1.2+0.6j, 0.8j]:
    G = lambda x: gp(dirac_kernel(x, np.array([0.1,-0.2,0.05]), k), c)
    Jf = lambda x: np.stack([gp(gp(dagger(G(x)), blade(b)), G(x))[...,0] for b in [1,2,4]], -1)
    x = np.array([[0.9,0.4,-0.6]]); h = 1e-5
    div = sum((Jf(x+h*np.eye(3)[j])[0,j]-Jf(x-h*np.eye(3)[j])[0,j])/(2*h) for j in range(3))
    rhs = -2*np.imag(k)*np.sum(np.abs(G(x))**2)
    print(f"k={k}: div J = {div.real:+.8f}{div.imag:+.1e}i,  -2 Im k |G|^2 = {rhs:+.8f}")
