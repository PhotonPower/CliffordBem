import sympy as sp
from cl3 import *

# 1. Algebraische Grundidentitaeten
I = PSEUDO
print("I^2 =", (I*I).c[0].real)
for e in (E1,E2,E3):
    assert ((I*e) - (e*I)).is_zero()
iI = 1j*I
print("(iI)^2 =", (iI*iI).c[0])
Pp = 0.5*(1+iI); Pm = 0.5*(1-iI)
assert (Pp*Pp - Pp).is_zero() and (Pm*Pm-Pm).is_zero() and (Pp*Pm).is_zero()
print("P+, P- orthogonale Projektoren: OK")
# Tangential-/Normalprojektion T(F) = (F - nFn)/2
import numpy as np
rng = np.random.default_rng(0)
n = rng.normal(size=3); n/=np.linalg.norm(n); N = MV.vec(n)
F = MV(rng.normal(size=8)+1j*rng.normal(size=8))
T = 0.5*(F - N*F*N); Nn = 0.5*(F + N*F*N)
assert ((N*T) + (T*N)).is_zero(); assert ((N*Nn) - (Nn*N)).is_zero()
print("T_n antikommutiert, N_n kommutiert mit n: OK")
