"""Mie-Referenz (Bohren/Huffman), nichtmagnetisch, Kugelradius 1, aeusseres Medium eps2."""
import numpy as np
from scipy.special import spherical_jn, spherical_yn
def psi(n, z):  return z*spherical_jn(n, z)
def dpsi(n, z): return spherical_jn(n, z) + z*spherical_jn(n, z, derivative=True)
def xi(n, z):   return z*(spherical_jn(n, z) + 1j*spherical_yn(n, z))
def dxi(n, z):  return (spherical_jn(n, z)+1j*spherical_yn(n, z)) + z*(spherical_jn(n, z, True)+1j*spherical_yn(n, z, True))
def coeffs(om, e1, e2=1.0, nmax=30):
    x = om*np.sqrt(e2); m = np.sqrt(complex(e1/e2)); n = np.arange(1, nmax+1)
    Da = m*psi(n, m*x)*dxi(n, x) - xi(n, x)*dpsi(n, m*x)
    Db = psi(n, m*x)*dxi(n, x) - m*xi(n, x)*dpsi(n, m*x)
    a = (m*psi(n, m*x)*dpsi(n, x) - psi(n, x)*dpsi(n, m*x))/Da
    b = (psi(n, m*x)*dpsi(n, x) - m*psi(n, x)*dpsi(n, m*x))/Db
    return n, a, b, Da, Db
def qext(om, e1, e2=1.0):
    n, a, b, _, _ = coeffs(om, e1, e2); x = om*np.sqrt(e2)
    return 2/x**2*np.sum((2*n+1)*np.real(a+b))
