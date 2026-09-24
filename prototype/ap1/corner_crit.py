"""Kritische Kontraste kappa = eps1/eps2 an einer Ecke: det sigma_Mellin(T_1)(tau; kappa) = 0.
Standardwahl a = 1 (mu1=mu2), b = sqrt(kappa): T_1 ist Laurent-Polynom in s = sqrt(kappa),
s*T_1(s) = s^2 A2 + s A1 + A0  -> quadratisches Eigenwertproblem."""
import numpy as np, scipy.linalg as sl
from corner import M_E, Jfull

def J_parts(alpha):
    # J(s) = J0 + s Jp + Jm / s   (Standardwahl, rm = 1, a = 1, b = s)
    s1, s2 = 2.0, 3.0
    J1, J2, J3 = [Jfull(alpha, s**2, 1.0, 1.0, s) for s in (1.0, s1, s2)]
    # J(s) linear in (1, s, 1/s): loese 3x3 System elementweise
    A = np.array([[1, 1, 1], [1, s1, 1/s1], [1, s2, 1/s2]])
    X = np.linalg.solve(A, np.stack([J1.ravel(), J2.ravel(), J3.ravel()]))
    return [x.reshape(16, 16) for x in X]

def crit_kappas(alpha, tau, M=None, c=0.5):
    M = M_E(alpha, tau, c) if M is None else M
    J0, Jp, Jm = J_parts(alpha); I = np.eye(16)
    A2 = 0.5*(I - M)@Jp; A1 = 0.5*(I + M) + 0.5*(I - M)@J0; A0 = 0.5*(I - M)@Jm
    # Linearisierung: [0 I; -A0 -A1] v = s [I 0; 0 A2] v
    Z = np.zeros((16, 16))
    L = np.block([[Z, I], [-A0, -A1]]); R = np.block([[I, Z], [Z, A2]])
    ev = sl.eigvals(L, R)
    ev = ev[np.isfinite(ev) & (np.abs(ev) > 1e-8) & (np.abs(ev) < 1e8)]
    return ev**2

if __name__ == '__main__':
    for alpha in [np.pi/2, 3*np.pi/2, np.pi]:
        print(f"alpha = {alpha/np.pi:.2f} pi")
        for tau in [0.0, 0.25, 0.5, 1.0, 2.0]:
            k = crit_kappas(alpha, tau)
            k = k[np.abs(k - 1) > 1e-6]          # kappa = 1 ist trivial (J = 1) -- filtern
            print(f"   tau={tau:4.2f}: ", np.round(np.sort_complex(k), 4))
