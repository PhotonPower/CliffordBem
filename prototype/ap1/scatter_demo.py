"""Streuung einer ebenen Welle (d = e3, p = e1) an der Einheitskugel, geloest
(a) im VOLLEN Multivektorraum (Dirac-Kerne mit allen 8 Blades, Hilfsgrade 0/3 frei) und
(b) im Maxwell-Raum (elektrische + magnetische Dipole).
Prueft: Hilfskomponenten? Uebereinstimmung (a) vs (b)? Extinktion vs. Mie?"""
import numpy as np
from mfs import *
from fields import gp, vec, PSEUDO, blade
import mie

def plane_wave(X, k, eps, d=np.array([0,0,1.]), p=np.array([1,0,0.])):
    ph = np.exp(1j*k*X@d)[:,None]
    return vec(np.sqrt(eps)*p[None,:]*ph) + gp(PSEUDO, vec(np.sqrt(eps)*np.cross(d,p)[None,:]*ph))

def basis(x, kind, k, eps, S):
    cols = []
    for xs in S:
        if kind == 'dirac':
            K = dirac_kernel(x, xs, k); cols += [gp(K, blade(b)) for b in range(8)]
        else:
            for p in np.eye(3):
                F = maxwell_dipole(x, xs, p, k, eps); cols += [F, gp(PSEUDO, F)]
    return np.stack(cols, -1)

def solve(s, kind, om, e1, e2=1.0, b=None, a=1.0):
    k1, k2 = om*np.sqrt(e1), om*np.sqrt(e2); re_ = e1/e2
    b = np.sqrt(re_) if b is None else b
    if a == 'energy': a, b = np.conj(re_), np.sqrt(re_)      # mu1 = mu2
    A = basis(s.X, kind, k1, e1, s.src_out); B = basis(s.X, kind, k2, e2, s.src_in)
    hin = plane_wave(s.X, k2, e2)
    JB = s.J(B, re_, 1.0, a, b); Jh = s.J(hin[:,:,None], re_, 1.0, a, b)[:,:,0]
    w = np.sqrt(s.W)[:,None]
    M = np.concatenate([(w[:,:,None]*A).reshape(-1,A.shape[-1]), -(w[:,:,None]*JB).reshape(-1,B.shape[-1])], 1)
    rhs = (w*Jh).reshape(-1)
    U, sv, Vh = np.linalg.svd(M, full_matrices=False); keep = sv > 1e-12*sv[0]
    c = Vh[keep].conj().T @ ((U[:,keep].conj().T @ rhs)/sv[keep])
    nA = A.shape[-1]
    h1, h2 = A@c[:nA], hin + B@c[nA:]
    L2 = lambda H: np.sqrt(np.sum(s.W[:,None]*np.abs(H)**2))
    out = dict(res=np.linalg.norm(M@c-rhs)/np.linalg.norm(rhs),
               aux_ext=L2(h2[:,[0,7]])/L2(h2), aux_int=L2(h1[:,[0,7]])/L2(h1))
    out['F_int'] = lambda x: basis(x, kind, k1, e1, s.src_out)@c[:nA]
    out['F_sca'] = lambda x: basis(x, kind, k2, e2, s.src_in)@c[nA:]
    # optisches Theorem: sigma_ext = 4pi/k Im(p . E_inf(d)),  E_s ~ e^{ikr}/r E_inf
    def q(Rf):
        Fs = out['F_sca'](np.array([[0,0,Rf]]))[0]
        Einf = np.array([Fs[1], Fs[2], Fs[4]])/np.sqrt(e2)*Rf*np.exp(-1j*k2*Rf)
        return (4*np.pi/k2*np.imag(Einf[0]))/np.pi
    q1, q2 = q(2e4), q(2e5)
    out['Qext'] = q2 - (q1 - q2)/9          # Richardson-Extrapolation in 1/R
    return out

if __name__ == '__main__':
    import sys
    AUX = 'energy' if 'energy' in sys.argv else 1.0
    print('Hilfsparameter:', 'Energiewahl' if AUX == 'energy' else 'Standardwahl')
    s = Sphere(20, 40, 70)
    print(f"{'Fall':22s} {'Res D':>8s} {'Res M':>8s} {'Hilfs a.':>9s} {'Hilfs i.':>9s} {'|D-M| i.':>9s} {'|D-M| a.':>9s} {'Qext D':>9s} {'Qext M':>9s} {'Qext Mie':>9s}")
    for om, e1, name in [(1.0, 2.25, 'Glas, ka=1'), (2.0, 2.25, 'Glas, ka=2'), (np.pi/1.5, 2.25, 'Glas, k1a=pi'),
                         (0.5, -11+1.2j, 'Gold, ka=0.5'), (0.3, -2.2+0.3j, 'nahe Dipolplasmon')]:
        D, M = solve(s, 'dirac', om, e1, a=AUX), solve(s, 'maxwell', om, e1)
        xi_, xo = np.array([[0.2,-0.3,0.4]]), np.array([[1.5,0.9,-1.1]])
        di = np.linalg.norm(D['F_int'](xi_)-M['F_int'](xi_))/np.linalg.norm(M['F_int'](xi_))
        do = np.linalg.norm(D['F_sca'](xo)-M['F_sca'](xo))/np.linalg.norm(M['F_sca'](xo))
        print(f"{name:22s} {D['res']:8.1e} {M['res']:8.1e} {D['aux_ext']:9.1e} {D['aux_int']:9.1e} {di:9.1e} {do:9.1e} "
              f"{D['Qext']:9.5f} {M['Qext']:9.5f} {mie.qext(om, e1):9.5f}", flush=True)
