"""Geschlossene Form des kritischen Satzes und gewichtete Raeume.
kappa kritisch auf der Mellin-Geraden Re a = c  <=>  q(kappa) = (1+kappa)/(1-kappa) = +-S(c+i tau),
S(a) = sin((pi-alpha) a)/sin(pi a).  Innen-Test: Windung von S(a)^2 - q^2 ueber tau."""
import numpy as np
def S(alpha, a): return np.sin((np.pi - alpha)*a)/np.sin(np.pi*a)
taus = np.linspace(-40, 40, 40001)
def wind_scalar(kappa, alpha=np.pi/2, c=0.5):
    q = (1 + kappa)/(1 - kappa); g = S(alpha, c + 1j*taus)**2 - q**2
    ph = np.unwrap(np.angle(g)); return (ph[-1] - ph[0])/(2*np.pi)
def c_star(kappa, alpha=np.pi/2):
    """groesstes c in (0, 1/2], fuer das kappa ausserhalb liegt (Windung 0)."""
    cs = np.linspace(0.5, 0.005, 400)
    for c in cs:
        if abs(wind_scalar(kappa, alpha, c)) < 0.5: return c
    return np.nan
if __name__ == '__main__':
    from corner_index import winding
    print("Abgleich Windung (skalar, geschl. Form) vs. Dirac-Symbol (Standardwahl), c = 1/2:")
    for k in [-3+0.3j, -2+0.3j, -0.5+0.1j, -5+0.5j, -6.215+0.678j, -11+1.2j, -2+3.5j, 2.25]:
        print(f"   kappa={k}: skalar {wind_scalar(k):+.2f}   Dirac {winding(k)[0]:+.2f}")
    print("\nNoetiges Gewicht r^(2 beta) an einer 90-Grad-Kante, beta = 1/2 - c*:")
    for k in [-5+0.5j, -3+0.3j, -2+0.3j, -1.5+0.2j, -1.2+0.1j, -0.5+0.1j, -6.215+0.678j]:
        cs = c_star(k); print(f"   kappa={k}: c* = {cs:.3f}  ->  beta >= {0.5 - cs:.3f}")
    # Energiewahl: ist die physikalische Kurve auch dort singulaer?
    from corner import M_E, Jfull
    kc = -2.009 - 1.3615j   # auf der Kurve bei tau = 0.4
    for t in [0.4, -0.4]:
        T = 0.5*(np.eye(16) + M_E(np.pi/2, t)) + 0.5*(np.eye(16) - M_E(np.pi/2, t))@Jfull(np.pi/2, kc, 1.0, np.conj(kc), np.sqrt(kc))
        Ts = 0.5*(np.eye(16) + M_E(np.pi/2, t)) + 0.5*(np.eye(16) - M_E(np.pi/2, t))@Jfull(np.pi/2, kc)
        sv = np.linalg.svd(T, compute_uv=False); svs = np.linalg.svd(Ts, compute_uv=False)
        print(f"kappa auf der Kurve, tau={t}: kleinster SW Energiewahl {sv[-1]:.1e}, Standardwahl {svs[-1]:.1e}")
