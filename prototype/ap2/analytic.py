"""Analytische Dreiecksintegrale (vektorisiert fuer viele Punkte x, ein Dreieck p0,p1,p2 mit Normale n):
  I_grad(x) = int_T (x - y)/|x - y|^3 dS_y   (= int grad_y 1/R),   I_inv(x) = int_T 1/|x - y| dS_y."""
import numpy as np
def tri_integrals(X, p, n):
    """X: (M,3); p: (3,3) Ecken (gegen den Uhrzeigersinn um n); n: Einheitsnormale. Rueckgabe (M,3), (M,)."""
    w = (X - p[0])@n                                   # Hoehe
    rho = X - w[:, None]*n[None, :]                    # Projektion in die Ebene
    Igrad_par = np.zeros_like(X); Iinv = np.zeros(len(X)); beta_sum = np.zeros(len(X))
    for e in range(3):
        q1, q2 = p[e], p[(e+1) % 3]
        L = np.linalg.norm(q2 - q1); t = (q2 - q1)/L; m = np.cross(t, n)      # aeussere Kantennormale in der Ebene
        lm = (q1[None, :] - rho)@t; lp = (q2[None, :] - rho)@t              # l^-, l^+
        P0 = (q1[None, :] - rho)@m                                         # vorzeichenbehafteter Kantenabstand
        Rm = np.linalg.norm(X - q1[None, :], axis=1); Rp = np.linalg.norm(X - q2[None, :], axis=1)
        R0sq = P0**2 + w**2
        # f_e = int_e dl/R = ln((R+ + l+)/(R- + l-)), stabil fuer l < 0 ueber ln((R- - l-)/(R+ - l+))
        with np.errstate(divide='ignore', invalid='ignore'):
            f1 = np.log((Rp + lp)/(Rm + lm)); f2 = np.log((Rm - lm)/(Rp - lp))
        f = np.where(lm + lp >= 0, f1, f2)
        Igrad_par += f[:, None]*m[None, :]
        aw = np.abs(w)
        with np.errstate(divide='ignore', invalid='ignore'):
            beta = np.arctan2(P0*lp, R0sq + aw*Rp) - np.arctan2(P0*lm, R0sq + aw*Rm)
        beta_sum += beta
        Iinv += P0*f
    Iinv -= np.abs(w)*beta_sum
    # Normalkomponente: int w/R^3 dS = sign(w) * Raumwinkel = |w|-gewichtet ueber beta
    Inorm = np.sign(w)*beta_sum
    return Igrad_par + Inorm[:, None]*n[None, :], Iinv
if __name__ == '__main__':
    from mesh import BARY, WQ
    rng = np.random.default_rng(0)
    p = np.array([[0.1, 0.0, 0.2], [1.0, 0.2, 0.1], [0.3, 0.9, 0.3]])
    n = np.cross(p[1]-p[0], p[2]-p[0]); A = 0.5*np.linalg.norm(n); n /= 2*A
    # Referenz: feine Unterteilung
    def ref(x, sub=60):
        acc_g = np.zeros(3); acc_i = 0
        for i in range(sub):
            for j in range(sub - i):
                for up in [0, 1]:
                    if up == 0: tri = np.array([[i, j], [i+1, j], [i, j+1]])/sub
                    elif j < sub - i - 1: tri = np.array([[i+1, j], [i+1, j+1], [i, j+1]])/sub
                    else: continue
                    for bq, wq in zip(BARY, WQ):
                        uv = bq@tri; y = p[0] + uv[0]*(p[1]-p[0]) + uv[1]*(p[2]-p[0]); d = x - y; r = np.linalg.norm(d)
                        acc_g += wq/sub**2*A*d/r**3; acc_i += wq/sub**2*A/r
        return acc_g, acc_i
    X = np.array([[0.5, 0.4, 1.0], [0.4, 0.3, 0.26], [1.5, -0.5, -0.3], [0.45, 0.35, 0.2 + 0.02]])
    G, Iv = tri_integrals(X, p, n)
    for k, x in enumerate(X):
        g, i = ref(x); print(f"x={x}: |dI_grad|/|I_grad| = {np.linalg.norm(G[k]-g)/np.linalg.norm(g):.1e},  dI_inv/I_inv = {abs(Iv[k]-i)/abs(i):.1e}")

def tri_integrals_multi(X, P, Nn):
    """Vektorisiert ueber Dreiecke: X (M,3), P (T,3,3), Nn (T,3) -> I_grad (T,M,3), I_inv (T,M)."""
    Xb = X[None, :, :]                                     # (1,M,3)
    w = np.einsum('tmd,td->tm', Xb - P[:, None, 0, :], Nn)  # (T,M)
    rho = Xb - w[..., None]*Nn[:, None, :]
    Ig = np.zeros(w.shape + (3,)); Iinv = np.zeros(w.shape); bsum = np.zeros(w.shape)
    for e in range(3):
        q1, q2 = P[:, e, :], P[:, (e+1) % 3, :]
        Lv = np.linalg.norm(q2 - q1, axis=1); t = (q2 - q1)/Lv[:, None]; mm = np.cross(t, Nn)
        lm = np.einsum('tmd,td->tm', q1[:, None, :] - rho, t); lp = np.einsum('tmd,td->tm', q2[:, None, :] - rho, t)
        P0 = np.einsum('tmd,td->tm', q1[:, None, :] - rho, mm)
        Rm = np.linalg.norm(Xb - q1[:, None, :], axis=2); Rp = np.linalg.norm(Xb - q2[:, None, :], axis=2)
        R0sq = P0**2 + w**2
        with np.errstate(divide='ignore', invalid='ignore'):
            f = np.where(lm + lp >= 0, np.log((Rp + lp)/(Rm + lm)), np.log((Rm - lm)/(Rp - lp)))
            aw = np.abs(w)
            beta = np.arctan2(P0*lp, R0sq + aw*Rp) - np.arctan2(P0*lm, R0sq + aw*Rm)
        Ig += f[..., None]*mm[:, None, :]; bsum += beta; Iinv += P0*f
    Iinv -= np.abs(w)*bsum
    return Ig + (np.sign(w)*bsum)[..., None]*Nn[:, None, :], Iinv
