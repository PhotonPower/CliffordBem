"""Teil III: gewichteter Energiefluss
S = E x conj(H) + conj(alpha) E / conj(sqrt mu) + beta conj(H) / sqrt(eps),
mit E = <G>_1/sqrt(eps), H = -I<G>_2/sqrt(mu), alpha=<G>_0, beta = -I<G>_3 (Koeffizient von I).
Behauptung 1: Re div S = -w [eps''|E|^2 + mu''|H|^2 + sin(arg mu)|alpha|^2 + sin(arg eps)|beta|^2]
Behauptung 2: n.S stetig ueber Gamma, falls h1 = J h2 mit a = conj(r_eps) sqrt(r_mu), b = conj(r_mu) sqrt(r_eps)."""
import numpy as np
from fields import gp, dirac_kernel, blade, vec, PSEUDO
def parts(G, eps, mu):
    alpha = G[..., 0]; beta = G[..., 7]              # G_3 = beta * e123 = beta I
    E = np.stack([G[..., 1], G[..., 2], G[..., 4]], -1)/np.sqrt(eps)
    # grade 2 = I b: e23 = I e1, e13 = -I e2, e12 = I e3
    b = np.stack([G[..., 6], -G[..., 5], G[..., 3]], -1); H = b/np.sqrt(mu)
    return alpha, E, H, beta
def S(G, eps, mu):
    al, E, H, be = parts(G, eps, mu)
    return np.cross(E, np.conj(H)) + (np.conj(al)/np.conj(np.sqrt(mu)))[..., None]*E + (be/np.sqrt(eps))[..., None]*np.conj(H)
rng = np.random.default_rng(7); c = rng.normal(size=8) + 1j*rng.normal(size=8)
om = 1.3
print("Behauptung 1 (Divergenz, finite Differenzen):")
for eps, mu in [(2.25, 1.0), (-11+1.2j, 1.0), (3+0.4j, 1.5+0.2j), (-4+0j, 1.0)]:
    k = om*np.sqrt(eps)*np.sqrt(mu)
    G = lambda x: gp(dirac_kernel(x, np.array([0.1, -0.2, 0.05]), k), c)
    x = np.array([[0.9, 0.4, -0.6]]); h = 1e-5
    div = sum((S(G(x+h*np.eye(3)[j]), eps, mu)[0, j] - S(G(x-h*np.eye(3)[j]), eps, mu)[0, j])/(2*h) for j in range(3))
    al, E, H, be = parts(G(x), eps, mu)
    rhs = -om*(np.imag(eps)*np.sum(abs(E)**2) + np.imag(mu)*np.sum(abs(H)**2)
               + np.sin(np.angle(mu))*abs(al[0])**2 + np.sin(np.angle(eps))*abs(be[0])**2)
    print(f"  eps={eps}, mu={mu}:  Re div S = {div.real:+.8f},  Vorhersage = {rhs:+.8f}")
print("Behauptung 2 (Stetigkeit von n.S unter J):")
n = rng.normal(size=3); n /= np.linalg.norm(n)
for (e1, m1, e2, m2) in [(2.25, 1, 1, 1), (-11+1.2j, 1, 1.77, 1), (3+0.4j, 1.5+0.2j, 1.2, 0.8)]:
    se, sm = np.sqrt(e1)/np.sqrt(e2), np.sqrt(m1)/np.sqrt(m2); re_, rm = e1/e2, m1/m2
    a, b = np.conj(re_)*sm, np.conj(rm)*se
    h2 = rng.normal(size=8) + 1j*rng.normal(size=8)
    v = np.array([h2[1], h2[2], h2[4]]); bb = np.array([h2[6], -h2[5], h2[3]])
    vn, bn = v@n, bb@n; vT, bT = v - vn*n, bb - bn*n
    v1 = se*vT + vn/se*n; b1 = sm*bT + bn/sm*n
    h1 = np.zeros(8, complex); h1[0] = a*h2[0]; h1[7] = b*h2[7]
    h1[1], h1[2], h1[4] = v1; h1[6], h1[5], h1[3] = b1[0], -b1[1], b1[2]
    s1, s2 = S(h1, e1, m1)@n, S(h2, e2, m2)@n
    print(f"  eps1={e1}, mu1={m1}: n.S innen = {s1:.10f}, aussen = {s2:.10f}")
