"""Abbildungen fuer docs/results_coated.md (aus dem Wurzelverzeichnis aufrufen: python3 tools/plot_coated.py).

1. fig_coated_thin.png: Fehler der Schichtwirkung Delta sigma gegen d/h, Differenz zur Rechnung ohne Schicht bzw. zur
   neutralen Rechnung (results/coated_thin.csv, Goldkern mit Glasschale, omega a = 0,5).
2. fig_coated_agsphere.png: Silberkugel (20 nm) in Wasser mit 2 nm Oxid (n = 1,7): Mie gegen BEM, Extinktion und
   Phase der Vorwaertsamplitude.
3. fig_coated_agcube.png: Silberwuerfel (50 nm, Rundung 10 nm) mit 2 nm Oxid: BEM, Phasenaenderung.
Ausgabe der Tabellen auf stdout.
"""
import csv, sys, os
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, 'tools')
import mie_coated as mc
from mie_spectrum import eps_tab


def rd(fn): return list(csv.DictReader(open(fn)))
def col(rows, k): return np.array([float(r[k]) for r in rows])


# 1. Duennschichtstudie
rows = rd('results/coated_thin.csv')
core = complex(-11, 1.2); om = 0.5
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
bare = {int(r['n']): r for r in rows if r['coat'] == 'none'}
cur = None; res = {}
for r in rows:
    if r['coat'] == 'none': continue
    if r['coat'] == 'neutral': res[cur][int(r['n'])]['neutral'] = r; continue
    cur = r['coat']; res.setdefault(cur, {})[int(r['n'])] = {'coated': r}
for i, (c, byn) in enumerate(res.items()):
    d = float(c.split(',')[0]); ec = complex(float(c.split(',')[1]), float(c.split(',')[2]))
    S1 = mc.forward_amplitude(om, [1.0, 1.0 + d], [core, ec]); S0 = mc.forward_amplitude(om, [1.0], [core])
    ds = 4 * np.pi / om ** 2 * (S1 - S0).real; dp = np.angle(S1) - np.angle(S0)
    ns = sorted(n for n in byn if 'neutral' in byn[n] and n in bare)
    dh = np.array([d / (1.05 / n) for n in ns])
    eb = [abs((float(byn[n]['coated']['sigma_ext']) - float(bare[n]['sigma_ext'])) / ds - 1) for n in ns]
    en = [abs((float(byn[n]['coated']['sigma_ext']) - float(byn[n]['neutral']['sigma_ext'])) / ds - 1) for n in ns]
    pn = [abs((float(byn[n]['coated']['S_arg']) - float(byn[n]['neutral']['S_arg'])) / dp - 1) for n in ns]
    hs = [1.05 / n for n in ns]
    ax[0].loglog(hs, eb, 's--', color=f'C{i}', mfc='none', label=f'd = {d}: gegen ohne Schicht')
    ax[0].loglog(hs, en, 'o-', color=f'C{i}', label=f'd = {d}: gegen neutral')
    ax[1].loglog(hs, pn, 'o-', color=f'C{i}', label=f'd = {d}')
for a in ax:
    hh = np.array([0.08, 0.27]); a.set_xlabel('Elementgröße h (Kernradius 1)'); a.grid(True, which='both', alpha=0.3)
ax[0].loglog(hh, 0.3 * (hh / 0.27) ** 2, 'k:', lw=1, label='Ordnung 2')
ax[1].loglog(hh, 0.3 * (hh / 0.27) ** 2, 'k:', lw=1, label='Ordnung 2')
ax[0].set_ylabel(r'rel. Fehler von $\Delta\sigma_{ext}$'); ax[0].set_title('Extinktionsänderung durch die Schicht', fontsize=10)
ax[1].set_ylabel(r'rel. Fehler von $\Delta\arg S(0)$'); ax[1].set_title('Phasenänderung (Differenz zu neutral)', fontsize=10)
ax[0].legend(fontsize=6.5); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig('docs/fig_coated_thin.png', dpi=130)


# 2. Silberkugel
def spectrum_pair(fn):
    r = rd(fn); return col(r, 'lambda_nm'), col(r, 'sigma_nm2'), col(r, 'S_re') + 1j * col(r, 'S_im')
if os.path.exists('results/agsphere20_ox2_neutral.csv'):
    L, sb, Sb = spectrum_pair('results/agsphere20_bare.csv')
    _, sc, Sc = spectrum_pair('results/agsphere20_ox2.csv')
    _, sn, Sn = spectrum_pair('results/agsphere20_ox2_neutral.csv')
    a, t, nb = 20.0, 2.0, 1.33
    Lf = np.linspace(L.min(), L.max(), 201)
    def mie(lams, coat):
        out = []
        for l in lams:
            e = eps_tab('data/materials/Ag_Johnson.yml', l); o = 2 * np.pi / l
            if coat: S = mc.forward_amplitude(o, [a, a + t], [e, 2.89], nb ** 2)
            else: S = mc.forward_amplitude(o, [a], [e], nb ** 2)
            out.append(S)
        return np.array(out)
    k2 = lambda l: (2 * np.pi * nb / l) ** 2
    M0, M1 = mie(Lf, False), mie(Lf, True); m0, m1 = mie(L, False), mie(L, True)
    sig = lambda S, l: 4 * np.pi / k2(l) * S.real
    corr = sb + (sc - sn); pc = np.angle(Sb) + (np.angle(Sc) - np.angle(Sn))
    fig, ax = plt.subplots(2, 1, figsize=(5.8, 5.6), sharex=True)
    ax[0].plot(Lf, sig(M0, Lf) / 1e3, 'k-', lw=1, label='Mie, ohne Schicht')
    ax[0].plot(Lf, sig(M1, Lf) / 1e3, 'C3-', lw=1, label='Aden-Kerker, 2 nm Oxid')
    ax[0].plot(L, sb / 1e3, 'ko', mfc='none', ms=5, label='BEM ohne Schicht')
    ax[0].plot(L, corr / 1e3, 'C3s', ms=4, label='BEM: ohne + (beschichtet − neutral)')
    if os.path.exists('results/agsphere20_ox2_thin.csv'):
        _, st, St = spectrum_pair('results/agsphere20_ox2_thin.csv')
        ax[0].plot(L, st / 1e3, 'C2d', ms=4, mfc='none', label='BEM: Dünnschicht 1. Ordnung (Mitte)')
    if os.path.exists('results/agsphere20_ox2_thin2.csv'):
        _, s2, S2 = spectrum_pair('results/agsphere20_ox2_thin2.csv')
        ax[0].plot(L, s2 / 1e3, 'C1x', ms=5, label='BEM: Dünnschicht 2. Ordnung')
    ax[0].set_ylabel(r'$\sigma_{ext}$ (10$^3$ nm$^2$)'); ax[0].legend(fontsize=7); ax[0].set_title('Silberkugel, Radius 20 nm, Wasser; Oxid n = 1,7', fontsize=10)
    ax[1].plot(Lf, np.angle(M1) - np.angle(M0), 'C3-', lw=1, label='Aden-Kerker − Mie')
    ax[1].plot(L, np.angle(Sc) - np.angle(Sn), 'C3s', ms=4, label='BEM: beschichtet − neutral')
    ax[1].plot(L, np.angle(Sc) - np.angle(Sb), 'C0^', ms=4, mfc='none', label='BEM: beschichtet − ohne Schicht (hier d/h ≈ 0,8)')
    if os.path.exists('results/agsphere20_ox2_thin.csv'):
        ax[1].plot(L, np.angle(St) - np.angle(Sb), 'C2d', ms=4, mfc='none', label='BEM: Dünnschicht 1. Ordnung − ohne')
    if os.path.exists('results/agsphere20_ox2_thin2.csv'):
        ax[1].plot(L, np.angle(S2) - np.angle(Sb), 'C1x', ms=5, label='BEM: Dünnschicht 2. Ordnung − ohne')
    ax[1].axhline(0, color='k', lw=0.5); ax[1].set_ylabel(r'$\Delta\arg S(0)$ (rad)'); ax[1].set_xlabel('Wellenlänge (nm)'); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig('docs/fig_coated_agsphere.png', dpi=130)
    print('\nSilberkugel 20 nm, 2 nm Oxid (sigma in nm^2, Phase in rad)')
    print(f"{'lam':>5} {'Mie0':>8} {'BEM0':>8} {'Mie1':>8} {'korr.':>8} {'roh':>8} | {'dS_Mie':>9} {'dS_neut':>9} {'dS_ohne':>9} | {'dArg_Mie':>9} {'dArg_n':>9} {'dArg_o':>9}")
    if os.path.exists('results/agsphere20_ox2_thin.csv'):
        print('Duennschicht 1. Ordnung (Mitte) / 2. Ordnung: lam, sigma1, sigma2, sigma Aden-Kerker, korr. zwei Flaechen | dArg1, dArg2, dArg exakt, dArg zwei Fl.')
        s2p = s2 if os.path.exists('results/agsphere20_ox2_thin2.csv') else st; S2p = S2 if os.path.exists('results/agsphere20_ox2_thin2.csv') else St
        for i, l in enumerate(L):
            print(f"  {l:5.0f} {st[i]:8.0f} {s2p[i]:8.0f} {sig(m1[i], l):8.0f} {corr[i]:8.0f} | {np.angle(St[i]) - np.angle(Sb[i]):+8.4f} {np.angle(S2p[i]) - np.angle(Sb[i]):+8.4f}"
                  f" {np.angle(m1[i]) - np.angle(m0[i]):+8.4f} {np.angle(Sc[i]) - np.angle(Sn[i]):+8.4f}")
    for i, l in enumerate(L):
        print(f"{l:5.0f} {sig(m0[i], l):8.0f} {sb[i]:8.0f} {sig(m1[i], l):8.0f} {corr[i]:8.0f} {sc[i]:8.0f} | {sig(m1[i], l) - sig(m0[i], l):+9.0f} {sc[i] - sn[i]:+9.0f} {sc[i] - sb[i]:+9.0f} |"
              f" {np.angle(m1[i]) - np.angle(m0[i]):+9.4f} {np.angle(Sc[i]) - np.angle(Sn[i]):+9.4f} {np.angle(Sc[i]) - np.angle(Sb[i]):+9.4f}")

# 3. Silberwuerfel
if os.path.exists('results/agcube_rho0.4_ox2_neutral.csv'):
    rb = sorted(rd('results/agcube_rho0.4.csv'), key=lambda r: float(r['lambda_nm']))
    L0 = col(rb, 'lambda_nm'); s0 = col(rb, 'sigma_nm2')
    L, sc, Sc = spectrum_pair('results/agcube_rho0.4_ox2.csv'); _, sn, Sn = spectrum_pair('results/agcube_rho0.4_ox2_neutral.csv')
    s0i = np.interp(L, L0, s0); corr = s0i + (sc - sn)
    fig, ax = plt.subplots(2, 1, figsize=(5.8, 5.6), sharex=True)
    ax[0].plot(L, s0i / 1e3, 'ko-', mfc='none', label='ohne Schicht (v0.12)')
    ax[0].plot(L, corr / 1e3, 'C3s-', label='2 nm Oxid: ohne + (beschichtet − neutral)')
    ax[0].set_ylabel(r'$\sigma_{ext}$ (10$^3$ nm$^2$)'); ax[0].legend(fontsize=7)
    ax[0].set_title('Silberwürfel 50 nm, Rundung 10 nm, Wasser; Oxid n = 1,7', fontsize=10)
    ax[1].plot(L, np.angle(Sc) - np.angle(Sn), 'C3s-', label='beschichtet − neutral'); ax[1].axhline(0, color='k', lw=0.5)
    if os.path.exists('results/agcube_rho0.4_ox2_thin.csv'):
        _, sth, Sth = spectrum_pair('results/agcube_rho0.4_ox2_thin.csv')
        ax[0].plot(L, sth / 1e3, 'C1x--', label='2 nm Oxid: Dünnschicht 2. Ordnung')
        ax[0].legend(fontsize=7)
        S0c = np.array([complex(float(r['S_re']), float(r['S_im'])) for r in rb]) if 'S_re' in rb[0] else None
    ax[1].set_ylabel(r'$\Delta\arg S(0)$ (rad)'); ax[1].set_xlabel('Wellenlänge (nm)'); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig('docs/fig_coated_agcube.png', dpi=130)
    if os.path.exists('results/agcube_rho0.4_ox2_thin.csv'):
        print('\nSilberwuerfel: Schichtwirkung zwei Flaechen (gegen neutral) und Duennschicht 2. Ordnung (gegen v0.12 ohne Schicht)')
        for i, l in enumerate(L):
            d2 = sc[i] - sn[i]; dt = sth[i] - s0i[i]
            print(f"  {l:5.0f}  ohne {s0i[i]:8.0f}  zwei Fl. {corr[i]:8.0f} (Delta {d2:+7.0f})  Duennschicht {sth[i]:8.0f} (Delta {dt:+7.0f}, {dt / d2 - 1:+6.1%})")
    print('\nSilberwuerfel 50 nm, Rundung 10 nm, 2 nm Oxid')
    print(f"{'lam':>5} {'ohne':>8} {'mit Oxid':>9} {'Delta':>8} {'dArg':>9} {'arg S neutral':>13}")
    for i, l in enumerate(L):
        print(f"{l:5.0f} {s0i[i]:8.0f} {corr[i]:9.0f} {sc[i] - sn[i]:+8.0f} {np.angle(Sc[i]) - np.angle(Sn[i]):+9.4f} {np.angle(Sn[i]):13.4f}")

# 4. Duennschicht-Naeherung 1. und 2. Ordnung: Fehler der Schichtwirkung gegen d (feinstes Netz je Variante)
def thin_series(fn, want_model):
    rows = rd(fn); bare = {int(r['n']): r for r in rows if r['coat'] == 'none'}
    S0 = mc.forward_amplitude(om, [1.0], [core]); series = {}
    for r in rows:
        c = r['coat']
        if not c.startswith('thin') or r['inward'] == '1': continue
        var, spec = c.split(':', 1)
        if ';' in spec: continue
        d, er, ei = [float(v) for v in spec.split(',')]
        if (er, ei) != (2.25, 0.0): continue
        model = var.split('@')[0].replace('thin', '').lstrip('-') or 'jump'
        if model != want_model: continue
        f = float(var.split('@')[1]) if '@' in var else 0.0; n = int(r['n'])
        if n not in bare: continue
        S1 = mc.forward_amplitude(om, [1.0, 1.0 + d], [core, 2.25])
        ds = 4 * np.pi / om ** 2 * (S1 - S0).real; dp = np.angle(S1) - np.angle(S0)
        es = (float(r['sigma_ext']) - float(bare[n]['sigma_ext'])) / ds - 1; ep = (float(r['S_arg']) - float(bare[n]['S_arg'])) / dp - 1
        key = (f, d)
        if key not in series or series[key][0] < n: series[key] = (n, es, ep)
    return series
if os.path.exists('results/thin_layer.csv'):
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.8))
    curves = [('results/thin_layer.csv', 'jump', 0.0, 'o-', 'C0', '1. Ordnung, Referenz Kernoberfläche'),
              ('results/thin_layer.csv', 'jump', 0.5, 's-', 'C1', '1. Ordnung, Referenz Schichtmitte'),
              ('results/thin_layer2.csv', 'dirac2', 0.0, 'D-', 'C2', '2. Ordnung (Dirac-Form), Referenz Kernoberfläche')]
    for fn, mdl, f, st, colr, lab in curves:
        if not os.path.exists(fn): continue
        ser = thin_series(fn, mdl); ks = sorted(k for k in ser if k[0] == f)
        ds = [k[1] for k in ks]
        ax[0].loglog(ds, [max(abs(ser[k][1]), 1e-4) for k in ks], st, color=colr, label=lab)
        ax[1].loglog(ds, [max(abs(ser[k][2]), 1e-4) for k in ks], st, color=colr, label=lab)
    if os.path.exists('results/coated_thin.csv'):
        tw = list(csv.DictReader(open('results/coated_thin.csv'))); pts = []
        S0 = mc.forward_amplitude(om, [1.0], [core])
        for i, r in enumerate(tw):
            if r['coat'] in ('none', 'neutral') or int(r['n']) != 12: continue
            d = float(r['coat'].split(',')[0]); nr = tw[i + 1]
            S1 = mc.forward_amplitude(om, [1.0, 1.0 + d], [core, 2.25]); ds_ = 4 * np.pi / om ** 2 * (S1 - S0).real; dp_ = np.angle(S1) - np.angle(S0)
            pts.append((d, abs((float(r['sigma_ext']) - float(nr['sigma_ext'])) / ds_ - 1), abs((float(r['S_arg']) - float(nr['S_arg'])) / dp_ - 1)))
        pts.sort()
        for j in range(2):
            ax[j].loglog([p[0] for p in pts], [p[1 + j] for p in pts], '^--', color='C3', mfc='none', label='zwei Flächen, gegen neutral (n = 12)')
    dd = np.array([0.005, 0.2])
    for a in ax:
        a.loglog(dd, 1.4 * dd, 'k:', lw=1, label='1,4 d/a'); a.loglog(dd, 4 * dd ** 2, 'k--', lw=0.8, label='4 (d/a)²')
        a.set_xlabel('Schichtdicke d (Kernradius a = 1)'); a.grid(True, which='both', alpha=0.3); a.set_ylim(5e-4, 0.5)
    ax[0].set_ylabel(r'rel. Fehler von $\Delta\sigma_{ext}$'); ax[1].set_ylabel(r'rel. Fehler von $\Delta\arg S(0)$')
    ax[0].set_title('Extinktionsänderung (feinstes Netz)', fontsize=10); ax[1].set_title('Phasenänderung (feinstes Netz)', fontsize=10)
    ax[0].legend(fontsize=6, loc='upper left'); fig.tight_layout(); fig.savefig('docs/fig_coated_thinlayer.png', dpi=130)

# 5. Chirale Schichten: Fehler des CD gegen das Netz (results/thin_chiral.csv, tools/analyze_chiral_thin.py)
if os.path.exists('results/thin_chiral.csv'):
    import mie_chiral_layered as mcl
    from collections import defaultdict
    g = defaultdict(dict)
    for r in rd('results/thin_chiral.csv'):
        if float(r['core_re']) != -11.0: continue
        g[(int(r['n']), r['coat'])][int(r['pol'])] = float(r['sigma_ext'])
    gold = (-11 + 1.2j, 1.0, 0.0)
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    for j, d in enumerate([0.02, 0.05, 0.1]):
        m = mcl.cross_sections(om, [1, 1 + d], [gold, (2.25, 1.0, 0.1)], +1)[0] - mcl.cross_sections(om, [1, 1 + d], [gold, (2.25, 1.0, 0.1)], -1)[0]
        for meth, st in [('thin-dirac2fit:', 'o-'), ('', 's--')]:
            key = f"{meth}{d},2.25,0,0.1"
            pts = sorted((n, abs((v[+1] - v[-1]) / m - 1)) for (n, c), v in g.items() if c == key and +1 in v and -1 in v)
            if not pts: continue
            lab = f"d = {d}: {'Dünnschicht 2. Ordnung' if meth else 'zwei Flächen'}"
            ax.loglog([1.05 / p[0] for p in pts], [max(p[1], 1e-4) for p in pts], st, color=f'C{j}', mfc='none' if not meth else None, label=lab)
    hh = np.array([0.08, 0.27]); ax.loglog(hh, 0.1 * (hh / 0.27) ** 2, 'k:', lw=1, label='Ordnung 2')
    ax.set_xlabel('Elementgröße h (Kernradius 1)'); ax.set_ylabel('rel. Fehler des CD'); ax.grid(True, which='both', alpha=0.3)
    ax.set_title('Goldkern mit chiraler Schale (χ = 0,1), CD gegen Mie', fontsize=10); ax.legend(fontsize=6.5)
    fig.tight_layout(); fig.savefig('docs/fig_coated_chiral.png', dpi=130)

# 6. Goldkugel (20 nm) in Wasser mit 1 nm chiraler Huelle: Extinktion und CD gegen Mie (results/ausphere20_chiral1.csv)
if os.path.exists('results/ausphere20_chiral1.csv'):
    import mie_chiral_layered as mcl
    cf = lambda rr, k: np.array([float(r[k]) for r in rr])          # 'col' ist oben als Farbvariable belegt
    rows = rd('results/ausphere20_chiral1.csv'); L = cf(rows, 'lambda_nm'); sp = cf(rows, 'sigma_plus_nm2'); cd = cf(rows, 'CD_nm2')
    nb, a, t, chi = 1.33, 20.0, 1.0, 0.01
    Lf = np.linspace(L.min(), L.max(), 101)
    def mie_pm(lams):
        P, M = [], []
        for l in lams:
            e = eps_tab('data/materials/Au_Johnson.yml', l); o = 2 * np.pi / l
            P.append(mcl.cross_sections(o, [a, a + t], [(e, 1, 0), (2.25, 1, chi)], +1, eps2=nb ** 2)[0])
            M.append(mcl.cross_sections(o, [a, a + t], [(e, 1, 0), (2.25, 1, chi)], -1, eps2=nb ** 2)[0])
        return np.array(P), np.array(M)
    Pf, Mf = mie_pm(Lf); Pm, Mm = mie_pm(L)
    fig, ax = plt.subplots(2, 1, figsize=(5.8, 5.4), sharex=True)
    ax[0].plot(Lf, Pf / 1e3, 'k-', lw=1, label='Mie (geschichtet, chiral)'); ax[0].plot(L, sp / 1e3, 'C3o', ms=4, mfc='none', label='BEM Dünnschicht 2. Ordnung')
    ax[0].set_ylabel(r'$\sigma_{ext}^{+}$ (10$^3$ nm$^2$)'); ax[0].legend(fontsize=7)
    ax[0].set_title('Goldkugel 20 nm in Wasser, 1 nm chirale Hülle (n = 1,5, χ = 0,01)', fontsize=10)
    ax[1].plot(Lf, Pf - Mf, 'k-', lw=1, label='Mie'); ax[1].plot(L, cd, 'C3o', ms=4, mfc='none', label='BEM'); ax[1].axhline(0, color='k', lw=0.5)
    ax[1].set_ylabel(r'CD = $\sigma^+ - \sigma^-$ (nm$^2$)'); ax[1].set_xlabel('Wellenlänge (nm)'); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig('docs/fig_coated_auchiral.png', dpi=130)
    print('\nGoldkugel mit chiraler Huelle: lam, sigma+ BEM, Mie, CD BEM, CD Mie, Fehler, g = CD/sigma (Mie)')
    for i, l in enumerate(L):
        print(f"  {l:5.0f} {sp[i]:8.1f} {Pm[i]:8.1f} {cd[i]:+10.4e} {Pm[i] - Mm[i]:+10.4e} {cd[i] / (Pm[i] - Mm[i]) - 1:+7.1%} {(Pm[i] - Mm[i]) / Pm[i]:+.2e}")

# 7. Zweitor (S-Matrix-Formulierung): absoluter Fehler von sigma gegen d und CD der chiralen Huelle (results/twoport*.csv)
if os.path.exists('results/twoport.csv'):
    from collections import defaultdict
    import mie_chiral_layered as mcl
    rows = rd('results/twoport.csv'); S0 = mc.forward_amplitude(om, [1.0], [core])
    bare = {int(r['n']): float(r['sigma_ext']) for r in rows if r['coat'] == 'none'}
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.8))
    for meth, pre, st in [('Zweitor', 'twoport:', 'o-'), ('Dünnschicht 2. Ordnung', 'thin-dirac2fit:', 's--')]:
        for j, n in enumerate([8, 12]):
            pts = []
            for r in rows:
                if int(r['n']) != n or not r['coat'].startswith(pre): continue
                dd = float(r['coat'].split(':')[1].split(',')[0])
                sx = 16 * np.pi * mc.forward_amplitude(om, [1, 1 + dd], [core, 2.25]).real
                pts.append((dd, abs(float(r['sigma_ext']) / sx - 1)))
            if pts:
                pts.sort(); ax[0].loglog([p[0] for p in pts], [p[1] for p in pts], st, color=f'C{j}', mfc='none' if pre.startswith('thin') else None, label=f'{meth}, n = {n}')
    for j, n in enumerate([8, 12]):
        if n in bare: ax[0].axhline(abs(bare[n] / (16 * np.pi * S0.real) - 1), color=f'C{j}', lw=0.8, ls=':', label=f'ohne Schicht, n = {n}')
    ax[0].set_xlabel('Schichtdicke d (Kernradius a = 1)'); ax[0].set_ylabel(r'rel. Fehler von $\sigma_{ext}$ (absolut)')
    ax[0].set_title('Goldkern mit Glasschale', fontsize=10); ax[0].grid(True, which='both', alpha=0.3); ax[0].legend(fontsize=6.5)
    g = defaultdict(dict)
    for r in rd('results/twoport_chiral.csv'): g[(int(r['n']), r['coat'])][int(r['pol'])] = float(r['sigma_ext'])
    gold = (core, 1.0, 0.0)
    for j, dd in enumerate([0.02, 0.05, 0.1, 0.2]):
        m = mcl.cross_sections(om, [1, 1 + dd], [gold, (2.25, 1.0, 0.1)], +1)[0] - mcl.cross_sections(om, [1, 1 + dd], [gold, (2.25, 1.0, 0.1)], -1)[0]
        pts = sorted((n, abs((v[+1] - v[-1]) / m - 1)) for (n, c), v in g.items() if c == f'twoport:{dd},2.25,0,0.1' and +1 in v and -1 in v)
        ax[1].loglog([1.05 / p[0] for p in pts], [p[1] for p in pts], 'o-', color=f'C{j}', label=f'd = {dd}')
    hh = np.array([0.08, 0.18]); ax[1].loglog(hh, 0.04 * (hh / 0.175) ** 2, 'k:', lw=1, label='Ordnung 2')
    ax[1].set_xlabel('Elementgröße h'); ax[1].set_ylabel('rel. Fehler des CD'); ax[1].grid(True, which='both', alpha=0.3)
    ax[1].set_title('Zweitor: chirale Hülle (χ = 0,1)', fontsize=10); ax[1].legend(fontsize=6.5)
    fig.tight_layout(); fig.savefig('docs/fig_coated_twoport.png', dpi=130)

# 8. Plasmoninduzierter CD gegen den Abstand der chiralen Schicht vom Gold (results/au_spacer_cd.csv, Zweitor)
if os.path.exists('results/au_spacer_cd.csv'):
    import mie_chiral_layered as mcl
    rows = rd('results/au_spacer_cd.csv')
    spacers = []
    for r in rows:
        c = r['coating'].split(' ')[0]; parts = c.split(';')
        spacers.append(float(parts[0].split(':')[0]) if len(parts) > 1 else 0.0)
    cdb = np.array([float(r['CD_nm2']) for r in rows]); spb = np.array([float(r['sigma_plus_nm2']) for r in rows])
    L = 530.0; e = eps_tab('data/materials/Au_Johnson.yml', L); o = 2 * np.pi / L; nb2 = 1.33 ** 2
    def mie_cd(s):
        R = [20.0] + ([20.0 + s] if s > 0 else []) + [21.0 + s]; M = [(e, 1, 0)] + ([(2.1025, 1, 0)] if s > 0 else []) + [(2.25, 1, 0.01)]
        p = mcl.cross_sections(o, R, M, +1, eps2=nb2)[0]; q = mcl.cross_sections(o, R, M, -1, eps2=nb2)[0]; return p - q, p
    sf = np.linspace(0, 10, 41); cdf = np.array([mie_cd(s)[0] for s in sf])
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.plot(sf, cdf, 'k-', lw=1, label='Mie (geschichtet, chiral)'); ax.plot(spacers, cdb, 'C3o', label='BEM Zweitor (n = 8)')
    ax.axhline(0, color='k', lw=0.5); ax.set_xlabel('Dicke des Abstandshalters (nm, n = 1,45)'); ax.set_ylabel(r'CD bei 530 nm (nm$^2$)')
    ax.set_title('Goldkugel 20 nm, 1 nm chirale Schicht (χ = 0,01), Wasser', fontsize=10); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig('docs/fig_coated_spacer.png', dpi=130)
    print('\nAbstandshalter: s, CD BEM, CD Mie, Fehler absolut, sigma+ BEM, Mie, Fehler, Iterationen')
    for s, c, sp, r in zip(spacers, cdb, spb, rows):
        m, p = mie_cd(s); print(f"  {s:4.0f} {c:+.5f} {m:+.5f} {c - m:+.4f} {sp:8.1f} {p:8.1f} {sp / p - 1:+6.1%} {r['iterations']}")
