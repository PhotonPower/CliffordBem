"""Auswertung der Duennschicht-Naeherung erster Ordnung (scatter_coated --thin, results/thin_layer.csv).

Schichtwirkung Delta = (mit Schicht) - (ohne Schicht, gleiches Netz) fuer sigma_ext und arg S(0), verglichen mit
  exakt:        Aden-Kerker (tools/mie_coated.py),
  linearisiert: Aden-Kerker in erster Ordnung in d (numerische Ableitung bei d -> 0); der Unterschied zu exakt ist der
                Abbruchfehler O(d^2) der Naeherung, der Unterschied BEM - linearisiert der Diskretisierungsfehler,
und mit der exakten Zwei-Flaechen-Rechnung (Differenz zur neutralen Rechnung, results/coated_thin.csv), soweit vorhanden.
Aufruf: python3 tools/analyze_thin.py [results/thin_layer.csv]
"""
import csv, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_coated as mc

fn = sys.argv[1] if len(sys.argv) > 1 else 'results/thin_layer.csv'
rows = list(csv.DictReader(open(fn)))
om = 0.5; core = complex(-11, 1.2)


def mie_delta(layers, inward=False, lin=False):
    """Delta sigma, Delta arg S gegen die Kugel mit Radius 1 (ohne Schicht); layers: [(d, eps)] von innen nach aussen."""
    S0 = mc.forward_amplitude(om, [1.0], [core])
    T = sum(d for d, _ in layers)
    def S(scale):
        if scale == 0: return S0
        r = 1.0 - scale * T if inward else 1.0; radii = [r]; eps = [core]
        for d, e in layers: r += scale * d; radii.append(r); eps.append(e)
        return mc.forward_amplitude(om, radii, eps)
    if lin:
        e = 1e-6; S1 = S0 + (S(e) - S0) / e
    else: S1 = S(1.0)
    return 4 * np.pi / om ** 2 * (S1 - S0).real, np.angle(S1) - np.angle(S0)


# Zwei-Flaechen-Ergebnisse (Differenz zur neutralen Rechnung)
two = {}
if os.path.exists('results/coated_thin.csv'):
    last = None
    for r in csv.DictReader(open('results/coated_thin.csv')):
        n = int(r['n'])
        if r['coat'] not in ('none', 'neutral'): last = (n, float(r['coat'].split(',')[0])); two[last] = {'c': r}
        elif r['coat'] == 'neutral' and last: two[last]['n'] = r

bare = {}; out = []
for r in rows:
    n = int(r['n']); c = r['coat']
    if c == 'none': bare[n] = r; continue
    var, spec = c.split(':', 1)
    layers = []
    for part in spec.split(';'):
        v = [float(x) for x in part.split(',')]; layers.append((v[0], complex(v[1], v[2] if len(v) > 2 else 0.0)))
    d = sum(l[0] for l in layers); ec = layers[0][1] if len(layers) == 1 else tuple(l[1] for l in layers)
    ref = float(var.split('@')[1]) if '@' in var else (1.0 if r['inward'] == '1' else 0.0)   # Anteil der Schicht innen
    model = var.split('@')[0].replace('thin', '').lstrip('-') or 'jump'
    out.append((n, d, ec, r['inward'] == '1', (model, ref), r, layers))

print(f"{'n':>3} {'d':>6} {'d/h':>5} | {'dSig BEM':>9} {'exakt':>9} {'linear':>9} {'Fehl.ex':>8} {'Fehl.lin':>8} {'2-Fl.':>7} |"
      f" {'dArg BEM':>9} {'exakt':>9} {'linear':>9} {'Fehl.ex':>8} {'Fehl.lin':>8} {'2-Fl.':>7} | It.")
prev = None
for n, d, ec, inw, mr, r, layers in sorted(out, key=lambda t: (t[3], str(t[2]), t[4], -t[1], t[0])):
    if n not in bare: continue
    key = (inw, str(ec), mr)
    if key != prev:
        print(f"-- Schicht eps = {ec}, {'nach innen' if inw else 'nach aussen'}, Modell {mr[0]}, Referenzflaeche {mr[1]:.2f} der Schicht von innen"); prev = key
    b = bare[n]
    ds = float(r['sigma_ext']) - float(b['sigma_ext']); dp = float(r['S_arg']) - float(b['S_arg'])
    es, ep_ = mie_delta(layers, inw); ls, lp = mie_delta(layers, inw, lin=True)
    t2s = t2p = ''
    if not inw and ec == 2.25 and mr[0] == 'jump' and (n, d) in two and 'n' in two[(n, d)]:
        c2, n2 = two[(n, d)]['c'], two[(n, d)]['n']
        t2s = f"{(float(c2['sigma_ext']) - float(n2['sigma_ext'])) / es - 1:+7.1%}"; t2p = f"{(float(c2['S_arg']) - float(n2['S_arg'])) / ep_ - 1:+7.1%}"
    print(f"{n:3d} {d:6.3f} {d / (1.05 / n):5.2f} | {ds:+9.5f} {es:+9.5f} {ls:+9.5f} {ds / es - 1:+8.1%} {ds / ls - 1:+8.1%} {t2s:>7} |"
          f" {dp:+9.5f} {ep_:+9.5f} {lp:+9.5f} {dp / ep_ - 1:+8.1%} {dp / lp - 1:+8.1%} {t2p:>7} | {r['iterations']}")
