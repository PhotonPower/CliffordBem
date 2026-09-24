"""Auswertung der Kompressions-Benchmarks: Wachstumsexponenten und Vergleich mit N log^a N."""
import csv, math, sys
fn = sys.argv[1] if len(sys.argv) > 1 else 'results/results_sphere.csv'
rows = list(csv.DictReader(open(fn)))
rows.sort(key=lambda r: int(r['N']))
print(f"{'N':>7} {'MB ges.':>9} {'MB NR':>9} {'MB/(N logN)':>12} {'MB/(N log^2N)':>14} {'Aufbau s':>9} {'Fehler':>8}")
for r in rows:
    N = int(r['N']); tot = float(r['MB_total']); lr = float(r['MB_lowrank'])
    print(f"{N:7d} {tot:9.1f} {lr:9.1f} {1e3*tot/(N*math.log(N)):12.4f} {1e3*tot/(N*math.log(N)**2):14.5f} {float(r['t_build_s']):9.1f} {float(r['rel_err']):8.1e}")
print("\nlokale Exponenten p (Speicher ~ N^p) zwischen aufeinanderfolgenden Netzen:")
for a, b in zip(rows[:-1], rows[1:]):
    Na, Nb = int(a['N']), int(b['N'])
    p = math.log(float(b['MB_total'])/float(a['MB_total']))/math.log(Nb/Na)
    pt = math.log(float(b['t_build_s'])/float(a['t_build_s']))/math.log(Nb/Na)
    print(f"  {Na:6d} -> {Nb:6d}: Speicher p = {p:.2f}, Aufbauzeit p = {pt:.2f}")
