"""Stats-model check: recompute the DP4 probabilities of every original .dp4 file from its own printed tables
(calculated and experimental shift per atom, pairing as printed) and compare with the printed probabilities.
Usage: check_model.py <res/core>"""
import sys, os, re, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dp4lib as L, dp4core as D
S = L.Small(sys.argv[1], 'dp4ai')
n = 0; worst = collections.defaultdict(float); bad = []; ncomp = collections.Counter()
for p in sorted(x for x in S.paths() if x.endswith('.dp4')):
    d = L.parse_dp4(S.read(p)); n += 1
    logs = {}
    for nuc, sigma in (('C', D.SIGMA_C), ('H', D.SIGMA_H)):
        ll = []
        for k, tabs in d['isomers'].items():
            rows = tabs.get(nuc, [])
            if len(rows) < 3: ll = None; break
            calc = [r[1] for r in rows]; exp = [r[3] for r in rows]
            sc = D.scale(calc, exp)
            ll.append(D.log_prob([s - e for s, e in zip(sc, exp)], sigma))
        if ll: logs[nuc] = ll
    out = {k: D.normalise(v) for k, v in logs.items()}
    if logs: out['all'] = D.normalise([sum(logs[k][i] for k in logs) for i in range(len(next(iter(logs.values()))))])
    for key in ('C', 'H', 'all'):
        if key in out and key in d['probs'] and len(out[key]) == len(d['probs'][key]):
            dev = max(abs(100 * a - b) for a, b in zip(out[key], d['probs'][key]))
            worst[key] = max(worst[key], dev); ncomp[key] += 1
            if dev > 1.0: bad.append((p.split('!/')[1], key, round(dev, 1)))
        else:
            ncomp[key + ' not comparable'] += 1
print('files: %d | compared (C, H, combined): %s' % (n, dict(ncomp)))
print('largest difference to a printed probability, in percentage points:', {k: round(v, 2) for k, v in worst.items()})
print('files with a difference above 1 point: %d' % len(set(b[0] for b in bad)), bad[:12])
