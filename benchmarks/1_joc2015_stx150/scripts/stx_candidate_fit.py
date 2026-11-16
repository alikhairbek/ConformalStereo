"""Plausibility of every candidate printed in the STX150 SI (validation set): error of the candidate against the
experimental shifts of the case after linear scaling (RMS and largest, 13C and 1H), at each of the six levels.
A stereoisomer that is a real competitor fits within a few ppm; a candidate with errors of tens of ppm is not
the same molecule with the same atom order, and the case tests nothing at that level.
Usage: stx_candidate_fit.py <stx out dir> <out csv>"""
import sys, os, csv, math, collections
stx_out, out_csv = sys.argv[1:3]
LEVELS = ['MMFF_PBE', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9', 'B3L_PBE', 'B3L_wB9']
blocks = collections.OrderedDict()
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_val_blocks.csv'))): blocks.setdefault(r['block'], []).append(r)
cases = list(csv.DictReader(open(os.path.join(stx_out, 'stx_val_cases.csv'))))
def fit(exp, sig):
    n = len(exp)
    if n < 3: return None, None
    mx, my = sum(sig) / n, sum(exp) / n; sxx = sum((s - mx) ** 2 for s in sig); sxy = sum((s - mx) * (e - my) for s, e in zip(sig, exp))
    if sxx == 0: return None, None
    b = sxy / sxx; a = my - b * mx; err = [e - (a + b * s) for s, e in zip(sig, exp)]
    return math.sqrt(sum(x * x for x in err) / n), max(abs(x) for x in err)
out = []
for c in cases:
    case = c['case']; rows = blocks[case]
    if not any(r['exp'] != '' for r in rows): continue
    for cand in c['candidates'].split():
        if len(blocks[cand]) != len(rows): continue
        for lev in LEVELS:
            rec = dict(case=case, candidate=cand, is_correct=int(cand == case), level=lev)
            for nuc in 'CH':
                pairs = [(float(r['exp']), float(q[lev])) for r, q in zip(rows, blocks[cand]) if r['nucleus'] == nuc and r['exp'] != '' and q[lev] != '']
                rms, mx = fit([p[0] for p in pairs], [p[1] for p in pairs])
                rec['n_' + nuc] = len(pairs); rec['rms_' + nuc] = '' if rms is None else round(rms, 3); rec['max_' + nuc] = '' if mx is None else round(mx, 3)
            out.append(rec)
w = csv.DictWriter(open(out_csv, 'w', newline=''), fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
f = lambda v: None if v == '' else float(v)
def q(v, p):
    v = sorted(v); return v[min(len(v) - 1, int(p * len(v)))]
for grp, sel in (('correct isomers', [o for o in out if o['is_correct']]), ('alternatives', [o for o in out if not o['is_correct']])):
    c = [f(o['rms_C']) for o in sel if f(o['rms_C']) is not None]; h = [f(o['rms_H']) for o in sel if f(o['rms_H']) is not None]
    print('%-16s candidate-level entries %4d | 13C RMS: median %.2f, 90%% %.2f, 99%% %.2f, largest %.2f | 1H RMS: median %.3f, 90%% %.3f, 99%% %.3f, largest %.3f' % (grp, len(sel), q(c, .5), q(c, .9), q(c, .99), max(c), q(h, .5), q(h, .9), q(h, .99), max(h)))
bad = collections.OrderedDict()
for o in out:
    c, h = f(o['rms_C']), f(o['rms_H'])
    if (c is not None and c > 6.0) or (h is not None and h > 0.6): bad.setdefault((o['case'], o['candidate']), []).append((o['level'], c, h))
print('candidates with a 13C RMS error above 6 ppm or a 1H RMS error above 0.6 ppm at some level: %d' % len(bad))
for (case, cand), v in bad.items():
    print('  case %-5s candidate %-10s %s' % (case, cand + (' (correct)' if cand == case else ''), ' | '.join('%s %s/%s' % (l.replace('_', ''), '-' if c is None else '%.1f' % c, '-' if h is None else '%.2f' % h) for l, c, h in v)))
