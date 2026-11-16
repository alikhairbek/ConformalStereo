"""2023 SI, section 9 (compounds 101-121, three MMFF levels): which printed candidate fits the experimental column?
Each table has one or two experimental columns. For every candidate and level: RMS error after linear scaling,
13C and 1H. Next to it the published probability of the correct isomer for the STX case with these data.
The point: in which column does the correct isomer sit?
Usage: jnp_table_order.py <jnp_tables.json> <jnp_probabilities.csv> <stx out dir> <stx_jnp_identify csv (MMFF_PBE)>"""
import sys, os, csv, json, math, collections
T = json.load(open(sys.argv[1])); stx_out = sys.argv[3]
prob = {r['case']: r for r in csv.DictReader(open(sys.argv[2]))}
ident = {r['case']: r for r in csv.DictReader(open(sys.argv[4]))}
blocks = collections.OrderedDict()
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_val_blocks.csv'))): blocks.setdefault(r['block'], []).append(r)
THREE = ['B3L_A_GAS', 'MPW_B_PCM', 'WB9_C_SMD']
def structure(t):
    ncol = len(t['rows'][0]) - 2; lev = t['levels']
    if ncol == 13 and len(set(lev)) >= 10: return None
    if ncol == 13: return [0], [dict((l, 1 + 3 * k + i) for i, l in enumerate(THREE)) for k in range(4)]
    if ncol == 7: return [0], [dict((l, 1 + 3 * k + i) for i, l in enumerate(THREE)) for k in range(2)]
    if ncol == 8: return [0, 4], [dict((l, 1 + i) for i, l in enumerate(THREE)), dict((l, 5 + i) for i, l in enumerate(THREE))]
    if ncol == 14: return [0, 4], [dict((l, 1 + i) for i, l in enumerate(THREE)), dict((l, 5 + i) for i, l in enumerate(THREE))] + [dict((l, 8 + 3 * k + i) for i, l in enumerate(THREE)) for k in range(2)]
def scaled_rms(exp, sig):
    n = len(exp)
    if n < 3: return float('nan')
    mx, my = sum(sig) / n, sum(exp) / n; sxx = sum((s - mx) ** 2 for s in sig); sxy = sum((s - mx) * (e - my) for s, e in zip(sig, exp))
    if sxx == 0: return float('nan')
    b = sxy / sxx; a = my - b * mx
    return math.sqrt(sum((e - (a + b * s)) ** 2 for s, e in zip(sig, exp)) / n)
def key(v): return sorted(round(x, 1) for x in v)
stx_exp = {n: key([float(r['exp']) for r in rows if r['nucleus'] == 'C' and r['exp'] != '']) for n, rows in blocks.items() if any(r['exp'] != '' for r in rows)}
def overlap(k, ks):
    pool = collections.Counter(k); s = 0
    for a in ks:
        if pool.get(a, 0) > 0: pool[a] -= 1; s += 1
    return s
pcols = [c for c in next(iter(prob.values())) if 'ScTMS+Uns' in c and 'ScMSTD' not in c]
print('case   exp col | published P(correct) %s | best-fitting candidate by level (13C+1H rank sum), errors of all candidates at MPW_B_PCM (13C/1H) | STX order (13C,1H)' % [c.split('|')[0] for c in pcols])
for t in T:
    st = structure(t)
    if not st or t['line'] < 26498: continue
    ex, cands = st
    for e in ex:
        k = key([r[2 + e] for r in t['rows'] if r[1] == 'C'])
        hit = [n for n, ks in stx_exp.items() if overlap(k, ks) >= 0.85 * max(len(k), len(ks))]
        name = hit[0] if hit else '?'
        best = []
        for lev in THREE:
            score = []
            for c in cands:
                rc = scaled_rms([r[2 + e] for r in t['rows'] if r[1] == 'C'], [r[2 + c[lev]] for r in t['rows'] if r[1] == 'C'])
                rh = scaled_rms([r[2 + e] for r in t['rows'] if r[1] == 'H'], [r[2 + c[lev]] for r in t['rows'] if r[1] == 'H'])
                score.append((rc, rh))
            # likelihood-like score: sum of squared standardised errors with typical sigmas 2.0 and 0.15 ppm
            nC = sum(r[1] == 'C' for r in t['rows']); nH = sum(r[1] == 'H' for r in t['rows'])
            chi = [nC * (rc / 2.0) ** 2 + nH * (rh / 0.15) ** 2 for rc, rh in score]
            best.append(1 + min(range(len(chi)), key=lambda i: chi[i]))
            if lev == 'MPW_B_PCM': shown = ' '.join('%d:%.2f/%.3f' % (i + 1, rc, rh) for i, (rc, rh) in enumerate(score))
        pp = [prob[name][c] if name in prob else '-' for c in pcols]
        idn = ident.get(name, {})
        order = '; '.join('%s' % v.split(';')[0] for kk, v in idn.items() if kk.startswith('diff_') and v)
        print('%-6s %d | %-22s | best %s | %s | %s' % (name, e, ' '.join(str(x) for x in pp), best, shown, order))
