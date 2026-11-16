"""Which isomers of the 2023 paper are the candidates printed in the STX150 SI? Rows of the two SIs are aligned
through the experimental values (the STX tables list methyl protons one by one, the 2023 tables once), then
(1) every STX candidate is fitted against every 2023 candidate (RMS of the residuals of a straight line), and
(2) the difference between two STX candidates is correlated with the difference between two 2023 candidates,
    which removes the effect of functional and basis set.
Usage: stx_jnp_identify.py <jnp_tables.json> <stx out dir> <out csv> [stx level] [2023 level]"""
import sys, os, csv, json, math, collections, itertools
tabs = json.load(open(sys.argv[1])); stx_out = sys.argv[2]; out_csv = sys.argv[3]
SL = sys.argv[4] if len(sys.argv) > 4 else 'MMFF_PBE'; JL = sys.argv[5] if len(sys.argv) > 5 else 'MPW_B_PCM'
blocks = collections.OrderedDict()
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_val_blocks.csv'))): blocks.setdefault(r['block'], []).append(r)
cases = list(csv.DictReader(open(os.path.join(stx_out, 'stx_val_cases.csv'))))
THREE = ['B3L_A_GAS', 'MPW_B_PCM', 'WB9_C_SMD']
def structure(t):
    ncol = len(t['rows'][0]) - 2; lev = t['levels']
    if ncol == 13 and len(set(lev)) >= 10: return [0], [dict((l, 1 + i) for i, l in enumerate(lev[:12]))]
    if ncol == 13: return [0], [dict((l, 1 + 3 * k + i) for i, l in enumerate(THREE)) for k in range(4)]
    if ncol == 7: return [0], [dict((l, 1 + 3 * k + i) for i, l in enumerate(THREE)) for k in range(2)]
    if ncol == 8: return [0, 4], [dict((l, 1 + i) for i, l in enumerate(THREE)), dict((l, 5 + i) for i, l in enumerate(THREE))]
    if ncol == 14: return [0, 4], [dict((l, 1 + i) for i, l in enumerate(THREE)), dict((l, 5 + i) for i, l in enumerate(THREE))] + [dict((l, 8 + 3 * k + i) for i, l in enumerate(THREE)) for k in range(2)]
    return [], []
def groups(seq):
    out = []
    for i, k in enumerate(seq):
        if out and out[-1][0] == k: out[-1][1].append(i)
        else: out.append((k, [i]))
    return out
def align(jrows, e, srows):
    """-> [(row of the 2023 table, [rows of the STX table])]; rows paired through nucleus and experimental value"""
    gj = groups([(r[1], round(r[2 + e], 2)) for r in jrows]); gs = groups([(r['nucleus'], round(float(r['exp']), 2) if r['exp'] != '' else None) for r in srows])
    pairs = []; i = 0; lost = 0
    for k, idx in gj:
        j = i
        while j < len(gs) and gs[j][0] != k: j += 1
        if j == len(gs): lost += len(idx); continue
        sidx = gs[j][1]; i = j + 1
        if len(sidx) == len(idx): pairs += [(a, [b]) for a, b in zip(idx, sidx)]
        elif len(sidx) % len(idx) == 0:
            m = len(sidx) // len(idx); pairs += [(a, sidx[n * m:(n + 1) * m]) for n, a in enumerate(idx)]
        else: lost += len(idx)
    return pairs, lost
def fit_rms(x, y):
    n = len(x)
    if n < 4: return None
    mx, my = sum(x) / n, sum(y) / n; sxx = sum((a - mx) ** 2 for a in x); sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    if sxx == 0: return None
    s = sxy / sxx; i = my - s * mx
    return math.sqrt(sum((b - (s * a + i)) ** 2 for a, b in zip(x, y)) / n)
def corr(x, y):
    n = len(x)
    if n < 4: return None
    mx, my = sum(x) / n, sum(y) / n; sxx = sum((a - mx) ** 2 for a in x); syy = sum((b - my) ** 2 for b in y)
    if sxx < 1e-12 or syy < 1e-12: return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)
def key(vals): return sorted(round(v, 1) for v in vals)
def overlap(k, ks):
    pool = collections.Counter(k); same = 0
    for a in ks:
        if pool.get(a, 0) > 0: pool[a] -= 1; same += 1
    return same
sets = []
for ti, t in enumerate(tabs):
    ex, cand = structure(t)
    if len(cand) < 2: continue                 # section 9 tables only (all candidates in one table)
    for e in ex: sets.append((ti, e, key([r[2 + e] for r in t['rows'] if r[1] == 'C'])))
out = []; done = set()
for c in cases:
    name = c['case']; rows = blocks[name]
    ks = key([float(r['exp']) for r in rows if r['nucleus'] == 'C' and r['exp'] != ''])
    hits = [(ti, e) for ti, e, k in sets if overlap(k, ks) >= 0.85 * max(len(k), len(ks))]
    if not hits: continue
    ti, e = hits[0]; t = tabs[ti]; ex, cands = structure(t)
    pairs, lost = align(t['rows'], e, rows)
    cand_names = c['candidates'].split()
    if any(len(blocks[n]) != len(rows) for n in cand_names): print(name, 'candidate blocks differ in length'); continue
    def svec(block, nuc):
        v = []
        for a, bs in pairs:
            if t['rows'][a][1] != nuc: continue
            vals = [blocks[block][b][SL] for b in bs]
            v.append(None if any(x == '' for x in vals) else sum(float(x) for x in vals) / len(vals))
        return v
    def jvec(k, nuc): return [t['rows'][a][2 + cands[k][JL]] for a, bs in pairs if t['rows'][a][1] == nuc]
    rec = dict(case=name, stx_candidates=' '.join(cand_names), jnp_line=t['line'], jnp_candidates=len(cands), rows_2023=len(t['rows']), rows_paired=len(pairs), rows_lost=lost)
    # (1) fits
    for s in cand_names:
        best = {}
        for nuc in 'CH':
            for k in range(len(cands)):
                x = jvec(k, nuc); y = svec(s, nuc); xy = [(a, b) for a, b in zip(x, y) if b is not None]
                v = fit_rms([a for a, _ in xy], [b for _, b in xy]); best[(nuc, k)] = v
        rec['fit_' + ('correct' if s == name else s)] = ' '.join('%s%d=%s' % (nuc, k + 1, ('%.2f' % v) if v is not None else '-') for (nuc, k), v in best.items())
    # (2) differences, for the correct candidate against each alternative
    for s in cand_names:
        if s == name: continue
        res = []
        for nuc in 'CH':
            ds = [None if (a is None or b is None) else a - b for a, b in zip(svec(name, nuc), svec(s, nuc))]
            for k1, k2 in itertools.permutations(range(len(cands)), 2):
                dj = [a - b for a, b in zip(jvec(k1, nuc), jvec(k2, nuc))]
                xy = [(a, b) for a, b in zip(ds, dj) if a is not None]
                r = corr([a for a, _ in xy], [b for _, b in xy])
                if r is not None: res.append((r, nuc, k1 + 1, k2 + 1))
        for nuc in 'CH':
            rr = sorted([x for x in res if x[1] == nuc], reverse=True)
            if rr:
                rec['diff_%s_%s' % (nuc, s if s != name else 'x')] = 'best %d-%d r=%.2f; next %d-%d r=%.2f' % (rr[0][2], rr[0][3], rr[0][0], rr[1][2], rr[1][3], rr[1][0]) if len(rr) > 1 else 'best %d-%d r=%.2f' % (rr[0][2], rr[0][3], rr[0][0])
    out.append(rec)
keys = []
for r in out:
    for k in r:
        if k not in keys: keys.append(k)
w = csv.DictWriter(open(out_csv, 'w', newline=''), fieldnames=keys); w.writeheader(); w.writerows(out)
print('cases of the STX150 set found in section 9 of the 2023 SI: %d | STX level %s against 2023 level %s' % (len(out), SL, JL))
for r in out:
    d = [(k[5:], v) for k, v in r.items() if k.startswith('diff_')]
    print('%-5s 2023 candidates %d | rows paired %d of %d | %s' % (r['case'], r['jnp_candidates'], r['rows_paired'], r['rows_2023'], ' | '.join('%s: %s' % (k, v) for k, v in d)))
