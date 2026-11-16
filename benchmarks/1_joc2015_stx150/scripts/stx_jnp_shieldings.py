"""Compare the shieldings printed in the STX150 SI (validation set) with those printed in the SI of the 2023 paper
(MMFF geometries) for the same experimental data sets. For every STX candidate block and every candidate of the
2023 table: RMS of the residuals of a straight-line fit, 13C and 1H, rows paired in table order.
Usage: stx_jnp_shieldings.py <jnp_tables.json> <stx out dir> <out csv> [case ...]"""
import sys, os, csv, json, math, collections
tabs = json.load(open(sys.argv[1])); stx_out = sys.argv[2]; out_csv = sys.argv[3]; focus = sys.argv[4:]
LEVELS = ['MMFF_PBE', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9', 'B3L_PBE', 'B3L_wB9']
blocks = collections.OrderedDict()
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_val_blocks.csv'))): blocks.setdefault(r['block'], []).append(r)
cases = list(csv.DictReader(open(os.path.join(stx_out, 'stx_val_cases.csv'))))

def structure(t):
    """-> list of experimental columns, list of candidates [{level: column}] (columns counted after label and nucleus)"""
    ncol = len(t['rows'][0]) - 2; lev = t['levels']
    three = ['B3L_A_GAS', 'MPW_B_PCM', 'WB9_C_SMD']
    if ncol == 13 and len(set(lev)) >= 10: return [0], [dict((l, 1 + i) for i, l in enumerate(lev[:12]))]
    if ncol == 13: return [0], [dict((l, 1 + 3 * k + i) for i, l in enumerate(three)) for k in range(4)]
    if ncol == 7: return [0], [dict((l, 1 + 3 * k + i) for i, l in enumerate(three)) for k in range(2)]
    if ncol == 8: return [0, 4], [dict((l, 1 + i) for i, l in enumerate(three)), dict((l, 5 + i) for i, l in enumerate(three))]
    if ncol == 14: return [0, 4], [dict((l, 1 + i) for i, l in enumerate(three)), dict((l, 5 + i) for i, l in enumerate(three))] + [dict((l, 8 + 3 * k + i) for i, l in enumerate(three)) for k in range(2)]
    if ncol == 12: return [], [dict((l, 3 * k + i) for i, l in enumerate(three)) for k in range(4)]
    return [], []

def fit_rms(x, y):
    n = len(x)
    if n < 4: return None
    mx, my = sum(x) / n, sum(y) / n; sxx = sum((a - mx) ** 2 for a in x); sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    if sxx == 0: return None
    s = sxy / sxx; i = my - s * mx
    return math.sqrt(sum((b - (s * a + i)) ** 2 for a, b in zip(x, y)) / n)

def key(vals): return sorted(round(v, 1) for v in vals)
def overlap(k, ks):
    pool = collections.Counter(k); same = 0
    for a in ks:
        if pool.get(a, 0) > 0: pool[a] -= 1; same += 1
    return same
# experimental data sets of the 2023 SI
sets = []     # (table index, exp column, sorted 13C key)
for ti, t in enumerate(tabs):
    ex, cand = structure(t)
    for e in ex: sets.append((ti, e, key([r[2 + e] for r in t['rows'] if r[1] == 'C'])))
res = []; summary = collections.Counter()
for c in cases:
    name = c['case']; rows = blocks[name]
    ks = key([float(r['exp']) for r in rows if r['nucleus'] == 'C' and r['exp'] != ''])
    hits = [(ti, e) for ti, e, k in sets if overlap(k, ks) >= 0.85 * max(len(k), len(ks))]
    if not hits: summary['no table'] += 1; continue
    summary['matched'] += 1
    for cand_name in c['candidates'].split():
        b = blocks[cand_name]
        for ti, e in hits:
            t = tabs[ti]; ex, cands = structure(t)
            if len(t['rows']) != len(b): summary['row count differs'] += 1; continue
            if any(r[1] != q['nucleus'] for r, q in zip(t['rows'], b)): summary['nucleus order differs'] += 1; continue
            for k, cd in enumerate(cands):
                for jl, col in cd.items():
                    if jl not in ('MPW_B_PCM', 'B3L_A_GAS', 'WB9_C_SMD'): continue
                    for sl in LEVELS:
                        for nuc in 'CH':
                            x = []; y = []
                            for r, q in zip(t['rows'], b):
                                if r[1] == nuc and q[sl] != '': x.append(r[2 + col]); y.append(float(q[sl]))
                            v = fit_rms(x, y)
                            if v is not None: res.append(dict(case=name, stx_block=cand_name, stx_level=sl, jnp_line=t['line'], jnp_exp_col=e, jnp_candidate=k + 1, jnp_level=jl, nucleus=nuc, n=len(x), rms=round(v, 3)))
w = csv.DictWriter(open(out_csv, 'w', newline=''), fieldnames=list(res[0].keys())); w.writeheader(); w.writerows(res)
print(dict(summary), '| rows written', len(res))
def show(name):
    sel = [r for r in res if r['case'] == name and r['stx_level'] == 'MMFF_PBE' and r['jnp_level'] == 'MPW_B_PCM']
    for blk in dict.fromkeys(r['stx_block'] for r in sel):
        parts = []
        for r in sel:
            if r['stx_block'] == blk: parts.append('%s#%d %s=%.2f' % (r['jnp_line'], r['jnp_candidate'], r['nucleus'], r['rms']))
        print('  %-6s exp of %-5s STX block %-9s MMFF_PBE vs 2023 MPW_B_PCM: %s' % ('', name, blk, ' | '.join(parts)))
for name in focus: show(name)
