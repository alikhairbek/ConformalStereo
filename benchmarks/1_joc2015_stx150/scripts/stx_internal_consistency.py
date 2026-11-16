"""Internal test of the shielding tables printed in the STX150 SI (validation set). For two candidates of one set,
the difference between their shieldings is a property of the two structures: at two levels that share the geometries
(xTB_PBE, xTB_mPW, xTB_wB9; B3L_PBE, B3L_wB9) the two difference vectors must be almost the same. A correlation near -1
means the two candidates are exchanged at one of the levels; a correlation near 0 means that a column belongs to
something else. Usage: stx_internal_consistency.py <stx out dir> <out csv>"""
import sys, os, csv, math, collections, itertools
stx_out, out_csv = sys.argv[1:3]
LEVELS = ['MMFF_PBE', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9', 'B3L_PBE', 'B3L_wB9']
SAME_GEOM = [('xTB_PBE', 'xTB_mPW'), ('xTB_PBE', 'xTB_wB9'), ('xTB_mPW', 'xTB_wB9'), ('B3L_PBE', 'B3L_wB9')]
blocks = collections.OrderedDict()
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_val_blocks.csv'))): blocks.setdefault(r['block'], []).append(r)
cases = list(csv.DictReader(open(os.path.join(stx_out, 'stx_val_cases.csv'))))
def corr(x, y):
    n = len(x)
    if n < 5: return None
    mx, my = sum(x) / n, sum(y) / n; sxx = sum((a - mx) ** 2 for a in x); syy = sum((b - my) ** 2 for b in y)
    if sxx < 1e-9 or syy < 1e-9: return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)
def diff(a, b, lev, nuc):
    out = []
    for ra, rb in zip(blocks[a], blocks[b]):
        if ra['nucleus'] != nuc: continue
        out.append(None if ra[lev] == '' or rb[lev] == '' else float(ra[lev]) - float(rb[lev]))
    return out
rows = []; seen = set()
for c in cases:
    cand = c['candidates'].split()
    for a, b in itertools.combinations(cand, 2):
        if (a, b) in seen or len(blocks[a]) != len(blocks[b]): continue
        seen.add((a, b))
        for l1, l2 in itertools.combinations(LEVELS, 2):
            rec = dict(pair=a + ' / ' + b, level_1=l1, level_2=l2, same_geometries=int((l1, l2) in SAME_GEOM))
            for nuc in 'CH':
                d1, d2 = diff(a, b, l1, nuc), diff(a, b, l2, nuc)
                xy = [(x, y) for x, y in zip(d1, d2) if x is not None and y is not None]
                r = corr([x for x, _ in xy], [y for _, y in xy]); rec['r_' + nuc] = '' if r is None else round(r, 3); rec['n_' + nuc] = len(xy)
                rec['size_%s_1' % nuc] = round(math.sqrt(sum(x * x for x, _ in xy) / len(xy)), 3) if xy else ''
                rec['size_%s_2' % nuc] = round(math.sqrt(sum(y * y for _, y in xy) / len(xy)), 3) if xy else ''
            rows.append(rec)
w = csv.DictWriter(open(out_csv, 'w', newline=''), fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
sg = [r for r in rows if r['same_geometries']]
def cls(r):
    v = [r[k] for k in ('r_C', 'r_H') if r[k] != '']
    if not v: return 'no test'
    if all(x > 0.8 for x in v): return 'consistent (r > 0.8, 13C and 1H)'
    if all(x < -0.8 for x in v): return 'exchanged (r < -0.8, 13C and 1H)'
    if any(x < -0.5 for x in v): return 'partly exchanged (r < -0.5 for one nucleus)'
    if all(x > 0.5 for x in v): return 'weaker (r 0.5 to 0.8)'
    return 'unrelated (r below 0.5)'
cnt = collections.Counter(cls(r) for r in sg)
print('pairs of candidates: %d | tests at two levels with the same geometries: %d' % (len(seen), len(sg)))
for k, v in cnt.most_common(): print('   %-48s %d' % (k, v))
print('tests that are not consistent:')
for r in sg:
    k = cls(r)
    if not k.startswith('consistent') and k != 'no test':
        print('   %-18s %s / %s  r(13C) %s  r(1H) %s   [%s]' % (r['pair'], r['level_1'], r['level_2'], r['r_C'], r['r_H'], k.split(' (')[0]))
