"""2025 SI: the difference between two printed isomers at two levels must be nearly the same vector.
Rows with a cell that cannot be the shielding of its nucleus are left out, so that what is tested is the order of
the isomers and not the misprints. Usage: obc_internal.py <obc_compounds.json>"""
import sys, json, math, itertools
C = json.load(open(sys.argv[1]))
def corr(x, y):
    n = len(x)
    if n < 4: return None
    mx, my = sum(x) / n, sum(y) / n; sxx = sum((a - mx) ** 2 for a in x); syy = sum((b - my) ** 2 for b in y)
    if sxx < 1e-12 or syy < 1e-12: return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)
fmt = lambda rr: '/'.join('-' if v is None else '%.2f' % v for v in rr)
n = 0; flagged = []
print('difference between isomers a and b at two levels, r(13C)/r(1H); compounds with a correlation below -0.5:')
for num, c in C.items():
    bad = {(iso, l, i) for iso, l, i in c['bad_cells']}
    out = []; neg = False; strong = False
    for l1, l2 in itertools.combinations(c['levels'], 2):
        rr = []
        for nuc in 'CH':
            idx = [i for i, r in enumerate(c['rows']) if r['nucleus'] == nuc and not any((iso, l, i) in bad for iso in 'ab' for l in (l1, l2))]
            rr.append(corr([c['isomers']['a'][l1][i] - c['isomers']['b'][l1][i] for i in idx], [c['isomers']['a'][l2][i] - c['isomers']['b'][l2][i] for i in idx]))
        out.append('%s~%s %s' % (l1, l2, fmt(rr)))
        if any(v is not None and v < -0.5 for v in rr): neg = True
        if all(v is not None and v < -0.8 for v in rr): strong = True
    n += 1
    if neg: flagged.append((num, strong)); print('  compound %2s%s: %s' % (num, ' (pair with %s)' % c['pair_with'] if c['pair_with'] else '', ' | '.join(out)))
print('compounds tested: %d | with a correlation below -0.5: %d %s | below -0.8 for 13C and 1H at some pair of levels: %s' % (n, len(flagged), [x for x, _ in flagged], [x for x, s in flagged if s]))
