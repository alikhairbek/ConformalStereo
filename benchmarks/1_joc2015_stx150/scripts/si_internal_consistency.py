"""Internal test of the shielding tables of the 2023 SI (section 9, three MMFF levels), as for the STX150 SI: the
difference between two printed candidates at two levels must be nearly the same vector.
(The same test for the 2025 SI is obc_internal.py.)  Usage: si_internal_consistency.py jnp <jnp_tables.json>"""
import sys, json, math, itertools
kind, path = sys.argv[1:3]; T = json.load(open(path))
def corr(x, y):
    n = len(x)
    if n < 4: return None
    mx, my = sum(x) / n, sum(y) / n; sxx = sum((a - mx) ** 2 for a in x); syy = sum((b - my) ** 2 for b in y)
    if sxx < 1e-12 or syy < 1e-12: return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)
fmt = lambda rr: '/'.join('-' if v is None else '%.2f' % v for v in rr)
THREE = ['B3L_A_GAS', 'MPW_B_PCM', 'WB9_C_SMD']
def structure(t):
    ncol = len(t['rows'][0]) - 2; lev = t['levels']
    if ncol == 13 and len(set(lev)) >= 10: return None
    if ncol == 13: return [dict((l, 1 + 3 * k + i) for i, l in enumerate(THREE)) for k in range(4)]
    if ncol == 7: return [dict((l, 1 + 3 * k + i) for i, l in enumerate(THREE)) for k in range(2)]
    if ncol == 8: return [dict((l, 1 + i) for i, l in enumerate(THREE)), dict((l, 5 + i) for i, l in enumerate(THREE))]
    if ncol == 14: return [dict((l, 1 + i) for i, l in enumerate(THREE)), dict((l, 5 + i) for i, l in enumerate(THREE))] + [dict((l, 8 + 3 * k + i) for i, l in enumerate(THREE)) for k in range(2)]
print('2023 SI, section 9: difference between two printed candidates at two of the three MMFF levels, r(13C)/r(1H)')
n = bad = 0
for t in T:
    c = structure(t)
    if not c or t['line'] < 26498: continue
    for i, j in itertools.combinations(range(len(c)), 2):
        out = []; flag = False
        for l1, l2 in itertools.combinations(THREE, 2):
            rr = []
            for nuc in 'CH':
                rows = [r for r in t['rows'] if r[1] == nuc]
                rr.append(corr([r[2 + c[i][l1]] - r[2 + c[j][l1]] for r in rows], [r[2 + c[i][l2]] - r[2 + c[j][l2]] for r in rows]))
            out.append('%s~%s %s' % (l1[:3], l2[:3], fmt(rr)))
            if any(v is not None and v < 0.5 for v in rr): flag = True
        n += 1; bad += flag
        if flag: print('  table at line %5d (first 13C value %s) candidates %d-%d: %s' % (t['line'], t['rows'][0][2], i + 1, j + 1, ' | '.join(out)))
print('candidate pairs tested: %d | with r below 0.5 for some pair of levels: %d' % (n, bad))
