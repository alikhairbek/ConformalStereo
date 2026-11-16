"""Candidate pairs of the STX150 validation set whose two members are exchanged between levels: the difference
vector [candidate 1 - candidate 2] at one level is the negative of the one at another level (r < -0.8 for 13C and 1H).
Reads the table written by stx_internal_consistency.py. Usage: stx_exchanged_candidates.py <internal consistency csv>"""
import sys, csv, collections
rows = list(csv.DictReader(open(sys.argv[1])))
f = lambda v: None if v == '' else float(v)
LEV = ['MMFF_PBE', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9', 'B3L_PBE', 'B3L_wB9']
by = collections.defaultdict(list)
for r in rows: by[r['pair']].append(r)
print('candidate pairs tested: %d | level pairs per candidate pair: 15' % len(by))
neg = collections.defaultdict(list); one = collections.defaultdict(list)
for pair, rs in by.items():
    for r in rs:
        c, h = f(r['r_C']), f(r['r_H'])
        if c is not None and h is not None and c < -0.8 and h < -0.8: neg[pair].append((r['level_1'], r['level_2'], c, h))
        elif (c is not None and c < -0.8) or (h is not None and h < -0.8): one[pair].append((r['level_1'], r['level_2'], c, h))
print('pairs with the two candidates exchanged between levels (r < -0.8 for 13C and for 1H):', len(neg))
for pair, v in neg.items():
    g = {LEV[0]: 0}
    for _ in range(3):
        for r in by[pair]:
            c = f(r['r_C'])
            if c is None or abs(c) <= 0.8: continue
            for a, b in ((r['level_1'], r['level_2']), (r['level_2'], r['level_1'])):
                if a in g and b not in g: g[b] = g[a] if c > 0 else 1 - g[a]
    print('  %-18s group 1: %s | group 2: %s | r(13C), r(1H): %s' % (pair, ' '.join(k for k in LEV if g.get(k) == 0), ' '.join(k for k in LEV if g.get(k) == 1), '; '.join('%s~%s %.2f, %.2f' % x for x in v)))
print('pairs with r < -0.8 for one nucleus only, at some level pair:', len(one))
for pair, v in one.items(): print('  %-18s %s' % (pair, '; '.join('%s~%s %s, %s' % x for x in v[:6])))
