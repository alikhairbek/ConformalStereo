"""Exploration, reported for transparency: would a score that also divides by the number of candidates even out the
coverage further? Not adopted: there is no principled exponent, and the sets of molecules with many candidates grow
strongly. Usage: explore_candidates.py <pooled_candidates.csv>

Scores: log-likelihood below the best candidate divided by N**a * log2(K)**c, for a in (0, 1/2, 1) and c in (0, 1).
One family out, cut-off from the own workflow; coverage and mean set size by number of candidates and of shifts, and
mean set size per benchmark."""
import sys, os, math, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cslib as cs

recs = cs.load(sys.argv[1])
by_wf = collections.OrderedDict()
for r in recs:
    by_wf.setdefault(r['wf'], []).append(r)
SC = []
for a, an in ((0.0, ''), (0.5, ' / sqrt(N)'), (1.0, ' / N')):
    for c, cn in ((0, ''), (1, ' / log2(K)')):
        name = 'ratio' + an + cn
        SC.append(name)
        for r in recs:
            v = (r['ell'].max() - r['ell']) / r['N'] ** a / (math.log2(r['K']) if c else 1.0)
            r['t'][name], r['s'][name], r['tops'][name] = float(v[r['k']]), np.sort(v), float(v[r['top']])
nclass = lambda r: '<25' if r['N'] < 25 else '25-39' if r['N'] < 40 else '40+'
for target in (0.90, 0.80):
    print('target %d%%: coverage and mean set size' % (100 * target))
    for s in SC:
        thr = {}
        for wf, rr in by_wf.items():
            for f in {x['family'] for x in rr}:
                thr[(wf, f)] = cs.threshold([(x['t'][s], x['family']) for x in rr if x['family'] != f], target)
        cov, size = [], []
        for r in recs:
            c, z = cs.apply_set(r, s, thr[(r['wf'], r['family'])])
            cov.append(c)
            size.append(z)
        cells = []
        for title, fn, classes in [('K', lambda r: cs.kclass(r['K']), cs.KCLASSES), ('N', nclass, ['<25', '25-39', '40+'])]:
            for c in classes:
                idx = [i for i, r in enumerate(recs) if fn(r) == c]
                x = cs.summarize([recs[i] for i in idx], [cov[i] for i in idx], [size[i] for i in idx])
                cells.append('%s %s: %4.1f%% %5.2f' % (title, c, 100 * x['coverage'], x['size']))
        tot = []
        for b in ['JOC2015', 'STX150', 'DP4-AI', 'Cam2017', 'Cam2019']:
            idx = [i for i, r in enumerate(recs) if r['bench'] == b]
            tot.append('%.2f' % cs.summarize([recs[i] for i in idx], [cov[i] for i in idx], [size[i] for i in idx])['size'])
        print('  %-28s | %s | per benchmark %s' % (s, ' | '.join(cells), ' '.join(tot)))
