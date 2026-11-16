"""Does another error model order the candidates better? Per-nucleus terms other than the Gaussian tail of DP4.

From the scaled-shift errors of every assigned nucleus (table <prefix>_shifts.csv) the log-likelihood of each candidate is
rebuilt with a different term per nucleus, z = error / sigma (2.269 ppm for 13C, 0.187 ppm for 1H):
  gauss tail (DP4)   log(2 * Phi(-|z|))
  t tail, nu         log(2 * (1 - T_nu(|z|)))      heavier tails: one bad nucleus counts for less
  laplace            -|z|                           the candidates are ordered by a weighted mean absolute error
  capped gauss, c    -min(z, c)^2 / 2               errors beyond c sigma all count the same
For each model: how often the actual structure is first, the temperature that would make its probabilities honest,
and the size of the sets at matched coverage and under split conformal prediction, with two forms of the score.
Usage: error_models.py <shifts.csv> <truth.csv> <out.csv> [splits] [seed]       (needs numpy and scipy)
"""
import bisect
import collections
import csv
import math
import os
import random
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scores import fitted_temperature, softmax

SIGMA = {'C': 2.269372270818724, 'H': 0.18731058105269952}


def models():
    m = collections.OrderedDict()
    m['gauss tail (DP4)'] = lambda z: np.log(np.maximum(2 * stats.norm.sf(z), 1e-300))
    for nu in (10, 5, 3, 2):
        m['t tail, nu=%d' % nu] = (lambda nu: lambda z: np.log(np.maximum(2 * stats.t.sf(z, nu), 1e-300)))(nu)
    m['laplace'] = lambda z: -z
    for c in (3, 2):
        m['capped gauss, c=%d' % c] = (lambda c: lambda z: -np.minimum(z, c) ** 2 / 2)(c)
    return m


def evaluate(recs, target, splits, seed):
    """recs: dicts with wf, family, k, L (list), N -> metrics for two score forms"""
    out = {}
    fam = sorted({r['family'] for r in recs})
    by = collections.defaultdict(list)
    for r in recs:
        by[r['wf']].append(r)
    for form in ('probability', 'log ratio / sqrt(N)'):
        for r in recs:
            if form == 'probability':
                s = [1 - p for p in softmax(r['L'])]
            else:
                m = max(r['L'])
                s = [(m - x) / math.sqrt(r['N']) for x in r['L']]
            r['t'] = s[r['k']]
            r['s'] = sorted(s)
        sz = sg = 0
        for wf, rr in by.items():
            ts = sorted(r['t'] for r in rr)
            q = ts[int(math.ceil(target * len(ts))) - 1]
            for r in rr:
                n = bisect.bisect_right(r['s'], q + 1e-12)
                sz += n
                sg += n == 1
        rng = random.Random(seed)
        a = [0, 0, 0, 0]
        for _ in range(splits):
            f = fam[:]
            rng.shuffle(f)
            cal = set(f[:len(f) // 2])
            for wf, rr in by.items():
                ts = sorted(r['t'] for r in rr if r['family'] in cal)
                idx = int(math.ceil((len(ts) + 1) * target))
                q = ts[idx - 1] if 0 < idx <= len(ts) else float('inf')
                for r in rr:
                    if r['family'] in cal:
                        continue
                    n = bisect.bisect_right(r['s'], q + 1e-12)
                    a[0] += r['t'] <= q + 1e-12
                    a[1] += n
                    a[2] += n == 1
                    a[3] += 1
        out[form] = (sz / len(recs), sg / len(recs), a[0] / a[3], a[1] / a[3], a[2] / a[3])
    return out


def main(shifts, truth, out, splits=500, seed=20261114):
    splits, seed = int(splits), int(seed)
    T = {r['experiment']: int(r['actual_isomer']) for r in csv.DictReader(open(truth))}
    err = collections.defaultdict(lambda: {'C': [], 'H': []})
    fam = {}
    for r in csv.DictReader(open(shifts)):
        if r['experiment'] in T:
            err[(r['experiment'], r['level'], int(r['isomer']))][r['nucleus']].append(abs(float(r['shift_scaled']) - float(r['exp'])) / SIGMA[r['nucleus']])
            fam[r['experiment']] = r['archive']
    keys = sorted({k[:2] for k in err})
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['error_model', 'actual_first', 'records', 'fitted_temperature', 'target', 'score', 'matched_size', 'matched_single',
                'conformal_coverage', 'conformal_size', 'conformal_single'])
    print('records %d' % len(keys))
    print('%-20s %12s %6s | 90%% target, probability: matched size, conformal coverage and size | per-nucleus log ratio: the same | 80%% target: sizes (probability, log ratio)'
          % ('error model', 'actual first', 'T'))
    for name, g in models().items():
        recs = []
        for exp, level in keys:
            n = max(k[2] for k in err if k[:2] == (exp, level)) if False else None
            isos = sorted(k[2] for k in err_index[(exp, level)])
            LC = [float(np.sum(g(np.array(err[(exp, level, i)]['C'])))) for i in isos]
            LH = [float(np.sum(g(np.array(err[(exp, level, i)]['H'])))) for i in isos]
            first = err[(exp, level, isos[0])]
            recs.append(dict(wf=level, family=fam[exp], k=T[exp] - 1, LC=LC, LH=LH, L=[a + b for a, b in zip(LC, LH)],
                             N=len(first['C']) + len(first['H'])))
        top = sum(1 for r in recs if max(range(len(r['L'])), key=lambda i: r['L'][i]) == r['k'])
        temp = fitted_temperature(recs)
        e90 = evaluate(recs, 0.90, splits, seed)
        e80 = evaluate(recs, 0.80, splits, seed)
        for target, e in ((0.90, e90), (0.80, e80)):
            for form, v in e.items():
                w.writerow([name, top, len(recs), '%.2f' % temp, target, form] + ['%.4f' % x for x in v])
        p, l = e90['probability'], e90['log ratio / sqrt(N)']
        print('%-20s %4d (%4.1f%%) %6.2f | %5.2f  %5.1f%%  %5.2f | %5.2f  %5.1f%%  %5.2f | %5.2f  %5.2f'
              % (name, top, 100.0 * top / len(recs), temp, p[0], 100 * p[2], p[3], l[0], 100 * l[2], l[3],
                 e80['probability'][3], e80['log ratio / sqrt(N)'][3]))


if __name__ == '__main__':
    # index of the isomers of each record, built once
    err_index = collections.defaultdict(set)
    _T = {r['experiment'] for r in csv.DictReader(open(sys.argv[2]))}
    for _r in csv.DictReader(open(sys.argv[1])):
        if _r['experiment'] in _T:
            err_index[(_r['experiment'], _r['level'])].add((_r['experiment'], _r['level'], int(_r['isomer'])))
    main(*sys.argv[1:6])
