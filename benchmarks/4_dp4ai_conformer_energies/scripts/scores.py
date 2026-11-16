"""Conformal scores for stereochemical assignment, compared on one benchmark.

Every score orders the candidates of a record (one experiment under one workflow) from most to least plausible; a
prediction set keeps the candidates whose score is at most a threshold. A better score gives smaller sets at the same
coverage. Two evaluations:
  1. matched coverage: per workflow the smallest threshold that covers the target share of the actual structures
     (threshold chosen on the data themselves: a property of the score, no calibration noise)
  2. split conformal: threshold from half of the families (molecule archives), sets on the other half
Scores, all from the DP4 log-likelihoods L of the candidates (13C + 1H) and the number of nuclei N:
  probability          1 - softmax(L)                       the usual choice
  tempered T           1 - softmax(L / T)                   fixed T, no fitting
  log ratio            max(L) - L                            does not depend on T
  log ratio / N, / sqrt(N)                                   the same per nucleus
  cumulative (T)       probability mass of the candidates at least as probable (adaptive prediction sets)
  rank                 position in the ordering
  carbon only, proton only   1 - softmax(L_C), 1 - softmax(L_H)
Usage: scores.py <loglik.csv> <truth.csv> <out.csv> [splits] [seed] [workflow filter: all | own]
"""
import bisect
import collections
import csv
import math
import random
import sys


def softmax(v):
    m = max(v)
    w = [math.exp(x - m) for x in v]
    z = sum(w)
    return [x / z for x in w]


def make_scores(rec):
    L = [a + b for a, b in zip(rec['LC'], rec['LH'])]
    n = rec['N']
    out = {}
    out['probability'] = [1 - p for p in softmax(L)]
    for T in (2, 3, 5, 10):
        out['tempered T=%d' % T] = [1 - p for p in softmax([x / T for x in L])]
    m = max(L)
    out['log ratio'] = [m - x for x in L]
    out['log ratio / N'] = [(m - x) / n for x in L]
    out['log ratio / sqrt(N)'] = [(m - x) / math.sqrt(n) for x in L]
    for T in (1, 3):
        p = softmax([x / T for x in L])
        out['cumulative T=%d' % T] = [sum(q for q in p if q >= pk) for pk in p]
    order = sorted(range(len(L)), key=lambda i: -L[i])
    rank = [0] * len(L)
    for pos, i in enumerate(order):
        rank[i] = pos + 1
    out['rank'] = rank
    out['carbon only'] = [1 - p for p in softmax(rec['LC'])]
    out['proton only'] = [1 - p for p in softmax(rec['LH'])]
    return out


def fitted_temperature(recs):
    """T that maximises the likelihood of the actual structures (one number for all records given)"""
    def nll(T):
        s = 0.0
        for r in recs:
            L = [(a + b) / T for a, b in zip(r['LC'], r['LH'])]
            m = max(L)
            s -= L[r['k']] - m - math.log(sum(math.exp(x - m) for x in L))
        return s
    lo, hi = 0.2, 60.0
    for _ in range(60):
        a, b = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if nll(a) < nll(b):
            hi = b
        else:
            lo = a
    return (lo + hi) / 2


def main(loglik, truth, out, splits=1000, seed=20261114, which='all'):
    splits, seed = int(splits), int(seed)
    T = {r['experiment']: int(r['actual_isomer']) for r in csv.DictReader(open(truth))}
    raw = collections.defaultdict(list)
    for r in csv.DictReader(open(loglik)):
        if r['experiment'] in T and (which == 'all' or r['shielding_level'] == r['energy_level']):
            raw[(r['experiment'], r['shielding_level'], r['energy_level'])].append(r)
    recs = []
    for (exp, sl, el), rows in sorted(raw.items()):
        rows.sort(key=lambda r: int(r['isomer']))
        rec = dict(exp=exp, wf=(sl, el), family=rows[0]['archive'], geo=rows[0]['geometry'], k=T[exp] - 1,
                   LC=[float(r['logL_C']) for r in rows], LH=[float(r['logL_H']) for r in rows],
                   N=int(rows[0]['n_C']) + int(rows[0]['n_H']), K=len(rows))
        sc = make_scores(rec)
        rec['truth'] = {name: v[rec['k']] for name, v in sc.items()}
        rec['sorted'] = {name: sorted(v) for name, v in sc.items()}
        rec['best'] = {name: min(range(len(v)), key=lambda i: v[i]) for name, v in sc.items()}
        recs.append(rec)
    names = list(recs[0]['truth'])
    by = collections.defaultdict(list)
    for r in recs:
        by[r['wf']].append(r)
    fam = sorted({r['family'] for r in recs})
    print('records %d, workflows %d, experiments %d, families %d; candidates per record %.2f on average'
          % (len(recs), len(by), len({r['exp'] for r in recs}), len(fam), sum(r['K'] for r in recs) / len(recs)))
    print('temperature that maximises the likelihood of the actual structures: %.2f (all records)' % fitted_temperature(recs))
    for g in sorted({r['geo'] for r in recs}):
        rr = [r for r in recs if r['geo'] == g]
        if len({r['geo'] for r in recs}) > 1:
            print('   %-20s %.2f (%d records)' % (g, fitted_temperature(rr), len(rr)))
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['evaluation', 'target', 'score', 'coverage', 'mean_set_size', 'single_structure_sets', 'single_and_right'])

    def size(r, name, q):
        return bisect.bisect_right(r['sorted'][name], q + 1e-12)

    for target in (0.90, 0.80):
        print()
        print('target %.0f%%  %-22s | matched coverage: size, single, single and right | split conformal: coverage, size, single, single and right'
              % (100 * target, 'score'))
        # ---- 1. matched coverage
        matched = {}
        for name in names:
            cov = sz = sg = sr = 0
            for wf, rr in by.items():
                ts = sorted(r['truth'][name] for r in rr)
                q = ts[int(math.ceil(target * len(ts))) - 1]
                for r in rr:
                    s = size(r, name, q)
                    hit = r['truth'][name] <= q + 1e-12
                    cov += hit
                    sz += s
                    sg += s == 1
                    sr += s == 1 and hit
            n = len(recs)
            matched[name] = (cov / n, sz / n, sg / n, sr / n)
            w.writerow(['matched coverage', target, name] + ['%.4f' % x for x in matched[name]])
        # ---- 2. split conformal by family
        rng = random.Random(seed)
        acc = {name: [0, 0, 0, 0, 0] for name in names}
        for _ in range(splits):
            f = fam[:]
            rng.shuffle(f)
            cal = set(f[:len(f) // 2])
            for wf, rr in by.items():
                c = [r for r in rr if r['family'] in cal]
                t = [r for r in rr if r['family'] not in cal]
                for name in names:
                    ts = sorted(r['truth'][name] for r in c)
                    idx = int(math.ceil((len(ts) + 1) * target))
                    q = ts[idx - 1] if 0 < idx <= len(ts) else float('inf')
                    a = acc[name]
                    for r in t:
                        s = size(r, name, q)
                        hit = r['truth'][name] <= q + 1e-12
                        a[0] += hit
                        a[1] += s
                        a[2] += s == 1
                        a[3] += s == 1 and hit
                        a[4] += 1
        for name in names:
            a = acc[name]
            sp = (a[0] / a[4], a[1] / a[4], a[2] / a[4], a[3] / a[4])
            w.writerow(['split conformal', target, name] + ['%.4f' % x for x in sp])
            m = matched[name]
            print('            %-22s | %5.2f  %5.1f%%  %5.1f%% | %5.1f%%  %5.2f  %5.1f%%  %5.1f%%'
                  % (name, m[1], 100 * m[2], 100 * m[3], 100 * sp[0], sp[1], 100 * sp[2], 100 * sp[3]))


if __name__ == '__main__':
    main(*sys.argv[1:7])
