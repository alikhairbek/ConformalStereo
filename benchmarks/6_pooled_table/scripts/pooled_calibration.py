"""One calibration for everything, or one per kind of calculation? Conformal calibration across the five benchmarks.

Records come from pool.py. Score: 1 - probability of the candidate. Families (molecules that share candidates, or that
occur in two benchmarks) are the exchangeable units: a calibration score is weighted by one over the number of records
of its family in the calibration pool, and coverage and set size are averaged per test family (hierarchical conformal
prediction). Designs, each with random halves of the families as calibration and test:
  own level      one threshold per benchmark and level (the reference; needs calibration data of that very level)
  by benchmark   one threshold per benchmark, all its levels together
  by group       one threshold per group = method and kind of geometry, across benchmarks and levels
  one for all    a single threshold for every record
and without any splitting:
  transfer       threshold from all other benchmarks (families of the test benchmark left out), applied to one benchmark;
                 once from every other record, once from the records of the same group only
Usage: pooled_calibration.py <pooled_candidates.csv> <out dir> [splits] [seed] [column: log_ratio | log_ratio_energy_4]
       [only records that have this column] [tag for the output file names]
"""
import bisect
import collections
import csv
import math
import os
import random
import sys


def softmax(v):
    m = max(v)
    w = [math.exp(x - m) for x in v]
    z = sum(w)
    return [x / z for x in w]


def group_of(method, geo):
    if method != 'DP4':
        return '%s, %s geometries' % (method, geo)
    if geo == 'MMFF':
        return 'DP4, MMFF geometries'
    if geo in ('PM6', 'PM7'):
        return 'DP4, semi-empirical geometries'
    return 'DP4, HF or DFT geometries'


def threshold(items, target):
    """items: [(score, family)] -> hierarchical conformal threshold"""
    if not items:
        return float('inf')
    cnt = collections.Counter(f for _, f in items)
    G = len(cnt)
    acc = 0.0
    for s, f in sorted(items):
        acc += 1.0 / cnt[f] / (G + 1)
        if acc >= target - 1e-12:
            return s
    return float('inf')


def fitted_temperature(recs):
    def nll(T):
        s = 0.0
        for r in recs:
            L = [x / T for x in r['L']]
            m = max(L)
            s -= L[r['k']] - m - math.log(sum(math.exp(x - m) for x in L))
        return s
    lo, hi = 0.05, 60.0
    for _ in range(70):
        a, b = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if nll(a) < nll(b):
            hi = b
        else:
            lo = a
    return (lo + hi) / 2


def main(table, out, splits=500, seed=20261114, column='log_ratio', require='', tag=''):
    splits, seed = int(splits), int(seed)
    raw = collections.OrderedDict()
    for r in csv.DictReader(open(table)):
        if r[column] == '' or (require and r[require] == ''):
            continue
        raw.setdefault((r['benchmark'], r['method'], r['geometry'], r['level'], r['case']), []).append(r)
    recs = []
    for (bench, method, geo, level, case), rr in raw.items():
        L = [float(r[column]) for r in rr]
        p = softmax(L)
        k = [int(r['is_true']) for r in rr].index(1)
        recs.append(dict(bench=bench, group=group_of(method, geo), level=(bench, level), family=rr[0]['family'], L=L, k=k, K=len(rr),
                         t=1 - p[k], s=sorted(1 - x for x in p), ptop=max(p), first=p[k] >= max(p) - 1e-15, pa=p[k]))
    benches = list(collections.OrderedDict((r['bench'], 1) for r in recs))
    groups = sorted({r['group'] for r in recs})
    fams = sorted({r['family'] for r in recs})
    print('records %d, families %d, score column %s' % (len(recs), len(fams), column))
    print()
    # ---- A. overview
    print('%-9s %-32s %7s %8s %15s %20s %12s' % ('benchmark', 'group', 'records', 'families', 'actual first', 'calls >=99%: right', 'temperature'))
    w = csv.writer(open(os.path.join(out, 'pooled_overview%s.csv' % tag), 'w', newline=''))
    w.writerow(['benchmark', 'group', 'records', 'families', 'actual_first', 'calls_99', 'calls_99_right', 'actual_at_1pct_or_less', 'temperature'])
    for b in benches + ['all']:
        for g in groups:
            rr = [r for r in recs if (b == 'all' or r['bench'] == b) and r['group'] == g]
            if not rr or (b == 'all' and len({r['bench'] for r in rr}) < 2):
                continue
            conf = [r for r in rr if r['ptop'] >= 0.99]
            right = sum(r['first'] for r in conf)
            T = fitted_temperature(rr)
            first = sum(r['first'] for r in rr)
            w.writerow([b, g, len(rr), len({r['family'] for r in rr}), first, len(conf), right, sum(r['pa'] <= 0.01 for r in rr), '%.2f' % T])
            print('%-9s %-32s %7d %8d %5d (%5.1f%%) %8d / %-4d (%5.1f%%) %9.2f' % (b, g, len(rr), len({r['family'] for r in rr}), first, 100.0 * first / len(rr),
                                                                           right, len(conf), 100.0 * right / max(1, len(conf)), T))
    # ---- B. reliability of the top probability
    print()
    print('how often the top candidate is the actual structure, by its stated probability')
    bins = [(0.0, 0.5), (0.5, 0.9), (0.9, 0.99), (0.99, 0.999), (0.999, 1.0001)]
    print('%-32s ' % 'group' + ' '.join('%17s' % ('%.3g-%.3g' % (a, min(b, 1))) for a, b in bins))
    for g in groups:
        rr = [r for r in recs if r['group'] == g]
        cells = []
        for a, b in bins:
            x = [r for r in rr if a <= r['ptop'] < b]
            cells.append('%4d/%-4d (%5.1f%%)' % (sum(r['first'] for r in x), len(x), 100.0 * sum(r['first'] for r in x) / len(x)) if x else '%17s' % '-')
        print('%-32s ' % g + ' '.join(cells))
    # ---- C. designs
    def pool_key(design, r):
        return {'own level': r['level'], 'by benchmark': r['bench'], 'by group': r['group'], 'one for all': 'all'}[design]

    designs = ['own level', 'by benchmark', 'by group', 'one for all']
    idx = {d: collections.defaultdict(list) for d in designs}
    for i, r in enumerate(recs):
        for d in designs:
            idx[d][pool_key(d, r)].append(i)
    wr = csv.writer(open(os.path.join(out, 'pooled_designs%s.csv' % tag), 'w', newline=''))
    wr.writerow(['target', 'design', 'benchmark_or_group', 'coverage', 'mean_set_size', 'single_structure_sets'])
    for target in (0.90, 0.80):
        rng = random.Random(seed)
        # per design, per reporting unit: sums over (split, family) of family means
        acc = {d: collections.defaultdict(lambda: [0.0, 0.0, 0.0, 0]) for d in designs}
        for _ in range(splits):
            f = fams[:]
            rng.shuffle(f)
            cal = set(f[:len(f) // 2])
            for d in designs:
                q = {}
                for key, ii in idx[d].items():
                    q[key] = threshold([(recs[i]['t'], recs[i]['family']) for i in ii if recs[i]['family'] in cal], target)
                per = collections.defaultdict(lambda: [0.0, 0.0, 0.0, 0])      # (unit, family) -> sums over records
                for r in recs:
                    if r['family'] in cal:
                        continue
                    th = q[pool_key(d, r)]
                    n = bisect.bisect_right(r['s'], th + 1e-12)
                    hit = r['t'] <= th + 1e-12
                    for unit in (r['bench'], r['group'], 'all'):
                        c = per[(unit, r['family'])]
                        c[0] += hit
                        c[1] += n
                        c[2] += n == 1
                        c[3] += 1
                for (unit, fam), c in per.items():
                    a = acc[d][unit]
                    a[0] += c[0] / c[3]
                    a[1] += c[1] / c[3]
                    a[2] += c[2] / c[3]
                    a[3] += 1
        print()
        print('target %.0f%%: coverage / mean set size / share of single-structure sets (%d random halves of the families)' % (100 * target, splits))
        print('%-32s ' % '' + ' '.join('%-22s' % d for d in designs))
        for unit in benches + groups + ['all']:
            cells = []
            for d in designs:
                a = acc[d][unit]
                cells.append('%5.1f%% %5.2f %5.1f%%   ' % (100 * a[0] / a[3], a[1] / a[3], 100 * a[2] / a[3]))
                wr.writerow([target, d, unit, '%.4f' % (a[0] / a[3]), '%.3f' % (a[1] / a[3]), '%.4f' % (a[2] / a[3])])
            print('%-32s ' % unit + ' '.join(cells))
        # ---- D. transfer between benchmarks
        print()
        print('target %.0f%%, transfer: threshold from the other benchmarks, families of the test benchmark left out' % (100 * target))
        print('%-9s | from every other record: coverage, set size | from the same groups only: coverage, set size' % 'test on')
        for b in benches:
            test = [r for r in recs if r['bench'] == b]
            tf = {r['family'] for r in test}
            out_cells = []
            for same in (False, True):
                if same:
                    qs = {g: threshold([(r['t'], r['family']) for r in recs if r['bench'] != b and r['family'] not in tf and r['group'] == g], target)
                          for g in {r['group'] for r in test}}
                else:
                    q0 = threshold([(r['t'], r['family']) for r in recs if r['bench'] != b and r['family'] not in tf], target)
                    qs = collections.defaultdict(lambda: q0)
                per = collections.defaultdict(lambda: [0.0, 0.0, 0])
                usable = 0
                for r in test:
                    th = qs[r['group']]
                    if th == float('inf') and same and not any(x['bench'] != b and x['family'] not in tf and x['group'] == r['group'] for x in recs):
                        continue
                    usable += 1
                    c = per[r['family']]
                    c[0] += r['t'] <= th + 1e-12
                    c[1] += bisect.bisect_right(r['s'], th + 1e-12)
                    c[2] += 1
                if per:
                    cov = sum(c[0] / c[2] for c in per.values()) / len(per)
                    size = sum(c[1] / c[2] for c in per.values()) / len(per)
                    out_cells.append('%5.1f%%  %5.2f  (%d records)' % (100 * cov, size, usable))
                    wr.writerow([target, 'transfer, same groups' if same else 'transfer, every other record', b, '%.4f' % cov, '%.3f' % size, ''])
                else:
                    out_cells.append('no records of these groups elsewhere')
            print('%-9s | %-42s | %s' % (b, out_cells[0], out_cells[1]))


if __name__ == '__main__':
    main(*sys.argv[1:8])
