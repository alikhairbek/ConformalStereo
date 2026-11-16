"""ConformalStereo, final method: which score, which calibration design, what it guarantees and what it costs.
Usage: final_method.py <pooled_candidates.csv> <out folder> [bootstrap=4000] [splits=500] [seed=20261114]

The method. For every candidate of a record: s = (log-likelihood of the best candidate - log-likelihood of the
candidate) / number of shifts. A workflow (benchmark, method, kind of geometry, NMR level) gets a cut-off: the
family-weighted conformal quantile of s of the actual structures of its calibration molecules. The set of a new molecule
is every candidate with s at or below the cut-off; the best candidate is always in it.

Evaluation. Sets of a family are made with a cut-off calibrated on the other families of the workflow (one family
out). Numbers are averaged per family within a workflow and then over these units; intervals come from resampling
families. Coverage under this scheme and under random halves holds by construction: it checks the implementation.
What the data can say is how large the sets are, where the cut-off may come from, and who is covered; sections 4, 6
and 7 are about that. Tables go to the out folder, the log to standard output.
"""
import sys
import os
import csv
import math
import collections

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cslib as cs                                                                      # noqa: E402
from cslib import INF                                                                   # noqa: E402

table, out = sys.argv[1:3]
BOOT = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
SPLITS = int(sys.argv[4]) if len(sys.argv) > 4 else 500
SEED = int(sys.argv[5]) if len(sys.argv) > 5 else 20261114
TARGETS = (0.95, 0.90, 0.80)
FINAL = 'likelihood ratio per shift'
RATIOS = ['likelihood ratio per shift', 'likelihood ratio / sqrt(N)', 'likelihood ratio']
os.makedirs(out, exist_ok=True)

recs = cs.load(table)
BENCH = list(collections.OrderedDict((r['bench'], 1) for r in recs))
GROUPS = list(collections.OrderedDict((r['group'], 1) for r in recs))


def workflows(rr):
    d = collections.OrderedDict()
    for r in rr:
        d.setdefault(r['wf'], []).append(r)
    return d


by_wf = workflows(recs)
print('records %d | candidates %d | families %d | workflows %d' % (len(recs), sum(r['K'] for r in recs), len({r['family'] for r in recs}), len(by_wf)))


def writer(name, header):
    fh = open(os.path.join(out, name), 'w', newline='')
    w = csv.writer(fh)
    w.writerow(header)
    return w


def needed(target):
    """smallest number of calibration families with which the target can be reached"""
    g = 1
    while g / (g + 1.0) < target - 1e-12:
        g += 1
    return g


def reachable(G, target, held_out=1):
    return G - held_out >= needed(target)


def bench_reaches(rr, b, target):
    return reachable(min(len({x['family'] for x in v}) for wf, v in workflows(rr).items() if wf[0] == b), target)


# ------------------------------------------------------------------ 1. the benchmarks and their overconfidence
print()
print('1. THE BENCHMARKS: how often the top candidate is right, and how often a confident call is right')
w1 = writer('final_benchmarks.csv', ['benchmark', 'kind_of_calculation', 'workflows', 'records', 'families', 'mean_candidates', 'top_right_records', 'top_right_share_records',
                                    'top_right_share_families', 'calls_99', 'calls_99_right', 'calls_999', 'calls_999_right', 'temperature'])
print('%-9s %-32s %3s %7s %5s %5s | %-22s %7s | %-17s | %-17s | %5s' % ('benchmark', 'kind of calculation', 'wf', 'records', 'fam.', 'cand.', 'top candidate right', 'by fam.',
                                                                    'calls >= 99%', 'calls >= 99.9%', 'T'))
for b in BENCH:
    groups = list(collections.OrderedDict((r['group'], 1) for r in recs if r['bench'] == b))
    for g in groups + (['all'] if len(groups) > 1 else []):
        rr = [r for r in recs if r['bench'] == b and (g == 'all' or r['group'] == g)]
        um = cs.unit_means(rr, [(float(r['first']), float(r['K'])) for r in rr])
        topf = sum(v[0] for v in um.values()) / len(um)
        Kf = sum(v[1] for v in um.values()) / len(um)
        c99 = [r for r in rr if r['ptop'] >= 0.99]
        c999 = [r for r in rr if r['ptop'] >= 0.999]
        T = cs.Fit(rr).temperature()
        row = [b, g, len({r['wf'] for r in rr}), len(rr), len({r['family'] for r in rr}), '%.8f' % Kf, sum(r['first'] for r in rr),
               '%.8f' % (sum(r['first'] for r in rr) / len(rr)), '%.8f' % topf, len(c99), sum(r['first'] for r in c99), len(c999), sum(r['first'] for r in c999), '%.8f' % T]
        w1.writerow(row)
        print('%-9s %-32s %3d %7d %5d %5.2f | %4d of %-4d (%5.1f%%) %6.1f%% | %4d of %-4d %5.1f%% | %4d of %-4d %5.1f%% | %5.2f' % (
            b, g, row[2], len(rr), row[4], Kf, row[6], len(rr), 100 * row[6] / len(rr), 100 * topf, row[10], row[9], 100 * row[10] / max(row[9], 1),
            row[12], row[11], 100 * row[12] / max(row[11], 1), T))
c99 = [r for r in recs if r['ptop'] >= 0.99]
c999 = [r for r in recs if r['ptop'] >= 0.999]
print('all records: top candidate right %d of %d (%.1f%%) | calls >= 99%%: %d of %d right (%.1f%%) | calls >= 99.9%%: %d of %d right (%.1f%%)' % (
    sum(r['first'] for r in recs), len(recs), 100 * sum(r['first'] for r in recs) / len(recs), sum(r['first'] for r in c99), len(c99), 100 * sum(r['first'] for r in c99) / len(c99),
    sum(r['first'] for r in c999), len(c999), 100 * sum(r['first'] for r in c999) / len(c999)))
fb = collections.defaultdict(set)
for r in recs:
    fb[r['family']].add(r['bench'])
print('families by benchmark: ' + ', '.join('%s %d' % (' + '.join(k), v) for k, v in collections.Counter(tuple(sorted(v)) for v in fb.values()).most_common()))

# reliability: stated probability of the top candidate against how often it is right
BINS = [(0.0, 0.5, 'below 50%'), (0.5, 0.9, '50 to 90%'), (0.9, 0.99, '90 to 99%'), (0.99, 0.999, '99 to 99.9%'), (0.999, 0.99999, '99.9 to 99.999%'), (0.99999, 1.01, '99.999% or more')]
w1b = writer('final_reliability.csv', ['kind_of_calculation', 'stated_probability_of_top_candidate', 'records', 'families', 'mean_stated', 'right', 'share_right', 'share_right_lo', 'share_right_hi'])
print()
print('reliability of the stated probability of the top candidate (share right, with 95% interval from resampling families)')
print('%-32s ' % 'kind of calculation' + ' '.join('%-22s' % b[2] for b in BINS))
rng = np.random.default_rng(SEED + 1)
for g in GROUPS + ['all']:
    cells = []
    for lo, hi, name in BINS:
        rr = [r for r in recs if (g == 'all' or r['group'] == g) and lo <= r['ptop'] < hi]
        if not rr:
            cells.append('-')
            continue
        fam = collections.defaultdict(lambda: np.zeros(2))
        for r in rr:
            fam[r['family']] += np.array([float(r['first']), 1.0])
        A = np.array(list(fam.values()))
        S = A[rng.integers(0, len(A), size=(BOOT, len(A)))].sum(axis=1)
        q = np.quantile(S[:, 0] / S[:, 1], [0.025, 0.975])
        share = sum(r['first'] for r in rr) / len(rr)
        w1b.writerow([g, name, len(rr), len(fam), '%.6f' % (sum(r['ptop'] for r in rr) / len(rr)), sum(r['first'] for r in rr), '%.8f' % share, '%.8f' % q[0], '%.8f' % q[1]])
        cells.append('%3.0f%% (%2.0f-%3.0f) n=%-4d' % (100 * share, 100 * q[0], 100 * q[1], len(rr)))
    print('%-32s ' % g + ' '.join('%-22s' % c for c in cells))


# ------------------------------------------------------------------ 2. leaving one family out
POOLS = {'own workflow': lambda r: r['wf'], 'whole benchmark': lambda r: r['bench'], 'same kind of calculation': lambda r: r['group'], 'everything': lambda r: 0}


def weighted_quantile(items, cnt, G, target, skip=None):
    """items sorted by score: [(score, family)]; the family `skip` is left out; G families remain"""
    acc = 0.0
    for s, f in items:
        if f == skip:
            continue
        acc += 1.0 / cnt[f] / (G + 1)
        if acc >= target - 1e-12:
            return s
    return INF


def lofo_thresholds(rr, score, target, pool='own workflow'):
    """for every record the cut-off calibrated on the records of its pool that are not of its family"""
    thr = {}
    if pool in ('other benchmarks', 'same kind of calculation, other benchmarks'):
        same_kind = pool != 'other benchmarks'
        keys = collections.OrderedDict(((r['bench'], r['group'] if same_kind else None), 1) for r in rr)
        for b, g in keys:
            items = sorted(((r['t'][score], r['family']) for r in rr if r['bench'] != b and (g is None or r['group'] == g)), key=lambda x: x[0])
            cnt = collections.Counter(f for _, f in items)
            for f in {r['family'] for r in rr if r['bench'] == b and (g is None or r['group'] == g)}:
                thr[(b, g, f)] = weighted_quantile(items, cnt, len(cnt) - (1 if f in cnt else 0), target, skip=f) if items else None
        return [thr[(r['bench'], r['group'] if same_kind else None, r['family'])] for r in rr]
    key = POOLS[pool]
    pools = collections.defaultdict(list)
    for r in rr:
        pools[key(r)].append(r)
    for pk, pr in pools.items():
        items = sorted(((r['t'][score], r['family']) for r in pr), key=lambda x: x[0])
        cnt = collections.Counter(f for _, f in items)
        for f in cnt:
            thr[(pk, f)] = weighted_quantile(items, cnt, len(cnt) - 1, target, skip=f)
    return [thr[(key(r), r['family'])] for r in rr]


def sets(rr, score, thr):
    cov, size = [], []
    for r, t in zip(rr, thr):
        c, s = cs.apply_set(r, score, t)
        cov.append(c)
        size.append(s)
    return cov, size


def cell(x):
    if x is None:
        return '%-34s' % '    -'
    return '%5.1f%% %5.2f %4.0f%% (%3.0f%%) %4.0f%%     ' % (100 * x['coverage'], x['size'], 100 * x['single'], 100 * x['single_right'] if x['single'] > 0 else 0, 100 * x['every'])


def by_bench(rr, cov, size, target, boot=0, seed=SEED + 2):
    """summary per benchmark; a benchmark whose workflows cannot reach the target is left out"""
    res = {}
    for b in BENCH:
        idx = [i for i, r in enumerate(rr) if r['bench'] == b]
        if not idx:
            continue
        res[b] = cs.summarize([rr[i] for i in idx], [cov[i] for i in idx], [size[i] for i in idx], boot=boot, seed=seed) if bench_reaches(rr, b, target) else None
    return res


def paired(rr, a, b, boot=BOOT, seed=SEED + 3):
    """mean over units of (a - b) for per-record values a, b; interval from resampling families"""
    um = cs.unit_means(rr, [(float(x) - float(y),) for x, y in zip(a, b)])
    fam = collections.defaultdict(lambda: np.zeros(2))
    for (wf, f), v in um.items():
        fam[f] += np.array([v[0], 1.0])
    A = np.array(list(fam.values()))
    rg = np.random.default_rng(seed)
    S = A[rg.integers(0, len(A), size=(boot, len(A)))].sum(axis=1)
    q = np.quantile(S[:, 0] / S[:, 1], [0.025, 0.975])
    return A[:, 0].sum() / A[:, 1].sum(), float(q[0]), float(q[1])


# nested temperature: for every workflow and every family left out, the temperature of the other families
TEMP = {}
for wf, rr in by_wf.items():
    fit = cs.Fit(rr)
    fam = np.array([r['family'] for r in rr])
    for f in sorted(set(fam)):
        TEMP[(wf, f)] = fit.temperature(fam != f)


def tempered_stated(target):
    cov, size = [], []
    for r in recs:
        c, s = cs.stated_set(cs.softmax(r['ell'] / TEMP[(r['wf'], r['family'])]), r['k'], target)
        cov.append(c)
        size.append(s)
    return cov, size


def tempered_cumulative(target):
    """conformal sets with the cumulative-probability score at the temperature fitted on the calibration families"""
    cov, size = [None] * len(recs), [None] * len(recs)
    pos = {id(r): i for i, r in enumerate(recs)}
    for wf, rr in by_wf.items():
        for f in sorted({r['family'] for r in rr}):
            beta = 1.0 / TEMP[(wf, f)]
            thr = cs.threshold([(float(cs.tail_mass_score(r['ell'], beta)[r['k']]), r['family']) for r in rr if r['family'] != f], target)
            for r in rr:
                if r['family'] == f:
                    s = cs.tail_mass_score(r['ell'], beta)
                    cov[pos[id(r)]] = bool(s[r['k']] <= thr)
                    size[pos[id(r)]] = int((s <= thr).sum())
    return cov, size


# ------------------------------------------------------------------ 3. scores and baselines
print()
print('3. SCORES AND BASELINES, one family out, cut-off from the own workflow')
print('   cells: coverage | mean set size | single-structure sets (share of them right) | sets with every candidate')
w3 = writer('final_scores.csv', ['target', 'kind', 'method', 'benchmark', 'coverage', 'coverage_lo', 'coverage_hi', 'mean_set_size', 'size_lo', 'size_hi', 'single_structure_sets',
                                 'single_structure_sets_right', 'sets_with_every_candidate', 'mean_share_of_candidates', 'mean_candidates'])
RES = {}                                                 # (target, method) -> per record (cov, size)
OTHERS = []
for target in TARGETS:
    methods = collections.OrderedDict()
    methods[('baseline', 'top candidate alone')] = ([r['first'] for r in recs], [1] * len(recs))
    st = [cs.stated_set(r['p'], r['k'], target) for r in recs]
    methods[('baseline', 'stated probabilities add up to the target')] = ([x[0] for x in st], [x[1] for x in st])
    methods[('baseline', 'the same after fitting a temperature')] = tempered_stated(target)
    for s in cs.SCORES:
        methods[('conformal', s)] = sets(recs, s, lofo_thresholds(recs, s, target))
    methods[('conformal', 'cumulative probability, fitted temperature')] = tempered_cumulative(target)
    OTHERS = [n for (k, n) in methods if n not in (FINAL, 'top candidate alone', 'stated probabilities add up to the target')]
    print()
    print('target %d%%' % (100 * target))
    print('%-44s ' % 'method' + ' '.join('%-34s' % b for b in BENCH))
    for (kind, name), (cov, size) in methods.items():
        RES[(target, name)] = (cov, size)
        res = by_bench(recs, cov, size, target, boot=BOOT)
        print('%-44s ' % name[:44] + ' '.join(cell(res[b]) for b in BENCH))
        for b in BENCH:
            x = res[b]
            if x:
                w3.writerow([target, kind, name, b, '%.8f' % x['coverage'], '%.8f' % x['coverage_lo'], '%.8f' % x['coverage_hi'], '%.8f' % x['size'], '%.8f' % x['size_lo'],
                             '%.8f' % x['size_hi'], '%.8f' % x['single'], '%.8f' % x['single_right'], '%.8f' % x['every'], '%.8f' % x['share'], '%.8f' % x['K']])
    res = by_bench(recs, *RES[(target, FINAL)], target, boot=BOOT)
    print('   %s, 95%% intervals from resampling families: ' % FINAL + ' | '.join('%s size %.2f (%.2f to %.2f), single %.0f%% (%.0f to %.0f)' % (
        b, x['size'], x['size_lo'], x['size_hi'], 100 * x['single'], 100 * x['single_lo'], 100 * x['single_hi']) for b, x in res.items() if x))

print()
print('mean set size of each alternative minus that of the %s score, paired by family (95%% interval from resampling families; * = interval excludes zero)' % FINAL)
w3b = writer('final_scores_paired.csv', ['target', 'method', 'benchmark', 'difference_of_mean_set_size', 'lo', 'hi', 'difference_of_coverage', 'cov_lo', 'cov_hi'])
for target in TARGETS:
    print('target %d%%' % (100 * target))
    for name in OTHERS:
        cells = []
        for b in BENCH:
            if not bench_reaches(recs, b, target):
                cells.append('%-26s' % '   -')
                continue
            idx = [i for i, r in enumerate(recs) if r['bench'] == b]
            rr = [recs[i] for i in idx]
            d = paired(rr, [RES[(target, name)][1][i] for i in idx], [RES[(target, FINAL)][1][i] for i in idx])
            dc = paired(rr, [RES[(target, name)][0][i] for i in idx], [RES[(target, FINAL)][0][i] for i in idx])
            w3b.writerow([target, name, b, '%.8f' % d[0], '%.8f' % d[1], '%.8f' % d[2], '%.8f' % dc[0], '%.8f' % dc[1], '%.8f' % dc[2]])
            cells.append('%+5.2f (%+5.2f to %+5.2f)%s' % (d[0], d[1], d[2], '*' if d[1] > 0 or d[2] < 0 else ' '))
        print('   %-44s ' % name[:44] + ' '.join('%-9s %s' % (b, c) for b, c in zip(BENCH, cells)))

# ------------------------------------------------------------------ 4. where the cut-off comes from
print()
print('4. WHERE THE CUT-OFF COMES FROM (%s score; the family of the record is always left out): coverage (95%% interval) and mean set size' % FINAL)
w4 = writer('final_designs.csv', ['target', 'cutoff_from', 'benchmark', 'records', 'coverage', 'coverage_lo', 'coverage_hi', 'mean_set_size', 'single_structure_sets', 'sets_with_every_candidate'])
DESIGNS = ['own workflow', 'whole benchmark', 'same kind of calculation', 'same kind of calculation, other benchmarks', 'everything', 'other benchmarks']
DES = {}
for target in (0.90, 0.80):
    print('target %d%%' % (100 * target))
    print('%-44s ' % 'cut-off from' + ' '.join('%-24s' % b for b in BENCH))
    for pool in DESIGNS:
        thr = lofo_thresholds(recs, FINAL, target, pool)
        keep = [i for i, t in enumerate(thr) if t is not None]
        cov, size = sets([recs[i] for i in keep], FINAL, [thr[i] for i in keep])
        DES[(target, pool)] = (keep, cov, size)
        cells = []
        for b in BENCH:
            idx = [n for n, i in enumerate(keep) if recs[i]['bench'] == b]
            if not idx:
                cells.append('   -')
                continue
            x = cs.summarize([recs[keep[n]] for n in idx], [cov[n] for n in idx], [size[n] for n in idx], boot=BOOT, seed=SEED + 4)
            w4.writerow([target, pool, b, len(idx), '%.8f' % x['coverage'], '%.8f' % x['coverage_lo'], '%.8f' % x['coverage_hi'], '%.8f' % x['size'], '%.8f' % x['single'], '%.8f' % x['every']])
            cells.append('%5.1f%% (%2.0f-%3.0f) %5.2f' % (100 * x['coverage'], 100 * x['coverage_lo'], 100 * x['coverage_hi'], x['size']))
        print('%-44s ' % pool + ' '.join('%-24s' % c for c in cells))
print('coverage per workflow under each design: lowest, highest, and number of workflows more than five points off the target')
for target in (0.90, 0.80):
    for pool in DESIGNS:
        keep, cov, size = DES[(target, pool)]
        per = []
        g = collections.defaultdict(list)
        for n, i in enumerate(keep):
            g[recs[i]['wf']].append(n)
        for wf, idx in g.items():
            x = cs.summarize([recs[keep[n]] for n in idx], [cov[n] for n in idx], [size[n] for n in idx])
            per.append((x['coverage'], wf))
        per.sort()
        print('   %d%%  %-44s lowest %5.1f%% (%s) | highest %5.1f%% | below target by more than 5: %2d of %d | above by more than 5: %2d' % (
            100 * target, pool, 100 * per[0][0], ' / '.join(per[0][1]), 100 * per[-1][0], sum(1 for c, _ in per if c < target - 0.05), len(per), sum(1 for c, _ in per if c > target + 0.05)))

# ------------------------------------------------------------------ 5. families or records as the units
print()
print('5. FAMILIES OR RECORDS AS THE UNITS (%s score, cut-off from the own workflow): coverage by family, mean set size, coverage by record' % FINAL)
w5 = writer('final_units.csv', ['target', 'cutoff', 'left_out', 'benchmark', 'coverage_by_family', 'mean_set_size', 'coverage_by_record'])
pos = {id(r): i for i, r in enumerate(recs)}
for target in (0.90, 0.80):
    variants = collections.OrderedDict()
    variants[('families weighted', 'one family')] = RES[(target, FINAL)]
    cov_b, size_b, cov_c, size_c = [None] * len(recs), [None] * len(recs), [None] * len(recs), [None] * len(recs)
    for wf, rr in by_wf.items():
        ts = [r['t'][FINAL] for r in rr]
        for j, r in enumerate(rr):
            cov_b[pos[id(r)]], size_b[pos[id(r)]] = cs.apply_set(r, FINAL, cs.record_threshold(ts[:j] + ts[j + 1:], target))
            cov_c[pos[id(r)]], size_c[pos[id(r)]] = cs.apply_set(r, FINAL, cs.record_threshold([x['t'][FINAL] for x in rr if x['family'] != r['family']], target))
    variants[('records unweighted', 'one record')] = (cov_b, size_b)
    variants[('records unweighted', 'one family')] = (cov_c, size_c)
    print('target %d%%' % (100 * target))
    for (thr_kind, left), (cov, size) in variants.items():
        cells = []
        for b in BENCH:
            idx = [i for i, r in enumerate(recs) if r['bench'] == b]
            x = cs.summarize([recs[i] for i in idx], [cov[i] for i in idx], [size[i] for i in idx])
            rec_cov = sum(cov[i] for i in idx) / len(idx)
            w5.writerow([target, thr_kind, left, b, '%.8f' % x['coverage'], '%.8f' % x['size'], '%.8f' % rec_cov])
            cells.append('%5.1f%% %5.2f %5.1f%%' % (100 * x['coverage'], x['size'], 100 * rec_cov))
        print('   %-19s %-14s ' % (thr_kind, left + ' out') + ' '.join('%-8s %s' % (b, c) for b, c in zip(BENCH, cells)))

# ------------------------------------------------------------------ 6. who is covered
print()
print('6. WHO IS COVERED (cut-off from the own workflow, one family out; 95% intervals from resampling families)')
w6 = writer('final_conditional.csv', ['target', 'score', 'split_by', 'class', 'records', 'families', 'top_candidate_right', 'coverage', 'coverage_lo', 'coverage_hi', 'mean_set_size', 'mean_candidates'])
nclass = lambda r: 'fewer than 25 shifts' if r['N'] < 25 else '25 to 39 shifts' if r['N'] < 40 else '40 shifts or more'      # noqa: E731
NCLASSES = ['fewer than 25 shifts', '25 to 39 shifts', '40 shifts or more']
for target in (0.90, 0.80):
    ok = [i for i, r in enumerate(recs) if reachable(len({x['family'] for x in by_wf[r['wf']]}), target)]
    print('target %d%%: coverage and mean set size by class, for the three likelihood-ratio scores' % (100 * target))
    for title, fn, classes in [('number of candidates', lambda r: cs.kclass(r['K']), cs.KCLASSES), ('number of shifts', nclass, NCLASSES)]:
        for c in classes:
            idx = [i for i in ok if fn(recs[i]) == c]
            cells = []
            for s in RATIOS:
                cov, size = RES[(target, s)]
                x = cs.summarize([recs[i] for i in idx], [cov[i] for i in idx], [size[i] for i in idx], boot=BOOT, seed=SEED + 5)
                w6.writerow([target, s, title, c, x['records'], x['families'], '%.8f' % x['top'], '%.8f' % x['coverage'], '%.8f' % x['coverage_lo'], '%.8f' % x['coverage_hi'],
                             '%.8f' % x['size'], '%.8f' % x['K']])
                cells.append('%5.1f%% (%2.0f-%3.0f) %5.2f' % (100 * x['coverage'], 100 * x['coverage_lo'], 100 * x['coverage_hi'], x['size']))
            print('   %-21s %-21s rec %4d fam %3d top %5.1f%% of %5.2f | per shift %s | / sqrt(N) %s | plain %s' % (title, c, x['records'], x['families'], 100 * x['top'], x['K'], *cells))
    cov, size = RES[(target, FINAL)]
    print('target %d%%: %s score' % (100 * target, FINAL))
    for title, fn, classes in [('record without data problems', lambda r: 'yes' if r['clean'] else 'no', ['yes', 'no']), ('kind of calculation', lambda r: r['group'], GROUPS),
                               ('size of the set', None, ['one structure', 'more than one, not all', 'every candidate'])]:
        for c in classes:
            if fn is None:
                test = {'one structure': lambda r, s: s == 1, 'more than one, not all': lambda r, s: 1 < s < r['K'], 'every candidate': lambda r, s: s == r['K']}[c]
                idx = [i for i in ok if test(recs[i], size[i])]
            else:
                idx = [i for i in ok if fn(recs[i]) == c]
            if not idx:
                continue
            x = cs.summarize([recs[i] for i in idx], [cov[i] for i in idx], [size[i] for i in idx], boot=BOOT, seed=SEED + 6)
            w6.writerow([target, FINAL, title, c, x['records'], x['families'], '%.8f' % x['top'], '%.8f' % x['coverage'], '%.8f' % x['coverage_lo'], '%.8f' % x['coverage_hi'],
                         '%.8f' % x['size'], '%.8f' % x['K']])
            print('   %-28s %-32s rec %4d fam %3d | top right %5.1f%% | coverage %5.1f%% (%4.1f to %5.1f) | set %5.2f of %5.2f' % (
                title, c, x['records'], x['families'], 100 * x['top'], 100 * x['coverage'], 100 * x['coverage_lo'], 100 * x['coverage_hi'], x['size'], x['K']))

# ------------------------------------------------------------------ 7. stress tests of exchangeability
print()
print('7. STRESS TESTS: calibration molecules that differ from the test molecules in a known way (same workflow)')
print('   cells: coverage (95% interval from resampling test families), mean set size; workflows used / test families')
w7s = writer('final_stress.csv', ['test', 'calibrated_on', 'tested_on', 'target', 'score', 'workflows', 'test_families', 'test_records', 'coverage', 'coverage_lo', 'coverage_hi', 'mean_set_size',
                                  'mean_candidates'])
in_joc = {f for f, v in fb.items() if 'JOC2015' in v}
in_2019 = {f for f, v in fb.items() if 'Cam2019' in v}


def fam_stat(rr, what):
    d = collections.defaultdict(list)
    for r in rr:
        d[r['family']].append(r[what])
    return {f: float(np.median(v)) for f, v in d.items()}


def side_source_stx(wf, rr):
    return {f: (f in in_joc) for f in {r['family'] for r in rr}} if wf[0] == 'STX150' else None


def side_source_cam(wf, rr):
    return {f: (f not in in_2019) for f in {r['family'] for r in rr}} if wf[0] == 'Cam2017' else None


def side_size(wf, rr):
    n = fam_stat(rr, 'N')
    med = float(np.median(list(n.values())))
    return {f: v <= med for f, v in n.items()}


def side_cand(wf, rr):
    k = fam_stat(rr, 'K')
    return {f: v <= 4 for f, v in k.items()}


def stress(name, label_true, label_false, side_fn, reverse=False):
    for target in (0.90, 0.80):
        used, tested = 0, []
        res = {s: ([], []) for s in RATIOS}
        for wf, rr in by_wf.items():
            side = side_fn(wf, rr)
            if side is None:
                continue
            cal = {f for f, v in side.items() if v != reverse}
            test = [r for r in rr if r['family'] not in cal]
            if len(cal) < needed(target) or len({r['family'] for r in test}) < 3:
                continue
            used += 1
            tested += test
            for s in RATIOS:
                thr = cs.threshold([(r['t'][s], r['family']) for r in rr if r['family'] in cal], target)
                for r in test:
                    c, z = cs.apply_set(r, s, thr)
                    res[s][0].append(c)
                    res[s][1].append(z)
        a, b = (label_false, label_true) if reverse else (label_true, label_false)
        if not used:
            print('   %-22s %-30s -> %-30s %d%%: no workflow has enough calibration families' % (name, a, b, 100 * target))
            continue
        cells = []
        for s in RATIOS:
            x = cs.summarize(tested, res[s][0], res[s][1], boot=BOOT, seed=SEED + 9)
            w7s.writerow([name, a, b, target, s, used, x['families'], x['records'], '%.8f' % x['coverage'], '%.8f' % x['coverage_lo'], '%.8f' % x['coverage_hi'], '%.8f' % x['size'], '%.8f' % x['K']])
            cells.append('%5.1f%% (%2.0f-%3.0f) %5.2f' % (100 * x['coverage'], 100 * x['coverage_lo'], 100 * x['coverage_hi'], x['size']))
        print('   %-22s %-30s -> %-30s %d%%: per shift %s | / sqrt(N) %s | plain %s | %2d wf, %3d fam, of %.2f' % (name, a, b, 100 * target, *cells, used, x['families'], x['K']))


stress('collection (STX150)', 'the 2015 collection', 'molecules added later', side_source_stx)
stress('collection (STX150)', 'the 2015 collection', 'molecules added later', side_source_stx, reverse=True)
stress('collection (Cam2017)', 'molecules not taken up in 2019', 'molecules taken up in 2019', side_source_cam)
stress('collection (Cam2017)', 'molecules not taken up in 2019', 'molecules taken up in 2019', side_source_cam, reverse=True)
stress('size', 'the half with fewer shifts', 'the half with more shifts', side_size)
stress('size', 'the half with fewer shifts', 'the half with more shifts', side_size, reverse=True)
stress('candidates', 'at most four candidates', 'five or more candidates', side_cand)
stress('candidates', 'at most four candidates', 'five or more candidates', side_cand, reverse=True)

# ------------------------------------------------------------------ 8. the workflows
print()
print('8. THE WORKFLOWS: cut-offs for use (all families) and what the sets cost (one family out)')
w8 = writer('final_workflows.csv', ['benchmark', 'method', 'geometry', 'level', 'records', 'families', 'mean_candidates', 'median_shifts', 'top_candidate_right', 'temperature'] +
            [x for t in TARGETS for x in ('cutoff_per_shift_%d' % (100 * t), 'least_likelihood_ratio_per_shift_%d' % (100 * t), 'coverage_%d' % (100 * t), 'mean_set_size_%d' % (100 * t),
                                          'single_structure_sets_%d' % (100 * t), 'sets_with_every_candidate_%d' % (100 * t))])
WF = []
for wf, rr in by_wf.items():
    idx = [pos[id(r)] for r in rr]
    G = len({r['family'] for r in rr})
    row = dict(wf=wf, n=len(rr), G=G, T=cs.Fit(rr).temperature(), N=float(np.median([r['N'] for r in rr])))
    base = cs.summarize(rr, [r['first'] for r in rr], [1] * len(rr))
    row['top'], row['K'] = base['top'], base['K']
    for target in TARGETS:
        cut = cs.threshold([(r['t'][FINAL], r['family']) for r in rr], target)
        x = cs.summarize(rr, [RES[(target, FINAL)][0][i] for i in idx], [RES[(target, FINAL)][1][i] for i in idx]) if reachable(G, target) else None
        row[target] = (cut, x)
    WF.append(row)
    line = list(wf) + [row['n'], G, '%.8f' % row['K'], '%.0f' % row['N'], '%.8f' % row['top'], '%.8f' % row['T']]
    for target in TARGETS:
        cut, x = row[target]
        line += ['' if cut == INF else '%.8f' % cut, '' if cut == INF else '%.8f' % math.exp(-cut)] + (['', '', '', ''] if x is None else
                                                                                                     ['%.8f' % x['coverage'], '%.8f' % x['size'], '%.8f' % x['single'], '%.8f' % x['every']])
    w8.writerow(line)
print('%d workflows. Cut-off = largest allowed distance of the log-likelihood from the best candidate, per shift.' % len(WF))
for target in TARGETS:
    cuts = [r[target][0] for r in WF]
    fin = [c for c in cuts if c != INF]
    posv = [c for c in fin if c > 0]
    print('   %d%%: cannot be reached (fewer than %d families) %2d | cut-off of zero (the best candidate alone) %2d | positive %2d: median %.3f, range %.3f to %.3f' % (
        100 * target, needed(target), sum(1 for c in cuts if c == INF), sum(1 for c in fin if c == 0), len(posv),
        float(np.median(posv)) if posv else 0.0, min(posv) if posv else 0.0, max(posv) if posv else 0.0))


def spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)).astype(float), np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    den = math.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den > 0 else float('nan')


print('rank correlation between the share of top candidates right and the share of candidates kept, over the workflows of a benchmark')
for b in BENCH:
    rows = [r for r in WF if r['wf'][0] == b]
    cells = []
    for target in TARGETS:
        v = [(r['top'], r[target][1]['share']) for r in rows if r[target][1] is not None]
        if len(v) >= 4:
            cells.append('%d%%: %+.2f' % (100 * target, spearman([x[0] for x in v], [x[1] for x in v])))
    print('   %-8s %2d workflows | %s' % (b, len(rows), ' | '.join(cells)))

print('kinds of calculation: top candidate right, then per target coverage, mean set size (95% interval), single-structure sets and median cut-off')
w8b = writer('final_kinds.csv', ['benchmark', 'method', 'geometry', 'workflows', 'records', 'families', 'mean_candidates', 'top_candidate_right'] +
             [x for t in TARGETS for x in ('coverage_%d' % (100 * t), 'mean_set_size_%d' % (100 * t), 'size_lo_%d' % (100 * t), 'size_hi_%d' % (100 * t), 'single_structure_sets_%d' % (100 * t),
                                           'median_cutoff_per_shift_%d' % (100 * t))])
KINDS = list(collections.OrderedDict(((r['bench'], r['method'], r['geo']), 1) for r in recs))
for kd in KINDS:
    idx = [i for i, r in enumerate(recs) if (r['bench'], r['method'], r['geo']) == kd]
    rr = [recs[i] for i in idx]
    G = min(len({x['family'] for x in v}) for wf, v in by_wf.items() if wf[:3] == kd)
    base = cs.summarize(rr, [r['first'] for r in rr], [1] * len(rr))
    line = list(kd) + [len({r['wf'] for r in rr}), len(rr), len({r['family'] for r in rr}), '%.8f' % base['K'], '%.8f' % base['top']]
    txt = []
    for target in TARGETS:
        cuts = [r[target][0] for r in WF if r['wf'][:3] == kd and r[target][0] != INF]
        if reachable(G, target):
            x = cs.summarize(rr, [RES[(target, FINAL)][0][i] for i in idx], [RES[(target, FINAL)][1][i] for i in idx], boot=BOOT, seed=SEED + 7)
            line += ['%.8f' % x['coverage'], '%.8f' % x['size'], '%.8f' % x['size_lo'], '%.8f' % x['size_hi'], '%.8f' % x['single'], '%.8f' % float(np.median(cuts))]
            txt.append('%d%%: %5.1f%% %.2f (%.2f to %.2f) single %3.0f%% cut %.3f' % (100 * target, 100 * x['coverage'], x['size'], x['size_lo'], x['size_hi'], 100 * x['single'], float(np.median(cuts))))
        else:
            line += ['', '', '', '', '', '']
    w8b.writerow(line)
    print('   %-8s %-9s %-7s wf %2d rec %4d fam %2d cand %5.2f top %5.1f%% | %s' % (kd[0], kd[1], kd[2], line[3], line[4], line[5], base['K'], 100 * base['top'], ' | '.join(txt)))
print('paired differences of mean set size between kinds of calculation of one benchmark (families present in both; shown when the 95% interval excludes zero)')
w8c = writer('final_kinds_paired.csv', ['benchmark', 'kind_a', 'kind_b', 'target', 'families', 'top_right_a', 'top_right_b', 'size_a', 'size_b', 'difference', 'lo', 'hi'])
for b in BENCH:
    kinds = [k for k in KINDS if k[0] == b]
    for i in range(len(kinds)):
        for j in range(i + 1, len(kinds)):
            for target in TARGETS:
                fam_val = []
                for kd in (kinds[i], kinds[j]):
                    idx = [n for n, r in enumerate(recs) if (r['bench'], r['method'], r['geo']) == kd]
                    G = min(len({x['family'] for x in v}) for wf, v in by_wf.items() if wf[:3] == kd)
                    if not reachable(G, target):
                        fam_val = None
                        break
                    um = cs.unit_means([recs[n] for n in idx], [(float(RES[(target, FINAL)][1][n]), float(recs[n]['first'])) for n in idx])
                    fv = collections.defaultdict(list)
                    for (wf, f), v in um.items():
                        fv[f].append(v)
                    fam_val.append({f: (sum(x[0] for x in v) / len(v), sum(x[1] for x in v) / len(v)) for f, v in fv.items()})
                if not fam_val:
                    continue
                common = sorted(set(fam_val[0]) & set(fam_val[1]))
                if len(common) < 8:
                    continue
                d = np.array([fam_val[0][f][0] - fam_val[1][f][0] for f in common])
                rg = np.random.default_rng(SEED + 8)
                bt = d[rg.integers(0, len(d), size=(BOOT, len(d)))].mean(axis=1)
                lo, hi = np.quantile(bt, [0.025, 0.975])
                a_top, b_top = np.mean([fam_val[0][f][1] for f in common]), np.mean([fam_val[1][f][1] for f in common])
                a_size, b_size = np.mean([fam_val[0][f][0] for f in common]), np.mean([fam_val[1][f][0] for f in common])
                w8c.writerow([b, ' '.join(kinds[i][1:]), ' '.join(kinds[j][1:]), target, len(common), '%.8f' % a_top, '%.8f' % b_top, '%.8f' % a_size, '%.8f' % b_size,
                              '%.8f' % d.mean(), '%.8f' % lo, '%.8f' % hi])
                if lo > 0 or hi < 0:
                    print('   %-8s %d%%  %-16s top %5.1f%% set %.2f | %-16s top %5.1f%% set %.2f | difference %+.2f (%+.2f to %+.2f), %d families' % (
                        b, 100 * target, ' '.join(kinds[i][1:]), 100 * a_top, a_size, ' '.join(kinds[j][1:]), 100 * b_top, b_size, d.mean(), lo, hi, len(common)))

# ------------------------------------------------------------------ 9. confident calls
print()
print('9. CONFIDENT CALLS: what the calibrated set does with them (%s score, one family out)' % FINAL)
w9 = writer('final_confident_calls.csv', ['stated_probability_of_top_candidate', 'top_candidate', 'target', 'records', 'families', 'actual_in_set', 'single_structure_sets', 'mean_set_size', 'mean_candidates'])
for lo, name in [(0.99, '99% or more'), (0.999, '99.9% or more')]:
    for right in (False, True):
        for target in (0.90, 0.80):
            cov, size = RES[(target, FINAL)]
            idx = [i for i, r in enumerate(recs) if r['ptop'] >= lo and r['first'] == right]
            n = len(idx)
            w9.writerow([name, 'right' if right else 'wrong', target, n, len({recs[i]['family'] for i in idx}), sum(cov[i] for i in idx), sum(size[i] == 1 for i in idx),
                         '%.8f' % (sum(size[i] for i in idx) / n), '%.8f' % (sum(recs[i]['K'] for i in idx) / n)])
            print('   stated %-13s and %-5s: %4d records (%3d families) | %d%% set holds the actual structure in %4d (%5.1f%%) | one structure in %4d (%5.1f%%) | mean size %.2f of %.2f' % (
                name, 'right' if right else 'wrong', n, len({recs[i]['family'] for i in idx}), 100 * target, sum(cov[i] for i in idx), 100 * sum(cov[i] for i in idx) / n,
                sum(size[i] == 1 for i in idx), 100 * sum(size[i] == 1 for i in idx) / n, sum(size[i] for i in idx) / n, sum(recs[i]['K'] for i in idx) / n))
for target in (0.90, 0.80):
    cov, size = RES[(target, FINAL)]
    one = [i for i in range(len(recs)) if size[i] == 1]
    print('   %d%% sets with one structure: %d of %d records (%.1f%%), right in %d (%.1f%%) | stated 99%% or more: %d records (%.1f%%), right in %.1f%%' % (
        100 * target, len(one), len(recs), 100 * len(one) / len(recs), sum(recs[i]['first'] for i in one), 100 * sum(recs[i]['first'] for i in one) / len(one),
        len(c99), 100 * len(c99) / len(recs), 100 * sum(r['first'] for r in c99) / len(c99)))

# ------------------------------------------------------------------ 10. energies of the conformers
print()
print('10. UNCERTAIN CONFORMER ENERGIES: likelihood averaged over a noise of 4 kJ/mol on the conformer energies against the plain likelihood')
rec_e = cs.load(table, 'log_ratio_energy_4')
keys_e = {(r['wf'], r['case']) for r in rec_e}
rec_p = [r for r in recs if (r['wf'], r['case']) in keys_e]
assert [(r['wf'], r['case']) for r in rec_p] == [(r['wf'], r['case']) for r in rec_e]
print('   records with conformer data %d | families %d | workflows %d | top candidate right %d -> %d' % (
    len(rec_e), len({r['family'] for r in rec_e}), len(workflows(rec_e)), sum(r['first'] for r in rec_p), sum(r['first'] for r in rec_e)))
w10 = writer('final_energy.csv', ['target', 'benchmark', 'records', 'families', 'coverage_plain', 'coverage_averaged', 'size_plain', 'size_averaged', 'difference', 'lo', 'hi',
                                  'single_plain', 'single_averaged', 'families_smaller', 'families_larger'])
for target in (0.90, 0.80):
    sp = sets(rec_p, FINAL, lofo_thresholds(rec_p, FINAL, target))
    se = sets(rec_e, FINAL, lofo_thresholds(rec_e, FINAL, target))
    for b in [x for x in BENCH if any(r['bench'] == x for r in rec_e)] + ['all']:
        idx = [i for i, r in enumerate(rec_e) if b == 'all' or r['bench'] == b]
        if b != 'all' and not bench_reaches(rec_e, b, target):
            continue
        xp = cs.summarize([rec_p[i] for i in idx], [sp[0][i] for i in idx], [sp[1][i] for i in idx])
        xe = cs.summarize([rec_e[i] for i in idx], [se[0][i] for i in idx], [se[1][i] for i in idx])
        d = paired([rec_e[i] for i in idx], [se[1][i] for i in idx], [sp[1][i] for i in idx], seed=SEED + 10)
        um = cs.unit_means([rec_e[i] for i in idx], [(float(se[1][i]) - float(sp[1][i]),) for i in idx])
        fam = collections.defaultdict(list)
        for (wf, f), v in um.items():
            fam[f].append(v[0])
        fd = [sum(v) / len(v) for v in fam.values()]
        w10.writerow([target, b, len(idx), len(fam), '%.8f' % xp['coverage'], '%.8f' % xe['coverage'], '%.8f' % xp['size'], '%.8f' % xe['size'], '%.8f' % d[0], '%.8f' % d[1], '%.8f' % d[2],
                      '%.8f' % xp['single'], '%.8f' % xe['single'], sum(1 for x in fd if x < -1e-9), sum(1 for x in fd if x > 1e-9)])
        print('   %d%% %-8s records %4d families %2d | coverage %5.1f%% -> %5.1f%% | set %.2f -> %.2f, difference %+.2f (%+.2f to %+.2f) | single %2.0f%% -> %2.0f%% | families smaller %d, larger %d' % (
            100 * target, b, len(idx), len(fam), 100 * xp['coverage'], 100 * xe['coverage'], xp['size'], xe['size'], d[0], d[1], d[2], 100 * xp['single'], 100 * xe['single'],
            sum(1 for x in fd if x < -1e-9), sum(1 for x in fd if x > 1e-9)))

# ------------------------------------------------------------------ 11. every record
w11 = writer('final_records.csv', ['benchmark', 'method', 'geometry', 'level', 'case', 'family', 'candidates', 'shifts', 'top_candidate', 'actual_structure', 'top_right',
                                  'stated_probability_top', 'stated_probability_actual', 'distance_of_actual_per_shift'] +
             [x for t in TARGETS for x in ('cutoff_%d' % (100 * t), 'set_size_%d' % (100 * t), 'covered_%d' % (100 * t))])
THR = {t: lofo_thresholds(recs, FINAL, t) for t in TARGETS}
for i, r in enumerate(recs):
    line = [r['bench'], r['method'], r['geo'], r['level'], r['case'], r['family'], r['K'], r['N'], r['names'][r['top']], r['names'][r['k']], int(r['first']), '%.6g' % r['ptop'],
            '%.6g' % r['pa'], '%.5f' % r['t'][FINAL]]
    for t in TARGETS:
        line += ['' if THR[t][i] == INF else '%.5f' % THR[t][i], RES[(t, FINAL)][1][i], int(RES[(t, FINAL)][0][i])]
    w11.writerow(line)

# ------------------------------------------------------------------ 12. random halves of the families
print()
print('12. RANDOM HALVES of the families (%d splits): coverage (5%% and 95%% of the splits), mean set size, sets with every candidate' % SPLITS)
w12 = writer('final_random_halves.csv', ['target', 'score', 'benchmark', 'coverage', 'coverage_5', 'coverage_95', 'mean_set_size', 'sets_with_every_candidate'])
fams = sorted({r['family'] for r in recs})
for target in (0.90, 0.80):
    rg = np.random.default_rng(SEED)
    agg = collections.defaultdict(list)
    for _ in range(SPLITS):
        perm = rg.permutation(len(fams))
        cal = {fams[i] for i in perm[:len(fams) // 2]}
        for score in (FINAL, 'likelihood ratio'):
            cov, size, test = [], [], []
            for wf, rr in by_wf.items():
                thr = cs.threshold([(r['t'][score], r['family']) for r in rr if r['family'] in cal], target)
                for r in rr:
                    if r['family'] not in cal:
                        c, s = cs.apply_set(r, score, thr)
                        cov.append(c)
                        size.append(s)
                        test.append(r)
            for b in BENCH:
                idx = [i for i, r in enumerate(test) if r['bench'] == b]
                x = cs.summarize([test[i] for i in idx], [cov[i] for i in idx], [size[i] for i in idx])
                agg[(score, b)].append((x['coverage'], x['size'], x['every']))
    print('target %d%%' % (100 * target))
    for score in (FINAL, 'likelihood ratio'):
        cells = []
        for b in BENCH:
            a = np.array(agg[(score, b)])
            w12.writerow([target, score, b, '%.8f' % a[:, 0].mean(), '%.8f' % np.quantile(a[:, 0], 0.05), '%.8f' % np.quantile(a[:, 0], 0.95), '%.8f' % a[:, 1].mean(), '%.8f' % a[:, 2].mean()])
            cells.append('%5.1f%% (%2.0f-%3.0f) %5.2f, all %2.0f%%' % (100 * a[:, 0].mean(), 100 * np.quantile(a[:, 0], 0.05), 100 * np.quantile(a[:, 0], 0.95), a[:, 1].mean(), 100 * a[:, 2].mean()))
        print('   %-27s ' % score + ' '.join('%-8s %s' % (b, c) for b, c in zip(BENCH, cells)))
