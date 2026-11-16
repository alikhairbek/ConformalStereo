"""Benchmark table of the 2017 deposit and a first look at calibration.

One record = one compound of the paper (42) at one NMR level (12): the DP4 probability of every candidate,
recomputed from the raw files, and the candidate that is the actual structure (identified from the paper's SI).
Usage: benchmark17.py <dep17_candidates.csv> <dep17_truth.csv> <out dir> [splits] [seed]
"""
import collections
import csv
import math
import os
import random
import sys

LEVELS = ['BF%d' % i for i in range(1, 13)]
NAMES = {'BF1': 'B3LYP/6-31G(d,p)', 'BF2': 'B3LYP/6-311G(d)', 'BF3': 'mPW1PW91/6-31G(d,p)', 'BF4': 'mPW1PW91/6-311G(d)',
         'BF5': 'WP04/6-31G(d,p)', 'BF6': 'WP04/6-311G(d)', 'BF7': 'M06-2X/6-31G(d,p)', 'BF8': 'M06-2X/6-311G(d)',
         'BF9': 'B3LYP/gen', 'BF10': 'mPW1PW91/gen', 'BF11': 'WP04/gen', 'BF12': 'M06-2X/gen'}


def main(candidates, truth, out, splits=2000, seed=20261114):
    splits, seed = int(splits), int(seed)
    T = {r['experiment']: r for r in csv.DictReader(open(truth))}
    P = collections.defaultdict(dict)
    arch = {}
    for r in csv.DictReader(open(candidates)):
        if r['experiment'] in T:
            P[(r['experiment'], r['level'])][int(r['isomer'])] = {k: float(r[k]) for k in ('p_all', 'p_C', 'p_H')}
            arch[r['experiment']] = r['archive']
    rec = []
    w = csv.writer(open(os.path.join(out, 'dep17_benchmark.csv'), 'w', newline=''))
    w.writerow(['si_compound', 'experiment', 'family', 'level', 'level_name', 'candidates', 'actual_isomer',
                'p_actual', 'rank_actual', 'top_isomer', 'p_top', 'p_actual_C_only', 'p_actual_H_only'])
    for exp in sorted(T, key=lambda e: (arch[e], e)):
        k = int(T[exp]['actual_isomer'])
        for l in LEVELS:
            p = P[(exp, l)]
            n = len(p)
            pa = p[k]['p_all']
            rank = 1 + sum(1 for i in p if p[i]['p_all'] > pa)
            top = max(p, key=lambda i: p[i]['p_all'])
            w.writerow([T[exp]['si_compound'], exp, arch[exp], l, NAMES[l], n, k, '%.6g' % pa, rank, top, '%.6g' % p[top]['p_all'],
                        '%.6g' % p[k]['p_C'], '%.6g' % p[k]['p_H']])
            rec.append(dict(exp=exp, family=arch[exp], level=l, n=n, p=[p[i]['p_all'] for i in sorted(p)], k=k - 1, pa=pa, rank=rank,
                            ptop=p[top]['p_all'], pc=p[k]['p_C'], ph=p[k]['p_H'],
                            rank_c=1 + sum(1 for i in p if p[i]['p_C'] > p[k]['p_C']), rank_h=1 + sum(1 for i in p if p[i]['p_H'] > p[k]['p_H'])))
    fam = sorted({r['family'] for r in rec})
    print('compounds %d in %d families (molecule archives), levels %d, records %d; candidates per compound: %s'
          % (len(T), len(fam), len(LEVELS), len(rec),
             ', '.join('%d x %d' % (c, n) for n, c in sorted(collections.Counter(int(T[e]['candidates']) for e in T).items()))))
    print('chance of naming the actual structure by a random pick: %.1f%% on average'
          % (100 * sum(1.0 / int(T[e]['candidates']) for e in T) / len(T)))
    print()
    print('%-22s %9s %9s %9s %12s %14s %12s' % ('level', 'first', 'first(C)', 'first(H)', 'calls >=99%', 'of them right', 'actual <=1%'))
    tot = collections.Counter()
    for l in LEVELS + ['all']:
        rr = [r for r in rec if l == 'all' or r['level'] == l]
        first = sum(r['rank'] == 1 for r in rr)
        conf = [r for r in rr if r['ptop'] >= 0.99]
        right = sum(r['rank'] == 1 for r in conf)
        low = sum(r['pa'] <= 0.01 for r in rr)
        print('%-22s %4d/%-4d %4d/%-4d %4d/%-4d %8d %10d (%4.1f%%) %8d' % (NAMES.get(l, 'all levels'), first, len(rr), sum(r['rank_c'] == 1 for r in rr), len(rr),
              sum(r['rank_h'] == 1 for r in rr), len(rr), len(conf), right, 100.0 * right / max(1, len(conf)), low))
    # ---- sets built from the probabilities themselves: smallest set whose probabilities sum to the target
    print()
    for target in (0.90, 0.80):
        cov = size = 0
        for r in rec:
            order = sorted(range(r['n']), key=lambda i: -r['p'][i])
            s, chosen = 0.0, []
            for i in order:
                chosen.append(i)
                s += r['p'][i]
                if s >= target:
                    break
            cov += r['k'] in chosen
            size += len(chosen)
        print('sets taken from the DP4 probabilities (smallest set summing to %.0f%%): contain the actual structure in %d of %d records (%.1f%%), mean size %.2f'
              % (100 * target, cov, len(rec), 100.0 * cov / len(rec), size / len(rec)))
    # ---- split conformal, calibration and test separated by family, one threshold per level
    rng = random.Random(seed)
    by = collections.defaultdict(list)
    for r in rec:
        by[r['level']].append(r)
    for target in (0.90, 0.80):
        alpha = 1 - target
        cov_all, size_all, single_all, n_all = [], [], [], 0
        per_level = collections.defaultdict(list)
        for _ in range(splits):
            f = fam[:]
            rng.shuffle(f)
            cal = set(f[:len(f) // 2])
            c = s = sg = nt = 0
            for l in LEVELS:
                scores = sorted(1 - r['pa'] for r in by[l] if r['family'] in cal)
                idx = int(math.ceil((len(scores) + 1) * (1 - alpha)))
                q = scores[idx - 1] if idx <= len(scores) else 1.0
                cl = nl = 0
                for r in by[l]:
                    if r['family'] in cal:
                        continue
                    chosen = [i for i in range(r['n']) if 1 - r['p'][i] <= q + 1e-12]
                    hit = r['k'] in chosen
                    c += hit
                    cl += hit
                    s += len(chosen)
                    sg += len(chosen) == 1
                    nt += 1
                    nl += 1
                per_level[l].append(cl / nl)
            cov_all.append(c / nt)
            size_all.append(s / nt)
            single_all.append(sg / nt)
        m = sum(cov_all) / len(cov_all)
        lv = {l: sum(v) / len(v) for l, v in per_level.items()}
        print('split conformal, target %.0f%%: coverage %.1f%% (levels %.1f-%.1f%%), mean set size %.2f, single-structure sets %.1f%% (%d splits by family)'
              % (100 * target, 100 * m, 100 * min(lv.values()), 100 * max(lv.values()), sum(size_all) / len(size_all),
                 100 * sum(single_all) / len(single_all), splits))


if __name__ == '__main__':
    main(*sys.argv[1:6])
