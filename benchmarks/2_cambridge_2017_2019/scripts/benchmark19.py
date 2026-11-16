"""Benchmark table of the 2019 deposit and a first look at what better geometries and energies change.

One record = one experiment under one workflow (geometry, NMR level, energies). Probabilities from workflows19.py,
actual structures from truth19.py.
Usage: benchmark19.py <dep19_workflows.csv> <dep19_truth.csv> <out dir> [splits] [seed]
"""
import collections
import csv
import math
import os
import random
import sys

GEOS = ['MMFF', 'PM6', 'PM7', 'HF/6-31G(d,p)', 'B3LYP/6-31G(d,p)', 'M06-2X/6-31G(d,p)']
NMRS = ['B3LYP/6-31G(d,p)', 'mPW1PW91/6-311G(d)', 'M06-2X/6-31G(d,p)']


def line(label, rr):
    n = len(rr)
    first = sum(r['rank'] == 1 for r in rr)
    conf = [r for r in rr if r['ptop'] >= 0.99]
    right = sum(r['rank'] == 1 for r in conf)
    low = sum(r['pa'] <= 0.01 for r in rr)
    return '%-58s %3d  %3d (%5.1f%%)  %3d  %3d (%5.1f%%)  %3d' % (label, n, first, 100.0 * first / n, len(conf), right,
                                                             100.0 * right / max(1, len(conf)), low)


def main(workflows, truth, out, splits=2000, seed=20261114):
    splits, seed = int(splits), int(seed)
    T = {r['experiment']: r for r in csv.DictReader(open(truth))}
    P = collections.defaultdict(dict)
    meta = {}
    for r in csv.DictReader(open(workflows)):
        key = (r['experiment'], r['shielding_level'], r['energy_level'])
        P[key][int(r['isomer'])] = float(r['p_all'])
        meta[key] = r
    rec = []
    w = csv.writer(open(os.path.join(out, 'dep19_benchmark.csv'), 'w', newline=''))
    w.writerow(['si_compound', 'experiment', 'family', 'geometry', 'nmr', 'energies', 'shielding_level', 'energy_level', 'candidates',
                'actual_isomer', 'p_actual', 'rank_actual', 'top_isomer', 'p_top'])
    for key in sorted(P, key=lambda k: (meta[k]['archive'], k[0], GEOS.index(meta[k]['geometry']), NMRS.index(meta[k]['nmr']), meta[k]['energies'] != 'same level')):
        m = meta[key]
        if key[0] not in T:
            continue
        k = int(T[key[0]]['actual_isomer'])
        p = P[key]
        pa = p[k]
        rank = 1 + sum(1 for i in p if p[i] > pa)
        top = max(p, key=lambda i: p[i])
        w.writerow([T[key[0]]['si_compound'], key[0], m['archive'], m['geometry'], m['nmr'], m['energies'], key[1], key[2], len(p), k,
                    '%.6g' % pa, rank, top, '%.6g' % p[top]])
        rec.append(dict(exp=key[0], family=m['archive'], geo=m['geometry'], nmr=m['nmr'], en=m['energies'], wf=(key[1], key[2]),
                        n=len(p), p=[p[i] for i in sorted(p)], k=k - 1, pa=pa, rank=rank, ptop=p[top]))
    fam = sorted({r['family'] for r in rec})
    exps = sorted({r['exp'] for r in rec})
    print('experiments %d in %d families, workflows %d, records %d' % (len(exps), len(fam), len({r['wf'] for r in rec}), len(rec)))
    print()
    head = '%-58s %3s  %-13s  %3s  %-13s  %3s' % ('', 'n', 'actual first', '99%', 'of them right', '<=1%')
    print('by geometry (three NMR levels, energies of the NMR level)')
    print(head)
    for g in GEOS:
        print(line(g, [r for r in rec if r['geo'] == g and r['en'] == 'same level']))
    print()
    print('by geometry, NMR level and energies')
    print(head)
    for g in GEOS:
        for n in NMRS:
            for en in ('same level', 'M06-2X/6-31G(d,p)'):
                rr = [r for r in rec if r['geo'] == g and r['nmr'] == n and r['en'] == en]
                if rr:
                    print(line('%s // %s, energies: %s' % (n, g, en if en != 'same level' else 'NMR level'), rr))
    # ---- the same experiments on MMFF and on B3LYP geometries
    print()
    both = [e for e in exps if any(r['exp'] == e and r['geo'] == 'MMFF' for r in rec)]
    print('the %d experiments that have every kind of geometry, mPW1PW91 shifts with M06-2X energies:' % len(both))
    print(head)
    for g in GEOS:
        rr = [r for r in rec if r['exp'] in both and r['geo'] == g and r['nmr'] == NMRS[1] and r['en'] != 'same level']
        if rr:
            print(line(g, rr))
    # ---- sets
    print()
    for target in (0.90,):
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
    rng = random.Random(seed)
    by = collections.defaultdict(list)
    for r in rec:
        by[r['wf']].append(r)
    for target in (0.90, 0.80):
        alpha = 1 - target
        tot = collections.defaultdict(lambda: [0, 0, 0, 0])
        for _ in range(splits):
            f = fam[:]
            rng.shuffle(f)
            cal = set(f[:len(f) // 2])
            for wf, rr in by.items():
                scores = sorted(1 - r['pa'] for r in rr if r['family'] in cal)
                idx = int(math.ceil((len(scores) + 1) * (1 - alpha)))
                q = scores[idx - 1] if 0 < idx <= len(scores) else 1.0
                for r in rr:
                    if r['family'] in cal:
                        continue
                    chosen = [i for i in range(r['n']) if 1 - r['p'][i] <= q + 1e-12]
                    for g in (r['geo'], 'all'):
                        t = tot[g]
                        t[0] += r['k'] in chosen
                        t[1] += len(chosen)
                        t[2] += len(chosen) == 1
                        t[3] += 1
        print('split conformal, target %.0f%% (%d splits by family, one threshold per workflow):' % (100 * target, splits))
        for g in GEOS + ['all']:
            t = tot[g]
            print('   %-20s coverage %5.1f%%  mean set size %.2f  single-structure sets %5.1f%%' % (g, 100.0 * t[0] / t[3], t[1] / t[3], 100.0 * t[2] / t[3]))


if __name__ == '__main__':
    main(*sys.argv[1:6])
