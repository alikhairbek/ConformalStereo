"""One table for all five benchmarks: every candidate of every record with its log-likelihood relative to the best
candidate of the record, the method, the kind of geometry and the family (molecules that share candidates or that
occur in two benchmarks get one family, so that they are never split between calibration and test).

  JOC2015   48 cases x 24 NMR levels, DP4+ on B3LYP geometries (recomputed from the printed shieldings)
  STX150    the usable records of the 110 validation cases at six levels (MM-DP4+, xTB-DP4+, DP4+)
  DP4-AI    46 molecules x 4 protocols, DP4 (single Gaussian), unassigned shifts
  Cam2017   42 compounds x 12 NMR levels on MMFF geometries, DP4
  Cam2019   25 of the same experiments x 30 workflows on six kinds of geometry, DP4
Usage: pool.py <joc_loglik.csv> <stx_candidates.csv with log ratios> <dp4ai tables dir> <tables17 dir> <tables19 dir> <out dir>
"""
import collections
import csv
import os
import re
import sys


def main(joc, stx, ai, t17, t19, out):
    rows = []          # benchmark, method, geometry, level, case, family, clean, candidate, is_true, n, logL, logL_e4

    def add(bench, method, geo, level, case, family, clean, items):
        """items: [(candidate, is_true, n_nuclei, logL, logL_e4 or None)]"""
        m = max(x[3] for x in items)
        m4 = max(x[4] for x in items) if all(x[4] is not None for x in items) else None
        for cand, true, n, L, L4 in items:
            rows.append([bench, method, geo, level, case, family, clean, cand, true, n, '%.5f' % (L - m), '' if m4 is None else '%.5f' % (L4 - m4)])

    g = collections.OrderedDict()
    for r in csv.DictReader(open(joc)):
        g.setdefault((r['case'], r['level']), []).append(r)
    for (case, level), rr in g.items():
        add('JOC2015', 'DP4+', 'B3LYP', level, case, 'dp:' + rr[0]['family'], 1,
            [(r['candidate'], int(r['is_correct']), int(r['n_C']) + int(r['n_H']), float(r['logL_scaled']) + float(r['logL_unscaled']), None) for r in rr])
    g = collections.OrderedDict()
    for r in csv.DictReader(open(stx)):
        if r['status'] in ('verified', 'close', 'corrected'):
            g.setdefault((r['case'], r['level']), []).append(r)
    kind = {'MMFF': ('MM-DP4+', 'MMFF'), 'xTB': ('xTB-DP4+', 'xTB'), 'B3L': ('DP4+', 'B3LYP')}
    for (case, level), rr in g.items():
        method, geo = kind[level.split('_')[0]]
        add('STX150', method, geo, level, case, 'dp:' + re.match(r'\d+', case).group(), int(rr[0]['status'] != 'close'),
            [(r['candidate'], int(r['is_correct']), int(r['n_nuclei']), float(r['log_ratio']), None) for r in rr])
    T = {r['experiment']: int(r['actual_isomer']) for r in csv.DictReader(open(os.path.join(ai, 'dp4ai_truth.csv')))}
    g = collections.OrderedDict()
    for r in csv.DictReader(open(os.path.join(ai, 'dp4ai_energy_marginal.csv'))):
        g.setdefault((r['experiment'], r['shielding_level']), []).append(r)
    for (case, level), rr in g.items():
        add('DP4-AI', 'DP4', 'MMFF' if rr[0]['geometry'] == 'MMFF' else 'B3LYP', level, case, 'ai:' + rr[0]['archive'],
            int(rr[0]['exact_counts'] == '1' and rr[0]['qc_clean'] == '1'),
            [(r['isomer'], int(int(r['isomer']) == T[case]), int(r['n_C']) + int(r['n_H']), float(r['logL_noise_0']), float(r['logL_noise_4'])) for r in rr])
    T = {r['experiment']: int(r['actual_isomer']) for r in csv.DictReader(open(os.path.join(t17, 'dep17_truth.csv')))}
    E4 = {(r['experiment'], r['shielding_level'], r['isomer']): float(r['logL_noise_4']) for r in csv.DictReader(open(os.path.join(t17, 'dep17_energy_marginal.csv')))}
    fam17 = {}
    g = collections.OrderedDict()
    for r in csv.DictReader(open(os.path.join(t17, 'dep17_loglik.csv'))):
        if r['experiment'] in T:
            g.setdefault((r['experiment'], r['shielding_level']), []).append(r)
            fam17[r['experiment']] = 'cam:' + r['archive']
    for (case, level), rr in g.items():
        add('Cam2017', 'DP4', 'MMFF', level, case, fam17[case], 1,
            [(r['isomer'], int(int(r['isomer']) == T[case]), int(r['n_C']) + int(r['n_H']), float(r['logL_C']) + float(r['logL_H']),
              E4.get((case, level, r['isomer']))) for r in rr])
    T = {r['experiment']: int(r['actual_isomer']) for r in csv.DictReader(open(os.path.join(t19, 'dep19_truth.csv')))}
    E4 = {(r['experiment'], r['shielding_level'], r['isomer']): float(r['logL_noise_4']) for r in csv.DictReader(open(os.path.join(t19, 'dep19_energy_marginal.csv')))}
    short = {'MMFF': 'MMFF', 'PM6': 'PM6', 'PM7': 'PM7', 'HF/6-31G(d,p)': 'HF', 'B3LYP/6-31G(d,p)': 'B3LYP', 'M06-2X/6-31G(d,p)': 'M06-2X'}
    g = collections.OrderedDict()
    for r in csv.DictReader(open(os.path.join(t19, 'dep19_loglik.csv'))):
        if r['experiment'] in T:
            g.setdefault((r['experiment'], r['shielding_level'], r['energy_level']), []).append(r)
    for (case, sl, el), rr in g.items():
        own = sl == el
        add('Cam2019', 'DP4', short[rr[0]['geometry']], sl if own else '%s+E(%s)' % (sl, el), case, fam17[case], 1,
            [(r['isomer'], int(int(r['isomer']) == T[case]), int(r['n_C']) + int(r['n_H']), float(r['logL_C']) + float(r['logL_H']),
              E4.get((case, sl, r['isomer'])) if own else None) for r in rr])
    w = csv.writer(open(os.path.join(out, 'pooled_candidates.csv'), 'w', newline=''))
    w.writerow(['benchmark', 'method', 'geometry', 'level', 'case', 'family', 'clean', 'candidate', 'is_true', 'n_nuclei', 'log_ratio', 'log_ratio_energy_4'])
    w.writerows(rows)
    rec = collections.OrderedDict()
    for r in rows:
        rec.setdefault(tuple(r[:5]), []).append(r)
    bad = [k for k, v in rec.items() if sum(x[8] for x in v) != 1]
    print('candidates %d, records %d, records without exactly one actual structure %d' % (len(rows), len(rec), len(bad)))
    by = collections.OrderedDict()
    for k, v in rec.items():
        b = by.setdefault(k[0], dict(records=0, cases=set(), families=set(), levels=set(), cands=0))
        b['records'] += 1
        b['cases'].add(k[4])
        b['families'].add(v[0][5])
        b['levels'].add(k[3])
        b['cands'] += len(v)
    fams = set()
    for name, b in by.items():
        print('%-8s records %4d  cases %3d  families %3d  levels %2d  candidates per record %.2f'
              % (name, b['records'], len(b['cases']), len(b['families']), len(b['levels']), b['cands'] / b['records']))
        fams |= b['families']
    print('families in all: %d' % len(fams))


if __name__ == '__main__':
    main(*sys.argv[1:7])
