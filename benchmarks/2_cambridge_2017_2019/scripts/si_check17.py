"""The four per-compound tables of the paper's SI (section S-4) against the raw deposit.

Each table gives, for the actual structure of each of the 42 compounds: DP4 probability, rank, Pearson correlation
(13C, 1H), mean absolute error before and after linear correction (13C, 1H) for one workflow:
   shifts B3LYP/6-31G**     energies B3LYP/6-31G**
   shifts mPW1PW91/6-311G*  energies mPW1PW91/6-311G* | B3LYP/6-31G** | M06-2X/6-31G**
The six model-free numbers are recomputed here from the Gaussian outputs with the shieldings of one level and the
Boltzmann weights of another. The DP4 probability is recomputed with the single-Gaussian model of the 2015 code.
Usage: si_check17.py <res/core> <experiments.json> <energy_source.json> <structures.jsonl.gz> <dep17_truth.csv> <in dir> <out.csv>
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core17
import dep17
import raw17
from truth17 import KEYS, TOL, pearson

TABLES = [('b3lyp', 'BF1', 'BF1', 'B3LYP/6-31G**', 'B3LYP/6-31G**'),
          ('mpw_mpw', 'BF4', 'BF4', 'mPW1PW91/6-311G*', 'mPW1PW91/6-311G*'),
          ('mpw_b3lyp', 'BF4', 'BF1', 'mPW1PW91/6-311G*', 'B3LYP/6-31G**'),
          ('mpw_m062x', 'BF4', 'BF7', 'mPW1PW91/6-311G*', 'M06-2X/6-31G**')]


def six(row):
    out = {}
    for el, key in (('C', 'c'), ('H', 'h')):
        u = [p[2] for p in row[el]]
        s = [p[3] for p in row[el]]
        e = [p[4] for p in row[el]]
        out[key + '_corr'] = pearson(u, e)
        out[key + '_mae'] = sum(abs(a - b) for a, b in zip(u, e)) / len(e)
        out[key + '_cmae'] = sum(abs(a - b) for a, b in zip(s, e)) / len(e)
    return out


def main(core, experiments, source, structures, truth, indir, out):
    S = dep17.Small(core)
    E = {e['name']: e for e in json.load(open(experiments))}
    const = core17.reference_constants(json.load(open(source)), core17.solvents_of(structures))
    T = {r['si_compound']: r for r in csv.DictReader(open(truth))}
    data = raw17.load(core, {e['archive'] for e in E.values()})
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['table', 'shifts', 'energies', 'si_compound', 'experiment', 'actual_isomer'] + ['si_' + k for k in KEYS] + ['my_' + k for k in KEYS]
               + ['worst_deviation_in_printed_units', 'si_dp4_percent', 'my_dp4_percent_single_gaussian', 'si_rank', 'my_rank'])
    for key, sl, el, sname, ename in TABLES:
        SI = list(csv.DictReader(open(os.path.join(indir, 'obc2017_si_s4_%s.csv' % key))))
        exact = near = prob = rank = top = 0
        bad = []
        for s in SI:
            t = T[s['compound']]
            e = E[t['experiment']]
            k = int(t['actual_isomer'])
            desc, rows = core17.candidates(S, data, e, sl, el, const)
            mine = six(rows[k - 1])
            dev = max(abs(mine[q] - float(s[q])) / (TOL[q] - 1e-6) for q in KEYS)
            p = core17.normalise([r['log']['C'] + r['log']['H'] for r in rows])
            rk = 1 + sum(1 for x in p if x > p[k - 1])
            exact += dev <= 1.0
            near += dev <= 2.0
            prob += abs(100 * p[k - 1] - float(s['dp4_prob'])) <= 0.1001
            rank += rk == int(s['rank'])
            top += rk == 1
            if dev > 2.0:
                bad.append('%s (%.0f)' % (s['compound'], dev))
            w.writerow([key, sname, ename, s['compound'], t['experiment'], k] + [s[q] for q in KEYS]
                       + ['%.4f' % mine[q] if 'corr' in q else '%.2f' % mine[q] for q in KEYS]
                       + ['%.1f' % dev, s['dp4_prob'], '%.1f' % (100 * p[k - 1]), s['rank'], rk])
        n = len(SI)
        print('shifts %-17s energies %-17s | six numbers within the printed precision %2d/%d, within twice that %2d/%d'
              ' | single-Gaussian DP4 within 0.1 point %2d/%d, same rank %2d/%d | actual structure first: SI %d, mine %d'
              % (sname, ename, exact, n, near, n, prob, n, rank, n, sum(1 for s in SI if s['rank'] == '1'), top))
        if bad:
            print('      further off:', ', '.join(bad))


if __name__ == '__main__':
    main(*sys.argv[1:8])
