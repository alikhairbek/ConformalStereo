"""Does averaging over errors in the conformer energies order the candidates better?
Reads the table written by energy_marginal.py and reports, per noise level: how often the actual structure is first,
the temperature that would make the probabilities honest, and the set sizes at matched coverage and under split
conformal prediction. Usage: energy_marginal_eval.py <energy_marginal.csv> <truth.csv> [splits] [seed] [by geometry: 0|1]
"""
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from error_models import evaluate
from scores import fitted_temperature


def main(table, truth, splits=300, seed=20261114, bygeo=0):
    splits, seed, bygeo = int(splits), int(seed), int(bygeo)
    T = {r['experiment']: int(r['actual_isomer']) for r in csv.DictReader(open(truth))}
    rows = list(csv.DictReader(open(table)))
    cols = [c for c in rows[0] if c.startswith('logL_noise_')]
    groups = collections.defaultdict(list)
    for r in rows:
        groups[(r['experiment'], r['shielding_level'])].append(r)
    geos = sorted({r['geometry'] for r in rows}) if bygeo else ['']
    for g in geos:
        if g:
            print(g)
        print('%-22s %-18s %6s | 90%% target: matched size (probability, per-nucleus log ratio), conformal coverage and size | 80%%: conformal size'
              % ('noise on the energies', 'actual first', 'T'))
        for c in cols:
            recs = []
            for (exp, level), rr in sorted(groups.items()):
                if g and rr[0]['geometry'] != g:
                    continue
                rr = sorted(rr, key=lambda r: int(r['isomer']))
                L = [float(r[c]) for r in rr]
                recs.append(dict(wf=level, family=rr[0]['archive'], k=T[exp] - 1, L=L, LC=L, LH=[0.0] * len(L),
                                 N=int(rr[0]['n_C']) + int(rr[0]['n_H'])))
            top = sum(1 for r in recs if max(range(len(r['L'])), key=lambda i: r['L'][i]) == r['k'])
            e90 = evaluate(recs, 0.90, splits, seed)
            e80 = evaluate(recs, 0.80, splits, seed)
            p, l = e90['probability'], e90['log ratio / sqrt(N)']
            print('%-22s %3d of %3d %5.1f%% %6.2f | %5.2f  %5.2f  | %5.1f%%  %5.2f  %5.2f | %5.2f  %5.2f'
                  % (c.replace('logL_noise_', '') + ' kJ/mol', top, len(recs), 100.0 * top / len(recs), fitted_temperature(recs),
                     p[0], l[0], 100 * p[2], p[3], l[3], e80['probability'][3], e80['log ratio / sqrt(N)'][3]))


if __name__ == '__main__':
    main(*sys.argv[1:6])
