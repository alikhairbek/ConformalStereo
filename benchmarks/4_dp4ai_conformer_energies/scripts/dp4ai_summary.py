"""One line per DP4-AI record (molecule x protocol) from the table of dp4ai_energy.py, and what the energy-averaged
likelihood changes: how often the reported structure is first, how often a call of 99% or more is right.
Usage: dp4ai_summary.py <dp4ai_energy_marginal.csv> <dp4ai_truth.csv> <out.csv>
"""
import collections
import csv
import math
import sys


def softmax(v):
    m = max(v)
    w = [math.exp(x - m) for x in v]
    z = sum(w)
    return [x / z for x in w]


def main(table, truth, out):
    T = {r['experiment']: int(r['actual_isomer']) for r in csv.DictReader(open(truth))}
    g = collections.OrderedDict()
    for r in csv.DictReader(open(table)):
        g.setdefault((r['experiment'], r['shielding_level']), []).append(r)
    cols = [('2020_rule', 'logL_2020_rule'), ('0', 'logL_noise_0'), ('2', 'logL_noise_2'), ('4', 'logL_noise_4'), ('8', 'logL_noise_8'), ('16', 'logL_noise_16')]
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['case', 'family', 'protocol', 'geometry', 'energies', 'candidates', 'n_C', 'n_H', 'exact_counts', 'qc_clean', 'actual_isomer']
               + [x for name, _ in cols for x in ('p_actual_%s' % name, 'rank_actual_%s' % name, 'p_top_%s' % name)])
    recs = []
    for (case, prot), rr in g.items():
        rr = sorted(rr, key=lambda r: int(r['isomer']))
        k = T[case] - 1
        row = [case, rr[0]['archive'], prot, rr[0]['geometry'], rr[0]['energy_level'], len(rr), rr[0]['n_C'], rr[0]['n_H'], rr[0]['exact_counts'],
               rr[0]['qc_clean'], T[case]]
        rec = dict(case=case, clean=rr[0]['exact_counts'] == '1' and rr[0]['qc_clean'] == '1', exact=rr[0]['exact_counts'] == '1')
        for name, col in cols:
            p = softmax([float(r[col]) for r in rr])
            rank = 1 + sum(1 for x in p if x > p[k])
            row += ['%.6g' % p[k], rank, '%.6g' % max(p)]
            rec[name] = (p[k], rank, max(p))
        w.writerow(row)
        recs.append(rec)
    for label, keep in (('all records', lambda r: True), ('clean records with one value per atom', lambda r: r['clean']),
                        ('records with fewer or more values than atoms', lambda r: not r['exact'])):
        rr = [r for r in recs if keep(r)]
        print('%s (%d)' % (label, len(rr)))
        print('   %-28s %14s %22s %14s' % ('likelihood', 'reported first', 'calls >=99%: right/all', 'reported <=1%'))
        for name, _ in cols:
            first = sum(r[name][1] == 1 for r in rr)
            conf = [r for r in rr if r[name][2] >= 0.99]
            right = sum(r[name][1] == 1 for r in conf)
            lab = 'sorted from the top (2020)' if name == '2020_rule' else ('optimal pairing, no noise' if name == '0' else 'noise %s kJ/mol' % name)
            print('   %-28s %5d (%4.1f%%) %12d / %-3d (%5.1f%%) %8d' % (lab, first, 100.0 * first / len(rr), right, len(conf),
                                                                   100.0 * right / max(1, len(conf)), sum(r[name][0] <= 0.01 for r in rr)))


if __name__ == '__main__':
    main(*sys.argv[1:4])
