"""Which candidate is the reported ("actual") structure of each experiment of the 2017 deposit?

The supporting information of the paper (Org. Biomol. Chem. 2017, 15, 8998; section S-4) lists, per compound and for
the workflow "shifts B3LYP/6-31G**, energies B3LYP/6-31G**", six numbers of the actual structure that do not depend on
any probability model: Pearson correlation of calculated and measured shifts (13C, 1H), mean absolute error before and
after the linear correction (13C, 1H). The same six numbers are computed here for every candidate of every experiment
from the raw files; the candidate that shows the printed numbers is the actual structure.
Usage: truth17.py <si table.csv> <dep17_shifts.csv> <dep17_candidates.csv> <level> <out.csv>
"""
import collections
import csv
import math
import sys


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sxx * syy)


def stats(rows):
    out = {}
    for el, key in (('C', 'c'), ('H', 'h')):
        r = [x for x in rows if x['nucleus'] == el]
        u = [float(x['shift_unscaled']) for x in r]
        s = [float(x['shift_scaled']) for x in r]
        e = [float(x['exp']) for x in r]
        out[key + '_corr'] = pearson(u, e)
        out[key + '_mae'] = sum(abs(a - b) for a, b in zip(u, e)) / len(e)
        out[key + '_cmae'] = sum(abs(a - b) for a, b in zip(s, e)) / len(e)
    return out


KEYS = ['c_corr', 'h_corr', 'c_mae', 'h_mae', 'c_cmae', 'h_cmae']
TOL = {'c_corr': 0.000101, 'h_corr': 0.000101, 'c_mae': 0.0101, 'h_mae': 0.0101, 'c_cmae': 0.0101, 'h_cmae': 0.0101}


def main(si, shifts, candidates, level, out):
    SI = list(csv.DictReader(open(si)))
    rows = collections.defaultdict(list)
    for r in csv.DictReader(open(shifts)):
        if r['level'] == level:
            rows[(r['experiment'], int(r['isomer']))].append(r)
    cand = {}
    for r in csv.DictReader(open(candidates)):
        if r['level'] == level:
            cand[(r['experiment'], int(r['isomer']))] = r
    st = {k: stats(v) for k, v in rows.items()}
    exps = sorted({k[0] for k in st})
    used = {}
    res = []
    for s in SI:
        target = {k: float(s[k]) for k in KEYS}
        hits = [k for k, v in st.items() if all(abs(v[q] - target[q]) <= TOL[q] for q in KEYS)]
        # nearest candidate in units of the printed precision, for the rows without an exact hit
        dist = sorted((max(abs(v[q] - target[q]) / TOL[q] for q in KEYS), k) for k, v in st.items())
        res.append((s, hits, dist[:2]))
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['si_compound', 'experiment', 'actual_isomer', 'structure', 'candidates', 'candidates_showing_the_printed_numbers',
                'si_dp4_percent', 'my_dp4_percent', 'si_rank', 'my_rank', 'note'])
    n_unique = n_none = n_multi = 0
    agree_rank = agree_prob = agree_rank_r = 0
    for s, hits, dist in res:
        if len(hits) == 1:
            n_unique += 1
            exp, k = hits[0]
            note = ''
        elif not hits:
            n_none += 1
            exp, k = dist[0][1]
            note = 'no candidate within the printed precision; nearest is %.1f units away, the next %.1f' % (dist[0][0], dist[1][0])
        else:
            n_multi += 1
            exp, k = hits[0]
            note = 'several candidates fit: ' + '; '.join('%s isomer %d' % h for h in hits)
        n = max(q[1] for q in cand if q[0] == exp)
        p = [float(cand[(exp, i)]['p_all']) for i in range(1, n + 1)]
        rank = 1 + sum(1 for x in p if x > p[k - 1])
        rank_r = 1 + sum(1 for x in p if round(100 * x, 1) > round(100 * p[k - 1], 1))
        agree_rank_r += rank_r == int(s['rank'])
        used.setdefault(exp, []).append(s['compound'])
        w.writerow([s['compound'], exp, k, cand[(exp, k)]['structure'], n, len(hits), s['dp4_prob'], '%.1f' % (100 * p[k - 1]),
                    s['rank'], rank, note])
        agree_rank += rank == int(s['rank'])
        agree_prob += abs(100 * p[k - 1] - float(s['dp4_prob'])) <= 0.1001
        print('%-32s -> %-22s isomer %2d of %2d (%-14s) hits %d | DP4 SI %5s mine %5.1f | rank SI %s mine %d %s'
              % (s['compound'], exp, k, n, cand[(exp, k)]['structure'], len(hits), s['dp4_prob'], 100 * p[k - 1], s['rank'], rank,
                 ('| ' + note) if note else ''))
    print()
    print('rows of the SI table: %d; exactly one candidate shows the printed numbers: %d; none: %d; several: %d' % (len(SI), n_unique, n_none, n_multi))
    print('rank of the actual structure equal to the SI: %d (%d when ranks are taken from the probabilities rounded to 0.1 point); DP4 probability within 0.1 point of the SI: %d' % (agree_rank, agree_rank_r, agree_prob))
    print('experiments without a row in the SI table:', ', '.join(e for e in exps if e not in used) or 'none')
    print('experiments matched by more than one row:', ', '.join(e for e, v in used.items() if len(v) > 1) or 'none')


if __name__ == '__main__':
    main(*sys.argv[1:6])
