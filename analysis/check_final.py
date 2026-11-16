"""Checks of the final-method analysis that do not reuse its evaluation code.
Usage: check_final.py <pooled_candidates.csv> <folder with the final tables>
 1. the cut-off, the set and the coverage of every record recomputed from the pooled table with plain Python
 2. the numbers per benchmark recomputed from final_records.csv
 3. coverage when the actual structure is drawn at random: one family out must still reach the target (it holds by construction)
 4. who is missed: how the misses are spread over the families
"""
import sys, csv, math, collections, random
table, out = sys.argv[1:3]
TARGETS = (0.95, 0.90, 0.80)
raw = collections.OrderedDict()
for r in csv.DictReader(open(table)):
    raw.setdefault((r['benchmark'], r['method'], r['geometry'], r['level'], r['case']), []).append(r)
recs = []
for key, rr in raw.items():
    ell = [float(r['log_ratio']) for r in rr]
    n = int(rr[0]['n_nuclei'])
    k = [int(r['is_true']) for r in rr].index(1)
    s = [(max(ell) - x) / n for x in ell]
    recs.append(dict(key=key, wf=key[:4], fam=rr[0]['family'], s=s, k=k, K=len(rr)))
by_wf = collections.defaultdict(list)
for r in recs:
    by_wf[r['wf']].append(r)


def cutoff(cal, target, truth=lambda r: r['k']):
    """cal: records; family-weighted quantile of the scores of the actual structures"""
    n_f = collections.Counter(r['fam'] for r in cal)
    G = len(n_f)
    pts = sorted((r['s'][truth(r)], 1.0 / n_f[r['fam']] / (G + 1)) for r in cal)
    acc = 0.0
    for s, w in pts:
        acc += w
        if acc >= target - 1e-12:
            return s
    return float('inf')


def lofo(target, truth=lambda r: r['k']):
    res = {}
    for wf, rr in by_wf.items():
        for f in {r['fam'] for r in rr}:
            c = cutoff([r for r in rr if r['fam'] != f], target, truth)
            for r in rr:
                if r['fam'] == f:
                    res[r['key']] = (c, sum(1 for x in r['s'] if x <= c), r['s'][truth(r)] <= c)
    return res


def family_mean(pairs):
    """pairs: [((workflow, family), value)] -> mean over units of the unit means"""
    d = collections.defaultdict(list)
    for u, v in pairs:
        d[u].append(v)
    return sum(sum(v) / len(v) for v in d.values()) / len(d)


ok = True
tab = {tuple(r[c] for c in ('benchmark', 'method', 'geometry', 'level', 'case')): r for r in csv.DictReader(open(out + '/final_records.csv'))}
print('1. every record recomputed')
for target in TARGETS:
    res = lofo(target)
    t = str(int(100 * target))
    bad = 0
    for r in recs:
        c, size, cov = res[r['key']]
        row = tab[r['key']]
        if row['cutoff_' + t] == '':
            if c != float('inf'):
                bad += 1
            continue
        if abs(float(row['cutoff_' + t]) - c) > 1e-4 or int(row['set_size_' + t]) != size or int(row['covered_' + t]) != int(cov):
            bad += 1
    print('   %s%%: records that differ from final_records.csv: %d of %d' % (t, bad, len(recs)))
    ok &= bad == 0
print('2. per benchmark from the recomputed sets (coverage, mean set size) against final_scores.csv')
sc = {(r['target'], r['benchmark']): r for r in csv.DictReader(open(out + '/final_scores.csv')) if r['method'] == 'likelihood ratio per shift'}
for target in TARGETS:
    res = lofo(target)
    for b in ['JOC2015', 'STX150', 'DP4-AI', 'Cam2017', 'Cam2019']:
        row = sc.get((str(target), b))
        if row is None:
            continue
        rr = [r for r in recs if r['key'][0] == b]
        cov = family_mean([((r['wf'], r['fam']), float(res[r['key']][2])) for r in rr])
        size = family_mean([((r['wf'], r['fam']), float(res[r['key']][1])) for r in rr])
        same = abs(cov - float(row['coverage'])) < 1e-4 and abs(size - float(row['mean_set_size'])) < 1e-3
        ok &= same
        print('   %d%% %-8s coverage %.4f (table %s) size %.3f (table %s) %s' % (100 * target, b, cov, row['coverage'], size, row['mean_set_size'], 'same' if same else 'DIFFERENT'))
print('3. the actual structure drawn at random (20 draws): one family out, coverage per benchmark (lowest draw, mean) and mean share of candidates kept')
rng = random.Random(7)
for target in (0.90, 0.80):
    agg = collections.defaultdict(list)
    for _ in range(20):
        pick = {r['key']: rng.randrange(r['K']) for r in recs}
        res = lofo(target, truth=lambda r: pick[r['key']])
        for b in ['JOC2015', 'STX150', 'DP4-AI', 'Cam2017', 'Cam2019']:
            rr = [r for r in recs if r['key'][0] == b]
            agg[b].append((family_mean([((r['wf'], r['fam']), float(res[r['key']][2])) for r in rr]), family_mean([((r['wf'], r['fam']), res[r['key']][1] / r['K']) for r in rr])))
    print('   %d%%: ' % (100 * target) + ' | '.join('%s %.1f%% (lowest %.1f%%), kept %.0f%%' % (b, 100 * sum(x[0] for x in v) / len(v), 100 * min(x[0] for x in v), 100 * sum(x[1] for x in v) / len(v))
                                              for b, v in agg.items()))
    ok &= all(min(x[0] for x in v) >= target - 0.02 for v in agg.values())
print('4. who is missed (sets of the %s score, one family out)' % 'likelihood ratio per shift')
w = csv.writer(open(out + '/final_misses.csv', 'w', newline=''))
w.writerow(['target', 'benchmark', 'families', 'units', 'units_missed', 'families_never_missed', 'families_missed_in_half_their_records_or_more', 'share_of_misses_in_the_three_families_with_most',
            'families_with_most_misses'])
for target in (0.90, 0.80):
    res = lofo(target)
    for b in ['JOC2015', 'STX150', 'DP4-AI', 'Cam2017', 'Cam2019']:
        rr = [r for r in recs if r['key'][0] == b]
        fam = collections.defaultdict(list)
        for r in rr:
            fam[r['fam']].append(0.0 if res[r['key']][2] else 1.0)
        miss = {f: sum(v) for f, v in fam.items()}
        rate = {f: sum(v) / len(v) for f, v in fam.items()}
        tot = sum(miss.values())
        top3 = sorted(miss, key=lambda f: -miss[f])[:3]
        w.writerow([target, b, len(fam), len(rr), int(tot), sum(1 for f in fam if miss[f] == 0), sum(1 for f in fam if rate[f] >= 0.5), '%.8f' % (sum(miss[f] for f in top3) / tot if tot else 0),
                    '; '.join('%s (%d of %d)' % (f, miss[f], len(fam[f])) for f in top3 if miss[f] > 0)])
        print('   %d%% %-8s families %2d | records missed %3d of %4d | families never missed %2d | missed in half their records or more %2d | three families carry %3.0f%% of the misses: %s' % (
            100 * target, b, len(fam), tot, len(rr), sum(1 for f in fam if miss[f] == 0), sum(1 for f in fam if rate[f] >= 0.5), 100 * sum(miss[f] for f in top3) / tot if tot else 0,
            '; '.join('%s %d/%d' % (f.split(':')[-1], miss[f], len(fam[f])) for f in top3 if miss[f] > 0)))
print('ALL CHECKS PASSED' if ok else 'A CHECK FAILED')
