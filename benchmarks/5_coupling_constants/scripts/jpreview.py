"""First look at the coupling constants of the 2017 deposit, from the depositors' own PyDP4 outputs in the _J folders.
Usage: jpreview.py <res/core of Kaggle run 2> <folder with parse17.py> <dep17_experiments.json> <dep17_truth.csv> <out folder> [dep17_records.csv]

Every such output prints, for each candidate, "DP4 based on C / H / J" and the probabilities with and without the
couplings; the numbers are the depositors'. Here each output is tied to one of the 42 compounds through its experimental
carbon shifts, its candidates are tied to the candidates of the shift-only output through their calculated carbon shifts
after scaling (one-to-one assignment, so a different order of the candidates does no harm; scaled shifts, because some
shift-only outputs were made at another level of theory), and the place of the actual structure is counted with and
without the couplings."""
import sys, os, re, csv, gzip, json, math, zipfile, collections
import numpy as np
from scipy.optimize import linear_sum_assignment
core, lib, exp_json, truth_csv, out = sys.argv[1:6]
records_csv = sys.argv[6] if len(sys.argv) > 6 else None
sys.path.insert(0, lib)
import parse17

rows = [l.rstrip('\n').split('\t') for l in gzip.open(os.path.join(core, 'small_index_dp4x2.tsv.gz'), 'rt', encoding='utf-8')][1:]
index = {r[0]: r for r in rows}
zf = zipfile.ZipFile(os.path.join(core, 'small_dp4x2_001.zip'))


def text_of(path):
    r = index[path]
    return zf.read((path if r[3] == '1' else r[4]).replace('!/', '/')).decode('latin-1')


truth = {r['experiment']: r for r in csv.DictReader(open(truth_csv, encoding='utf-8'))}
main = {}
for e in json.load(open(exp_json)):
    p = e['levels'].get('BF1')
    if p and p in index:
        iso = parse17.parse(text_of(p))['isomers']
        if iso:
            main[(e['archive'], e['name'])] = dict(iso=iso, exp=sorted(row[3] for row in iso[min(iso)]['C']), stems=e['stems'])


def cost(a, b, by_label):
    if by_label:
        da, db = {r[0]: r[2] for r in a}, {r[0]: r[2] for r in b}
        common = sorted(set(da) & set(db))
        if len(common) < max(len(da), len(db)):
            return float('inf')
        return math.sqrt(sum((da[x] - db[x]) ** 2 for x in common) / len(common))
    va, vb = sorted(r[2] for r in a), sorted(r[2] for r in b)
    if len(va) != len(vb):
        return float('inf')
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(va, vb)) / len(va))


TOL = 0.25                                        # points; the per-candidate terms are printed with three digits
jfiles = sorted(p for p in index if '_J/' in p and p.endswith('.dp4'))
recs, skipped = [], []
for path in jfiles:
    d = parse17.parse(text_of(path))
    iso = d['isomers']
    ks = sorted(iso)
    stem = os.path.basename(path)[:-4]
    if not iso or not all('J' in iso[k]['p'] for k in ks):
        skipped.append((path, 'no coupling term printed'))
        continue
    n_j = max(sum(len(row[5]) for row in iso[k]['H']) for k in ks)
    if n_j == 0:
        skipped.append((path, 'description without experimental couplings'))
        continue
    archive = path.split('BasisFunc/')[-1].split('.zip!/')[0]
    mine = sorted(row[3] for row in iso[ks[0]]['C'])
    hits = [key for key, m in main.items() if key[0] == archive and len(m['exp']) == len(mine) and max(abs(x - y) for x, y in zip(m['exp'], mine)) < 0.006]
    rec = dict(path=path, file=stem, archive=archive, n=len(ks), J_values=n_j, protons_with_J=max(sum(1 for row in iso[k]['H'] if row[5]) for k in ks))
    logp = {name: [math.log(iso[k]['p'][name]) if iso[k]['p'][name] > 0 else -math.inf for k in ks] for name in 'CHJ'}

    def post(names):
        s = [sum(logp[n][i] for n in names) for i in range(len(ks))]
        top = max(s)
        w = [math.exp(x - top) for x in s]
        return [100 * x / sum(w) for x in w]
    rec['post'] = {'CH': post('CH'), 'CHJ': post('CHJ'), 'J': post('J')}
    printed = {k: d['results'].get(k) for k in ('all', 'CH', 'J')}
    rec['agree'] = all(printed[a] and len(printed[a]) == len(ks) and max(abs(x - y) for x, y in zip(rec['post'][b], printed[a])) <= TOL
                       for a, b in (('all', 'CHJ'), ('CH', 'CH'), ('J', 'J')))
    if len(hits) == 1 and len(main[hits[0]]['iso']) == len(ks):
        m = main[hits[0]]
        rec['experiment'] = hits[0][1]
        qs = sorted(m['iso'])
        best = None
        for by_label in (True, False):
            c = np.array([[cost(iso[k]['C'], m['iso'][q]['C'], by_label) for q in qs] for k in ks])
            if not np.isfinite(c).all():
                continue
            ri, ci = linear_sum_assignment(c)
            if best is None or c[ri, ci].sum() < best[0]:
                best = (c[ri, ci].sum(), c, ci, by_label)
        if best is not None:
            total, c, ci, by_label = best
            rec['map'] = {ks[i]: qs[ci[i]] for i in range(len(ks))}
            rec['same_order'] = all(ks[i] == qs[ci[i]] for i in range(len(ks)))
            rec['match_rmsd'] = float(max(c[i, ci[i]] for i in range(len(ks))))
            rec['match_by'] = 'label' if by_label else 'sorted shifts'
            if hits[0][1] in truth:
                t_main = int(truth[hits[0][1]]['actual_isomer'])
                inv = {v: k for k, v in rec['map'].items()}
                t = inv[t_main]
                i = ks.index(t)
                rec['actual'] = t
                rec['actual_margin'] = float(min(c[i, j] for j in range(len(qs)) if j != ci[i]) - c[i, ci[i]]) if len(qs) > 1 else float('inf')
                rec['actual_rmsd'] = float(c[i, ci[i]])
                rec['compound'] = truth[hits[0][1]]['si_compound']
                for name in ('CH', 'CHJ', 'J'):
                    p = rec['post'][name]
                    rec['rank_' + name] = 1 + sum(1 for x in p if x > p[i])
                    rec['p_' + name] = p[i]
    else:
        rec['note'] = 'experiments with these carbon shifts in the archive: %d' % len(hits)
    recs.append(rec)

os.makedirs(out, exist_ok=True)
cols = ['compound', 'experiment', 'file', 'n', 'actual', 'same_order', 'match_by', 'match_rmsd', 'actual_rmsd', 'actual_margin', 'agree', 'protons_with_J', 'J_values',
        'rank_CH', 'rank_CHJ', 'rank_J', 'p_CH', 'p_CHJ', 'p_J', 'note', 'path']
with open(os.path.join(out, 'jpreview_outputs.csv'), 'w', newline='', encoding='utf-8') as fh:
    w = csv.writer(fh)
    w.writerow(cols)
    for r in recs:
        w.writerow([('%.4g' % r[c]) if isinstance(r.get(c), float) else r.get(c, '') for c in cols])
print('PyDP4 outputs in the _J folders: %d | with experimental couplings: %d | tied to one of the 42 compounds: %d' % (
    len(jfiles), len(recs), sum(1 for r in recs if 'actual' in r)))
print('printed probabilities reproduced from the per-candidate terms (within %.2f point): %d of %d' % (TOL, sum(r['agree'] for r in recs), len(recs)))
for r in recs:
    if 'actual' not in r:
        print('   not tied: %-28s %s' % (r['file'], r.get('note', 'candidates not matched')))
known = [r for r in recs if 'actual' in r]
print('candidates in the same order as in the shift-only output: %d of %d | matched by label %d, by sorted shifts %d' % (
    sum(r['same_order'] for r in known), len(known), sum(r['match_by'] == 'label' for r in known), sum(r['match_by'] != 'label' for r in known)))
print('match of the actual structure: largest distance %.3f ppm, smallest margin to the next candidate %.3f ppm' % (
    max(r['actual_rmsd'] for r in known), min(r['actual_margin'] for r in known)))
weak = [r for r in known if r['actual_margin'] < 0.05 or r['actual_rmsd'] > 1.0]
for r in weak:
    print('   weak match: %-26s distance %.3f, margin %.3f (%d candidates)' % (r['file'], r['actual_rmsd'], r['actual_margin'], r['n']))
by_exp = collections.OrderedDict()
pref = lambda r: (0 if r['file'].endswith('NMRj') else 1 if r['file'].endswith('jreal') else 2)
for r in known:
    if r in weak:
        continue
    cur = by_exp.get(r['experiment'])
    if cur is None or pref(r) < pref(cur):
        by_exp[r['experiment']] = r
use = list(by_exp.values())
n = len(use)
print('experiments used: %d of %d with a coupling output (%d outputs are second versions of the same experiment)' % (n, len(set(r['experiment'] for r in known)), len(known) - len(weak) - n))
summary = {}
for name, label in (('CH', 'shifts only (C and H)'), ('CHJ', 'shifts and couplings'), ('J', 'couplings only')):
    first = sum(1 for r in use if r['rank_' + name] == 1)
    hi = [r for r in use if max(r['post'][name]) >= 99]
    hi_ok = sum(1 for r in hi if r['rank_' + name] == 1)
    low = sum(1 for r in use if r['p_' + name] <= 1)
    summary[name] = dict(first=first, calls99=len(hi), calls99_right=hi_ok, actual_le1=low, mean_p=sum(r['p_' + name] for r in use) / n)
    print('  %-24s actual structure first: %2d of %d | calls at 99%% or more: %2d, right %2d | actual at 1%% or less: %2d | mean probability of the actual %.1f%%' % (
        label, first, n, len(hi), hi_ok, low, summary[name]['mean_p']))
better = [r for r in use if r['rank_CHJ'] < r['rank_CH']]
worse = [r for r in use if r['rank_CHJ'] > r['rank_CH']]
print('place of the actual structure once couplings are added: better in %d, worse in %d, same in %d' % (len(better), len(worse), n - len(better) - len(worse)))
miss = [r for r in use if r['rank_CHJ'] > 1]
for tag, group in (('better', better), ('worse', worse)):
    for r in group:
        print('   %-6s %-28s %2d candidates | place %d -> %d | probability %5.1f%% -> %5.1f%% | %d couplings on %d protons' % (
            tag, r['compound'][:28], r['n'], r['rank_CH'], r['rank_CHJ'], r['p_CH'], r['p_CHJ'], r['J_values'], r['protons_with_J']))
print('actual structure not first even with couplings (%d): %s' % (len(miss), '; '.join('%s place %d of %d' % (r['compound'], r['rank_CHJ'], r['n']) for r in miss)))
wrong99 = [r for r in use if max(r['post']['CHJ']) >= 99 and r['rank_CHJ'] > 1]
print('calls at 99%% or more that are wrong with couplings: %s' % '; '.join('%s (%.3f%% for the actual)' % (r['compound'], r['p_CHJ']) for r in wrong99))
jonly = [r for r in use if r['rank_J'] == 1 and r['rank_CH'] > 1]
print('couplings alone put the actual structure first where the shifts do not: %s' % '; '.join('%s (%d candidates)' % (r['compound'], r['n']) for r in jonly))
up = sum(1 for r in use if r['p_CHJ'] > r['p_CH'] + 0.05)
down = sum(1 for r in use if r['p_CHJ'] < r['p_CH'] - 0.05)
print('probability of the actual structure: up in %d, down in %d, unchanged (within 0.05 point) in %d' % (up, down, n - up - down))
big_down = [r for r in use if r['p_CHJ'] < r['p_CH'] - 5]
for r in big_down:
    print('   down by more than 5 points: %-26s %5.1f%% -> %5.1f%% (place %d -> %d)' % (r['compound'][:26], r['p_CH'], r['p_CHJ'], r['rank_CH'], r['rank_CHJ']))
vals = sorted(r['J_values'] for r in use)
print('couplings per experiment: median %d (from %d to %d)' % (vals[n // 2], vals[0], vals[-1]))
json.dump(dict(n=n, experiments=[r['experiment'] for r in use], summary=summary, better=[r['compound'] for r in better], worse=[r['compound'] for r in worse]),
          open(os.path.join(out, 'jpreview_summary.json'), 'w'), indent=1)
if records_csv:                                   # the same compounds in the shift-only benchmark at the same level of theory
    used = set(r['experiment'] for r in use)
    lvl = [r for r in csv.DictReader(open(records_csv, encoding='utf-8')) if r['level'] == 'BF1' and r['experiment'] in used]
    print('same %d compounds in the shift-only benchmark at %s: actual structure first in %d' % (
        len(lvl), lvl[0]['level_name'], sum(1 for r in lvl if r['top_all'] == truth[r['experiment']]['actual_isomer'])))
    rest = sorted(truth[e]['si_compound'] for e in truth if e not in used)
    print('compounds without a coupling output (%d): %s' % (len(rest), '; '.join(rest)))
print('skipped: %d without a coupling term, %d without experimental couplings' % (sum(1 for _, w in skipped if w.startswith('no')), sum(1 for _, w in skipped if w.startswith('desc'))))
