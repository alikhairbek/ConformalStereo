"""What the coupling constants add: ranking, confidence and conformal sets on the compounds with experimental couplings.
Usage: jeval.py <jbench_candidates.csv> <jbench_couplings.csv> <jpreview_outputs.csv> <jmodel_recovered.json> <out folder> [splits] [seed]

Probabilities are the softmax of the candidates' log-likelihoods. Conformal score: 1 - probability. Families (molecule
archives of the deposit) are the exchangeable units, as in the pooled calibration: a calibration score weighs one over
the number of records of its family, coverage and set size are averaged per test family. Two schemes: random halves of
the families (calibration / test), and every family in turn against all the others."""
import sys, os, csv, math, random, collections
import json
table, couplings, preview, model_json, out = sys.argv[1:6]
splits = int(sys.argv[6]) if len(sys.argv) > 6 else 500
seed = int(sys.argv[7]) if len(sys.argv) > 7 else 20261114


def softmax(v):
    m = max(v)
    w = [math.exp(x - m) for x in v]
    z = sum(w)
    return [x / z for x in w]


def threshold(items, target):
    if not items:
        return float('inf')
    cnt = collections.Counter(f for _, f in items)
    G = len(cnt)
    acc = 0.0
    for s, f in sorted(items):
        acc += 1.0 / cnt[f] / (G + 1)
        if acc >= target - 1e-12:
            return s
    return float('inf')


def fitted_temperature(recs):
    def nll(T):
        s = 0.0
        for L, k in recs:
            L2 = [x / T for x in L]
            m = max(L2)
            s -= L2[k] - m - math.log(sum(math.exp(x - m) for x in L2))
        return s
    lo, hi = 0.05, 200.0
    for _ in range(90):
        a, b = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if nll(a) < nll(b):
            hi = b
        else:
            lo = a
    return (lo + hi) / 2


raw = collections.OrderedDict()
for r in csv.DictReader(open(table)):
    raw.setdefault(r['experiment'], []).append(r)
f = lambda r, k: float(r[k])
VARIANTS = collections.OrderedDict([
    ('shifts (C and H)', lambda r: f(r, 'logL_C') + f(r, 'logL_H')),
    ('shifts + J, Gaussian of the source', lambda r: f(r, 'logL_C') + f(r, 'logL_H') + f(r, 'logL_J_gauss')),
    ('shifts + J, Fermi-contact part', lambda r: f(r, 'logL_C') + f(r, 'logL_H') + f(r, 'logL_J_fc')),
    ('shifts + J, recovered density', lambda r: f(r, 'logL_C') + f(r, 'logL_H') + f(r, 'logL_J_rec')),
    ('J alone, Gaussian of the source', lambda r: f(r, 'logL_J_gauss')),
    ('J alone, recovered density', lambda r: f(r, 'logL_J_rec')),
    ('shifts, energies averaged (4 kJ/mol)', lambda r: f(r, 'E4_CH')),
    ('shifts + J Gaussian, energies averaged', lambda r: f(r, 'E4_CHJ_gauss')),
    ('shifts + J recovered, energies averaged', lambda r: f(r, 'E4_CHJ_rec')),
])
recs = {}
for name, fn in VARIANTS.items():
    lst = []
    for exp, rr in raw.items():
        L = [fn(r) for r in rr]
        p = softmax(L)
        k = [int(r['is_actual']) for r in rr].index(1)
        lst.append(dict(exp=exp, compound=rr[0]['compound'], family=rr[0]['archive'], L=L, k=k, K=len(rr), p=p, t=1 - p[k], s=sorted(1 - x for x in p),
                        ptop=max(p), rank=1 + sum(1 for x in p if x > p[k]), pa=p[k]))
    recs[name] = lst
n = len(raw)
fams = sorted({rr[0]['archive'] for rr in raw.values()})
print('compounds %d, candidates %d, families %d | couplings per compound: median %d' % (
    n, sum(len(rr) for rr in raw.values()), len(fams), sorted(int(rr[[int(r['is_actual']) for r in rr].index(1)]['couplings']) for rr in raw.values())[n // 2]))
os.makedirs(out, exist_ok=True)
w = csv.writer(open(os.path.join(out, 'jeval_overview.csv'), 'w', newline=''))
w.writerow(['data_and_model', 'compounds', 'actual_first', 'calls_99', 'calls_99_right', 'actual_at_1pct_or_less', 'mean_probability_of_actual', 'temperature'])
print()
print('%-42s %12s %18s %12s %10s %6s' % ('data and model', 'actual first', 'calls >=99%: right', 'actual <=1%', 'mean p', 'T'))
for name, lst in recs.items():
    first = sum(r['rank'] == 1 for r in lst)
    conf = [r for r in lst if r['ptop'] >= 0.99]
    right = sum(r['rank'] == 1 for r in conf)
    low = sum(r['pa'] <= 0.01 for r in lst)
    T = fitted_temperature([(r['L'], r['k']) for r in lst])
    mp = 100 * sum(r['pa'] for r in lst) / n
    w.writerow([name, n, first, len(conf), right, low, '%.1f' % mp, '%.2f' % T])
    print('%-42s %5d of %-3d %10d / %-5d %10d %9.1f%% %7.2f' % (name, first, n, right, len(conf), low, mp, T))
base = {r['exp']: r for r in recs['shifts (C and H)']}
print()
for name in ('shifts + J, Gaussian of the source', 'shifts + J, recovered density', 'shifts + J, Fermi-contact part', 'shifts, energies averaged (4 kJ/mol)',
             'shifts + J recovered, energies averaged'):
    better = [r for r in recs[name] if r['rank'] < base[r['exp']]['rank']]
    worse = [r for r in recs[name] if r['rank'] > base[r['exp']]['rank']]
    print('%-42s place of the actual structure against shifts alone: better in %d, worse in %d' % (name, len(better), len(worse)))
    for tag, group in (('better', better), ('worse', worse)):
        for r in group:
            print('      %-6s %-30s %2d candidates: place %d -> %d, probability %5.1f%% -> %5.1f%%' % (tag, r['compound'][:30], r['K'], base[r['exp']]['rank'], r['rank'],
                                                                                               100 * base[r['exp']]['pa'], 100 * r['pa']))
# how large are the coupling errors of the actual structures?
act = [r for r in csv.DictReader(open(couplings)) if r['is_actual'] == '1']
errs = sorted(abs(float(r['experimental_Hz']) - float(r['corrected_Hz'])) for r in act)
nop = sum(1 for r in act if r['partner'] == '')
q = lambda f_: errs[int(f_ * (len(errs) - 1))]
print()
print('couplings of the actual structures: %d (without a calculated partner: %d) | absolute error after the linear correction of the source:' % (len(errs), nop))
print('   median %.2f Hz, 75%% %.2f, 90%% %.2f, 95%% %.2f, largest %.2f | above 1.6 Hz (two standard deviations of the Gaussian of the source): %d (%.0f%%); above 3 Hz: %d (%.0f%%)' % (
    q(0.5), q(0.75), q(0.9), q(0.95), errs[-1], sum(e > 1.6 for e in errs), 100.0 * sum(e > 1.6 for e in errs) / len(errs), sum(e > 3 for e in errs), 100.0 * sum(e > 3 for e in errs) / len(errs)))
# printed numbers of the deposit outputs next to the recomputed ones
pv = {}
for r in csv.DictReader(open(preview)):
    if r['experiment'] and r['rank_CH']:
        cur = pv.get(r['experiment'])
        if cur is None or (not cur['file'].endswith('NMRj') and r['file'].endswith('NMRj')):
            pv[r['experiment']] = r
gaussian_outputs = set(json.load(open(model_json))['gaussian_outputs'])
ren = {exp: any(r['renumbered'] == '1' for r in rr) for exp, rr in raw.items()}
mine = {r['exp']: r for r in recs['shifts + J, recovered density']}
mine_g = {r['exp']: r for r in recs['shifts + J, Gaussian of the source']}
wp2 = csv.writer(open(os.path.join(out, 'jeval_printed_vs_recomputed.csv'), 'w', newline=''))
wp2.writerow(['compound', 'candidates', 'renumbered_candidates', 'model_of_the_deposit_output', 'printed_place_shifts', 'printed_place_with_couplings',
              'recomputed_place_shifts', 'recomputed_place_with_couplings_recovered_density', 'recomputed_place_with_couplings_gaussian'])
diff = []
pf = pw = 0
for exp in raw:
    if exp not in pv:
        continue
    model = 'Gaussian' if pv[exp]['file'] in gaussian_outputs else 'recovered density'
    row = [base[exp]['compound'], base[exp]['K'], int(ren[exp]), model, int(pv[exp]['rank_CH']), int(pv[exp]['rank_CHJ']), base[exp]['rank'], mine[exp]['rank'], mine_g[exp]['rank']]
    wp2.writerow(row)
    pf += row[4] == 1
    pw += row[5] == 1
    own = row[8] if model == 'Gaussian' else row[7]
    if row[5] != own:
        diff.append(row + [own])
print()
print('deposit outputs as printed, same %d compounds: actual structure first with shifts %d, with couplings %d' % (len(pv), pf, pw))
print('place of the actual structure with couplings, printed against recomputed with the model of that output: differs for %d compounds' % len(diff))
for row in diff:
    print('   %-30s %2d candidates | printed %d, recomputed %d | %s%s' % (row[0][:30], row[1], row[5], row[9], row[3], ', candidates numbered differently' if row[2] else ''))

# conformal sets
wr = csv.writer(open(os.path.join(out, 'jeval_conformal.csv'), 'w', newline=''))
wr.writerow(['scheme', 'target', 'data_and_model', 'coverage', 'mean_set_size', 'single_structure_sets', 'empty_sets'])
show = ['shifts (C and H)', 'shifts + J, Gaussian of the source', 'shifts + J, recovered density', 'shifts + J, Fermi-contact part',
        'shifts, energies averaged (4 kJ/mol)', 'shifts + J Gaussian, energies averaged', 'shifts + J recovered, energies averaged']


def summarise(test, thr):
    by = collections.defaultdict(list)
    for r in test:
        size = sum(1 for x in r['s'] if x <= thr)
        by[r['family']].append((r['t'] <= thr, size, size == 1, size == 0))
    G = len(by)
    return [sum(sum(v[i] for v in vs) / len(vs) for vs in by.values()) / G for i in range(4)]


for target in (0.90, 0.80):
    print()
    print('conformal sets at %d%%: random halves of the families (%d splits) | each family against all the others' % (100 * target, splits))
    print('%-42s %9s %8s %8s %7s | %9s %8s %8s %7s' % ('data and model', 'coverage', 'size', 'single', 'empty', 'coverage', 'size', 'single', 'empty'))
    rng = random.Random(seed)
    acc = {name: [0.0] * 4 for name in show}
    for _ in range(splits):
        fl = fams[:]
        rng.shuffle(fl)
        cal = set(fl[:len(fl) // 2])
        for name in show:
            lst = recs[name]
            thr = threshold([(r['t'], r['family']) for r in lst if r['family'] in cal], target)
            v = summarise([r for r in lst if r['family'] not in cal], thr)
            for i in range(4):
                acc[name][i] += v[i] / splits
    fam_size = {}
    for name in show:
        lst = recs[name]
        cov = size = single = empty = 0.0
        fam_size[name] = []
        for fam in fams:
            thr = threshold([(r['t'], r['family']) for r in lst if r['family'] != fam], target)
            v = summarise([r for r in lst if r['family'] == fam], thr)
            cov += v[0] / len(fams); size += v[1] / len(fams); single += v[2] / len(fams); empty += v[3] / len(fams)
            fam_size[name].append(v[1])
        a = acc[name]
        wr.writerow(['random halves', target, name, '%.4f' % a[0], '%.3f' % a[1], '%.4f' % a[2], '%.4f' % a[3]])
        wr.writerow(['one family out', target, name, '%.4f' % cov, '%.3f' % size, '%.4f' % single, '%.4f' % empty])
        print('%-42s %8.1f%% %8.2f %7.1f%% %6.1f%% | %8.1f%% %8.2f %7.1f%% %6.1f%%' % (name, 100 * a[0], a[1], 100 * a[2], 100 * a[3], 100 * cov, size, 100 * single, 100 * empty))
    rb = random.Random(seed + 1)
    for name in ('shifts + J, recovered density', 'shifts, energies averaged (4 kJ/mol)', 'shifts + J recovered, energies averaged', 'shifts + J, Gaussian of the source'):
        d = [x - y for x, y in zip(fam_size[name], fam_size['shifts (C and H)'])]
        boot = sorted(sum(rb.choice(d) for _ in d) / len(d) for _ in range(4000))
        print('   set size against shifts alone, each family against the others: %-42s %+.2f (95%% interval over families %+.2f to %+.2f)' % (
            name, sum(d) / len(d), boot[int(0.025 * len(boot))], boot[int(0.975 * len(boot))]))
# per compound
wp = csv.writer(open(os.path.join(out, 'jeval_compounds.csv'), 'w', newline=''))
wp.writerow(['compound', 'experiment', 'family', 'candidates'] + ['place: ' + v for v in VARIANTS] + ['probability: ' + v for v in VARIANTS])
for i, exp in enumerate(raw):
    wp.writerow([recs[show[0]][i]['compound'], exp, recs[show[0]][i]['family'], recs[show[0]][i]['K']] + [recs[v][i]['rank'] for v in VARIANTS] + ['%.4g' % recs[v][i]['pa'] for v in VARIANTS])
