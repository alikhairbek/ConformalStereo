"""Audit of the probabilities PUBLISHED with the DP4-AI paper (three-Gaussian model, cross-validated; raw-data route
'DP4-AI' and manual route 'pairwise'), and a comparison with my recomputation from the deposit (default single Gaussian).
The reported structure is isomer 1 of the deposit; for the two molecules without deposit labels the sheet's own mark is used.
Usage: dp4ai_published_audit.py <out_dir of this folder> <tables dir of the DP4-AI benchmark package> [n_splits]"""
import sys, os, csv, math, random, collections
out, bench = sys.argv[1:3]; n_splits = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
rows = list(csv.DictReader(open(os.path.join(out, 'dp4ai_published.csv'))))
P = collections.defaultdict(dict); TRUE = {}; marked = {}
for r in rows:
    for meth in ('ai', 'pair'):
        v = r[meth + '_all']
        if v != '': P[(r['level'], meth, r['molecule'])][r['stereo']] = float(v) / 100
    if r['level'] == 'MM' and r['is_true'] == '1': marked[r['molecule']] = r['stereo']
    if r['deposit_isomer'] == '1': TRUE[r['molecule']] = r['stereo']
for m, s in marked.items(): TRUE.setdefault(m, s)
fam = lambda m: {'JB3': 'JB2', 'JB13B': 'JB13A', 'NL1B': 'NL1A', 'NL2B': 'NL2A', 'NP3B': 'NP3A', 'TS3B': 'TS3A'}.get(m, m)
print('molecules: %d | reported structure taken from the deposit for %d, from the sheet mark for %s' % (len(TRUE), sum(1 for m in TRUE if m not in ('NP5', 'TP4')), 'NP5 TP4'))
print('\n=== 1. Published probabilities against outcome ===')
print('%-5s %-8s | cases | top-1 right | calls >= 99%%: right/all | 90-99%%: right/all | true <= 1%% | nominal 90%% set covers, mean size' % ('level', 'route'))
tot = collections.Counter()
for lvl in ('MM', 'Opt', 'Opt_E'):
    for meth, name in (('ai', 'DP4-AI'), ('pair', 'pairwise')):
        b = collections.Counter(); n = 0; cov = 0; size = 0
        for (l, me, m), pr in P.items():
            if l != lvl or me != meth: continue
            n += 1; t = TRUE[m]; top = max(pr, key=pr.get); ok = top == t
            k = 'hi' if pr[top] >= 0.99 else ('mid' if pr[top] >= 0.90 else 'lo')
            b[k] += 1; b[k + '_ok'] += ok; b['ok'] += ok; b['low'] += pr.get(t, 0) <= 0.01
            acc = 0; ch = []
            for s, p in sorted(pr.items(), key=lambda kv: -kv[1]):
                ch.append(s); acc += p
                if acc >= 0.9: break
            cov += t in ch; size += len(ch)
        print('%-5s %-8s | %5d | %3d (%2.0f%%)   | %10d/%-10d | %8d/%-8d | %10d | %3d/%d (%2.0f%%), %.2f' % (lvl, name, n, b['ok'], 100 * b['ok'] / n, b['hi_ok'], b['hi'], b['mid_ok'], b['mid'], b['low'], cov, n, 100 * cov / n, size / n))
        for k in ('hi', 'hi_ok', 'low'): tot[(meth, k)] += b[k]
        tot[(meth, 'n')] += n; tot[(meth, 'cov')] += cov
for meth, name in (('ai', 'DP4-AI'), ('pair', 'pairwise')):
    print('pooled over the three levels, %s: calls at >= 99%% right in %d of %d (%.1f%%); true isomer <= 1%% in %d of %d (%.1f%%); nominal 90%% sets cover %.1f%%' % (
        name, tot[(meth, 'hi_ok')], tot[(meth, 'hi')], 100 * tot[(meth, 'hi_ok')] / tot[(meth, 'hi')], tot[(meth, 'low')], tot[(meth, 'n')], 100 * tot[(meth, 'low')] / tot[(meth, 'n')], 100 * tot[(meth, 'cov')] / tot[(meth, 'n')]))
print('\n=== 2. Conformal sets on the published probabilities (score = 1 - P(reported structure); %d half/half splits by family) ===' % n_splits)
random.seed(20261114)
def conformal(lvl, meth, alpha):
    items = [(m, pr) for (l, me, m), pr in P.items() if l == lvl and me == meth]
    fams = sorted(set(fam(m) for m, _ in items)); member = {f: [it for it in items if fam(it[0]) == f] for f in fams}; t = collections.Counter()
    for _ in range(n_splits):
        random.shuffle(fams); h = len(fams) // 2
        cal = [it for f in fams[:h] for it in member[f]]; test = [it for f in fams[h:] for it in member[f]]
        s = sorted(1 - pr.get(TRUE[m], 0) for m, pr in cal); k = math.ceil((len(s) + 1) * (1 - alpha)); q = s[k - 1] if k <= len(s) else float('inf')
        for m, pr in test:
            kept = [x for x, p in pr.items() if 1 - p <= q + 1e-12]
            t['n'] += 1; t['cov'] += TRUE[m] in kept; t['frac'] += len(kept) / len(pr); t['size'] += len(kept); t['one_ok'] += kept == [TRUE[m]]; t['all'] += len(kept) == len(pr); t['none'] += not kept
    n = t['n']; return 100 * t['cov'] / n, t['size'] / n, 100 * t['frac'] / n, 100 * t['one_ok'] / n, 100 * t['all'] / n, 100 * t['none'] / n
for alpha in (0.2, 0.1):
    print('target %.0f%% | coverage %%, mean set size, candidates kept %%, one left and right %%, all kept %%, none kept %%' % (100 * (1 - alpha)))
    for lvl in ('MM', 'Opt', 'Opt_E'):
        for meth, name in (('ai', 'DP4-AI'), ('pair', 'pairwise')):
            print('   %-5s %-8s | %5.1f  %5.2f  %5.1f  %5.1f  %5.1f  %5.1f' % ((lvl, name) + conformal(lvl, meth, alpha)))
print('\n=== 3. Published pairwise route against my recomputation (different error models, same data) ===')
mine = collections.defaultdict(dict)
for r in csv.DictReader(open(os.path.join(bench, 'dp4ai_probabilities.csv'))):
    if r['data'] == 'all': mine[(r['protocol'], r['case'])][r['isomer']] = float(r['probability'])
num = collections.defaultdict(dict)
for r in rows:
    if r['deposit_isomer']: num[r['molecule']][r['stereo']] = r['deposit_isomer']
print('%-5s vs %-11s | molecules | same top isomer | mean |difference| of P(reported structure), points' % ('level', 'my protocol'))
for lvl, prots in (('MM', ['MM']), ('Opt', ['DFT_EB3LYP', 'DFT_EmPW', 'DFT_EM062X']), ('Opt_E', ['DFT_EB3LYP', 'DFT_EmPW', 'DFT_EM062X'])):
    for prot in prots:
        n = same = 0; d = 0.0
        for (l, me, m), pr in P.items():
            if l != lvl or me != 'pair' or (prot, m) not in mine or m not in num: continue
            mp = mine[(prot, m)]; pub = {num[m][s]: p for s, p in pr.items() if s in num[m]}
            if set(pub) != set(mp): continue
            n += 1; same += max(pub, key=pub.get) == max(mp, key=mp.get); d += abs(pub['1'] - mp['1'])
        print('%-5s vs %-11s | %9d | %4d (%2.0f%%)      | %.1f' % (lvl, prot, n, same, 100 * same / max(n, 1), 100 * d / max(n, 1)))
print("\n=== 4. Does the result depend on which isomer is taken as the reported structure? ===")
raw = list(csv.DictReader(open(os.path.join(out, 'dp4ai_published_raw.csv'))))
by = collections.defaultdict(lambda: collections.OrderedDict())
for r in raw: by[r['level']].setdefault(r['molecule'], []).append(r)
own = {}
for lvl in ('MM', 'Opt', 'Opt_E'):
    for m, v in by[lvl].items():
        mk = [x['stereo'] for x in v if x['marked'] == '1']
        own[(lvl, m)] = v[0]['stereo'] if lvl == 'Opt' else (mk[0] if mk else v[0]['stereo'])     # in the Opt sheet the mark column sits one row low
print('the mark in a sheet differs from the deposit label for:', sorted((l, m) for (l, m), s in own.items() if s != TRUE[m]))
for label, truth in (('deposit labels', lambda l, m: TRUE[m]), ("each sheet's own marks", lambda l, m: own[(l, m)])):
    for meth, name in (('ai_all', 'DP4-AI'), ('pair_all', 'pairwise')):
        hi = ok = 0
        for lvl in ('MM', 'Opt', 'Opt_E'):
            for m, v in by[lvl].items():
                pr = {x['stereo']: float(x[meth]) for x in v if x[meth] != ''}
                top = max(pr, key=pr.get)
                if pr[top] >= 99: hi += 1; ok += top == truth(lvl, m)
        print('   %-22s %-8s calls at >= 99%% right in %d of %d (%.1f%%)' % (label, name, ok, hi, 100 * ok / hi))
