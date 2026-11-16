"""Confidence audit and conformal preview on the DP4-AI deposit (recomputed DP4, reported structure = isomer 1).
Usage: dp4ai_audit.py <out_dir> [n_splits] [seed]"""
import sys, os, csv, math, random, collections
out = sys.argv[1]; n_splits = int(sys.argv[2]) if len(sys.argv) > 2 else 2000; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261113
PROT = ['MM', 'DFT_EB3LYP', 'DFT_EmPW', 'DFT_EM062X']
TRUE = {}                                           # case -> number of the reported structure among the candidates
def load(tag):
    cases = [r for r in csv.DictReader(open(os.path.join(out, 'dp4ai_cases%s.csv' % tag)))]
    probs = collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(out, 'dp4ai_probabilities%s.csv' % tag))):
        probs[(r['case'], r['protocol'], r['data'])][int(r['isomer'])] = float(r['probability'])
        if r['is_true'] == '1': TRUE[r['case']] = int(r['isomer'])
    return cases, probs
cases, probs = load(''); cases_old, probs_old = load('_NMRold')
clean = lambda r: r['qc_clean'] == '1' and r['exact_counts'] == '1'
print('=== 1. Stated confidence against outcome (combined 13C + 1H where both exist; description <case>NMR) ===')
print('%-11s %-6s | cases | top-1 right | top prob >= 99%%: right/all | 90-99%%: right/all | < 90%%: right/all | true isomer <= 1%% | nominal 90%% set: covers truth, mean size' % ('protocol', 'subset'))
pooled = collections.Counter()
for prot in PROT:
    for name, keep in (('all', lambda r: True), ('clean', clean)):
        rs = [r for r in cases if r['protocol'] == prot and r['data'] == 'all' and keep(r)]
        b = collections.Counter(); cov = 0; size = 0
        for r in rs:
            pt, top, ok = float(r['p_true']), float(r['p_top']), r['top1_correct'] == '1'
            k = 'hi' if top >= 0.99 else ('mid' if top >= 0.90 else 'lo')
            b[k] += 1; b[k + '_ok'] += ok
            pr = probs[(r['case'], prot, 'all')]; acc = 0.0; chosen = []
            for iso, p in sorted(pr.items(), key=lambda kv: -kv[1]):
                chosen.append(iso); acc += p
                if acc >= 0.9: break
            cov += TRUE[r['case']] in chosen; size += len(chosen)
        n = len(rs)
        print('%-11s %-6s | %5d | %3d (%2.0f%%)   | %14d/%-10d | %9d/%-7d | %8d/%-7d | %17d | %12d/%d (%2.0f%%), %.2f' % (prot, name, n, sum(r['top1_correct'] == '1' for r in rs), 100 * sum(r['top1_correct'] == '1' for r in rs) / n,
              b['hi_ok'], b['hi'], b['mid_ok'], b['mid'], b['lo_ok'], b['lo'], sum(float(r['p_true']) <= 0.01 for r in rs), cov, n, 100 * cov / n, size / n))
        if name == 'clean':
            pooled['hi'] += b['hi']; pooled['hi_ok'] += b['hi_ok']; pooled['n'] += n; pooled['low'] += sum(float(r['p_true']) <= 0.01 for r in rs); pooled['cov'] += cov
print('pooled over the four protocols, clean cases: calls at >= 99%% right in %d of %d (%.1f%%); true isomer <= 1%% in %d of %d case-protocol pairs (%.1f%%); nominal 90%% sets cover %d of %d (%.1f%%)' % (
      pooled['hi_ok'], pooled['hi'], 100 * pooled['hi_ok'] / pooled['hi'], pooled['low'], pooled['n'], 100 * pooled['low'] / pooled['n'], pooled['cov'], pooled['n'], 100 * pooled['cov'] / pooled['n']))

print('\n=== 2. One molecule, eight versions of the same analysis (4 protocols x 2 description files) ===')
ver = collections.defaultdict(dict)
for tag, cs in (('NMR', cases), ('NMRold', cases_old)):
    for r in cs:
        if r['data'] == 'all': ver[r['case']][(r['protocol'], tag)] = r
n_cases = len(ver); st = collections.Counter(); contradict = []; flips = []
for c, d in ver.items():
    tops = [(k, int(r['top_isomer']), float(r['p_top'])) for k, r in d.items()]
    conf = set(t for _, t, p in tops if p >= 0.99)
    st['versions'] += len(tops)
    st['always right'] += all(t == TRUE[c] for _, t, _ in tops)
    st['never right'] += all(t != TRUE[c] for _, t, _ in tops)
    st['some version puts >= 99%% on a wrong isomer'] += any(t != TRUE[c] and p >= 0.99 for _, t, p in tops)
    if len(conf) > 1: contradict.append('%s (%s)' % (c, ', '.join('isomer %d' % x for x in sorted(conf))))
    pts = [float(r['p_true']) for r in d.values()]
    if max(pts) >= 0.99 and min(pts) <= 0.01: flips.append(c)
print('molecules: %d, analyses: %d' % (n_cases, st['versions']))
print('reported structure ranked first in every version: %d | in no version: %d | in some but not all: %d' % (st['always right'], st['never right'], n_cases - st['always right'] - st['never right']))
print('molecules where at least one version puts >= 99%% on a wrong isomer: %d' % st['some version puts >= 99%% on a wrong isomer'])
print('molecules where two versions each put >= 99%% on a different isomer: %d -> %s' % (len(contradict), '; '.join(contradict)))
print('molecules where P(reported structure) is >= 99%% in one version and <= 1%% in another: %d -> %s' % (len(flips), ' '.join(flips)))
same_desc = collections.Counter()
for c, d in ver.items():
    for tag in ('NMR', 'NMRold'):
        tops = [(int(r['top_isomer']), float(r['p_top'])) for (p, t), r in d.items() if t == tag]
        if len(tops) >= 2:
            same_desc[tag + ' molecules'] += 1
            same_desc[tag + ' contradiction'] += len(set(t for t, p in tops if p >= 0.99)) > 1
print('same description file, protocols only: two protocols each >= 99%% on different isomers in %d of %d molecules (NMR) and %d of %d (NMRold)' % (
      same_desc['NMR contradiction'], same_desc['NMR molecules'], same_desc['NMRold contradiction'], same_desc['NMRold molecules']))

print('\n=== 3. Conformal preview (score = 1 - P(reported structure); %d random half/half splits by candidate-set family) ===' % n_splits)
random.seed(seed)
def conformal(rs, pr_of, alpha):
    fams = sorted(set(r['family'] for r in rs)); member = {f: [r for r in rs if r['family'] == f] for f in fams}
    tot = collections.Counter()
    for _ in range(n_splits):
        random.shuffle(fams); half = len(fams) // 2
        cal = [r for f in fams[:half] for r in member[f]]; test = [r for f in fams[half:] for r in member[f]]
        s = sorted(1 - float(r['p_true']) for r in cal); k = math.ceil((len(s) + 1) * (1 - alpha))
        q = s[k - 1] if k <= len(s) else float('inf')
        for r in test:
            pr = pr_of(r); kept = [i for i, p in pr.items() if 1 - p <= q + 1e-12]
            t = TRUE[r['case']]
            tot['n'] += 1; tot['cov'] += t in kept; tot['size'] += len(kept); tot['frac'] += len(kept) / len(pr)
            tot['single_right'] += kept == [t]; tot['single_wrong'] += (len(kept) == 1 and kept != [t]); tot['all'] += len(kept) == len(pr); tot['empty'] += not kept
    n = tot['n']
    return [100 * tot['cov'] / n, tot['size'] / n, 100 * tot['frac'] / n, 100 * tot['single_right'] / n, 100 * tot['single_wrong'] / n, 100 * tot['all'] / n, 100 * tot['empty'] / n]
for alpha in (0.20, 0.10):
    print('target coverage %.0f%% | columns: coverage %%, mean set size, share of candidates kept %%, one candidate left and right %%, one left and wrong %%, all candidates kept %%, none kept %%' % (100 * (1 - alpha)))
    for prot in PROT:
        for name, keep in (('all', lambda r: True), ('clean', clean)):
            rs = [r for r in cases if r['protocol'] == prot and r['data'] == 'all' and keep(r)]
            v = conformal(rs, lambda r: probs[(r['case'], prot, 'all')], alpha)
            print('   %-11s %-6s n=%2d fam=%2d | %5.1f  %5.2f  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f' % ((prot, name, len(rs), len(set(r['family'] for r in rs))) + tuple(v)))
