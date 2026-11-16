"""Checks of the coupling matrices written by the third Kaggle run (notebook 1.2).
Usage: jcheck.py <res/core of run 3> <res/core of run 2> <folder with jref.py> <out folder>"""
import sys, os, gzip, json, glob, collections, statistics, csv
core3, core2, reflib, out = sys.argv[1:5]
sys.path.insert(0, reflib)
import jref


def jl(folder, pattern):
    for p in sorted(glob.glob(os.path.join(folder, pattern))):
        with gzip.open(p, 'rt', encoding='utf-8') as fh:
            for line in fh:
                yield json.loads(line)


J = collections.defaultdict(dict)
for r in jl(core3, 'jcoup_dp4x2_*.jsonl.gz'):
    J[r['path']][r['job']] = r
G2 = collections.defaultdict(dict)
for r in jl(core2, 'gaussian_dp4x2_*.jsonl.gz'):
    if '_J/' in r['path']:
        G2[r['path']][r['job']] = r
res = []
ok = lambda name, cond, detail='': res.append((name, bool(cond), str(detail)))

ok('outputs with coupling matrices: %d' % len(J), len(J) == 3797)
spin2 = set(p for p, steps in G2.items() if any('spinspin' in (s.get('route') or '').lower() for s in steps.values()))
ok('the same %d outputs that run 2 recorded as spin-spin jobs' % len(spin2), set(J) == spin2, len(set(J) ^ spin2))
shape = collections.Counter((tuple(sorted(steps)), tuple(tuple(sorted(steps[k]['J'])) for k in sorted(steps)), tuple(len(steps[k]['jt']) for k in sorted(steps))) for steps in J.values())
ok('every output: step 1 with the Fermi-contact matrix, step 2 with Fermi-contact and total (ten titles)', list(shape) == [((0, 1), (('FC',), ('FC', 'total')), (2, 10))], shape.most_common(3))
ok('no matrix flagged for its layout', not any('bad' in s for steps in J.values() for s in steps.values()))
titles = collections.Counter(tuple(s['jt']) for steps in J.values() for s in steps.values())
ok('titles of the matrices as expected', set(titles) == {('FC:K', 'FC:J'), tuple(n + ':' + k for n in ('FC', 'SD', 'PSO', 'DSO', 'total') for k in 'KJ')}, titles.most_common(3))
# layout found in the deposit: with ReadAtoms (atoms=H) Gaussian 09 prints the matrix for every atom, with zeros in the
# rows and columns of the atoms left out; in step 1 the diagonal of the chosen atoms holds one large constant
bad_lab = bad_zero = bad_len = bad_diag1 = bad_diag2 = 0
sizes, diag_vals = collections.Counter(), collections.Counter()
for p, steps in J.items():
    el = steps[1]['el'].split()
    nat = len(el)
    for k, name in ((0, 'FC'), (1, 'FC'), (1, 'total')):
        b = steps[k]['J'][name]
        v = b['v']
        if b['lab'] != list(range(1, nat + 1)):
            bad_lab += 1
            continue
        if len(v) != nat * (nat + 1) // 2:
            bad_len += 1
            continue
        for i in range(nat):
            base = i * (i + 1) // 2
            if el[i] != 'H':
                if any(v[base + j] != 0.0 for j in range(i + 1)) or any(v[r * (r + 1) // 2 + i] != 0.0 for r in range(i, nat)):
                    bad_zero += 1
                    break
        d = [v[i * (i + 1) // 2 + i] for i in range(nat)]
        if k == 0:
            diag_vals.update(round(d[i] / 1e9) for i in range(nat) if el[i] == 'H')
            bad_diag1 += any(d[i] != 0.0 for i in range(nat) if el[i] != 'H')
        else:
            bad_diag2 += any(x != 0.0 for x in d)
    sizes[sum(1 for e in el if e == 'H')] += 1
ok('row numbers run over every atom of the molecule in all %d matrices' % (3 * len(J)), bad_lab == 0 and bad_len == 0, (bad_lab, bad_len))
ok('rows and columns of atoms other than hydrogen hold only zeros (%d to %d hydrogen atoms per molecule)' % (min(sizes), max(sizes)), bad_zero == 0, bad_zero)
ok('step 2: diagonals are zero in both matrices', bad_diag2 == 0, bad_diag2)
ok('step 1: the diagonal holds one constant on the hydrogen atoms (%s GHz) and zero elsewhere' % ', '.join('%d' % x for x in sorted(diag_vals)),
   bad_diag1 == 0 and len(diag_vals) == 1, diag_vals.most_common(3))
# Fermi-contact part printed twice
dmax, ndiff = 0.0, 0
for p, steps in J.items():
    a, b = steps[0]['J']['FC']['v'], steps[1]['J']['FC']['v']
    nat = len(steps[1]['J']['FC']['lab'])
    diag = set(i * (i + 1) // 2 + i for i in range(nat))
    d = max(abs(x - y) for i, (x, y) in enumerate(zip(a, b)) if i not in diag) if len(a) == len(b) else float('inf')
    dmax = max(dmax, d)
    ndiff += d > 0
ok('off the diagonal, the Fermi-contact matrix of step 2 repeats that of step 1 (largest difference %.2g Hz, %d outputs differ)' % (dmax, ndiff), dmax <= 1e-4, dmax)
# energies and termination against run 2
de = max(abs(steps[k].get('scf', 9e9) - G2[p][k].get('scf', -9e9)) for p, steps in J.items() for k in (0, 1))
norm = sum(1 for steps in J.values() if steps[0].get('norm') and steps[1].get('norm'))
ok('SCF energies equal to those recorded in run 2 (largest difference %.1e hartree)' % de, de < 1e-9, de)
ok('both steps ended normally in %d of %d outputs' % (norm, len(J)), norm == len(J), norm)
# text kept as printed, read with the reader of PyDP4
T = {'FC': 'Fermi Contact (FC) contribution to J (Hz)', 'total': 'Total nuclear spin-spin coupling J (Hz)'}
raw = list(jl(core3, 'jraw_dp4x2_*.jsonl.gz'))
n_ok = 0
folders = collections.Counter()
for r in raw:
    ref, lab = jref.read_cmatrix(r['text'] + '\n', T[r['name']])
    b = J[r['path']][r['job']]['J'][r['name']]
    n = len(b['lab'])
    m = [[0.0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1):
            m[i][j] = m[j][i] = b['v'][k]
            k += 1
    for i in range(n):                                # the reader of PyDP4 sets the diagonal to zero
        m[i][i] = 0.0
    n_ok += ([int(x) for x in lab] == b['lab'] and ref == m)
    folders[r['path'].rsplit('/', 1)[0]] += 1
ok('%d matrices kept as printed text (two outputs in each of %d folders): the reader of PyDP4 gives the stored numbers for %d' % (len(raw), len(folders), n_ok),
   n_ok == len(raw) == 150 and set(folders.values()) == {6} and all(r['why'] == 'sample' for r in raw), collections.Counter(folders.values()))
# what the numbers look like
tot, fc, rest = [], [], []
for p, steps in J.items():
    a, b = steps[1]['J']['total'], steps[1]['J']['FC']
    el = steps[1]['el'].split()
    hs = [i for i, e in enumerate(el) if e == 'H']
    for x, i in enumerate(hs):
        for j in hs[:x]:
            k = i * (i + 1) // 2 + j
            tot.append(a['v'][k]); fc.append(b['v'][k]); rest.append(a['v'][k] - b['v'][k])
q = lambda v, f: sorted(v)[int(f * (len(v) - 1))]
print('proton-proton couplings read: %d pairs in %d conformers' % (len(tot), len(J)))
print('total J (Hz): min %.2f, 1%% %.2f, median %.3f, 99%% %.2f, max %.2f | pairs above 1 Hz in size: %.1f%% (%d)' % (
    min(tot), q(tot, 0.01), q(tot, 0.5), q(tot, 0.99), max(tot), 100 * sum(abs(x) > 1 for x in tot) / len(tot), sum(abs(x) > 1 for x in tot)))
print('total minus Fermi contact (Hz): min %.2f, 1%% %.3f, median %.4f, 99%% %.3f, max %.2f' % (min(rest), q(rest, 0.01), q(rest, 0.5), q(rest, 0.99), max(rest)))
strong = [(t_, f_) for t_, f_ in zip(tot, fc) if abs(t_) > 1]
ma = sum(abs(t_ - f_) for t_, f_ in strong) / len(strong)
sx = sum(f_ for _, f_ in strong); sy = sum(t_ for t_, _ in strong); n = len(strong)
sxx = sum(f_ * f_ for _, f_ in strong); sxy = sum(t_ * f_ for t_, f_ in strong)
slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
print('pairs above 1 Hz: mean |total - FC| %.3f Hz, largest %.2f Hz | total = %.4f x FC %+.3f Hz by least squares' % (ma, max(abs(t_ - f_) for t_, f_ in strong), slope, (sy - slope * sx) / n))
arch = collections.Counter(p.split('BasisFunc/')[-1].split('.zip!/')[0] for p in J)
fold = collections.Counter(p.split('BasisFunc/')[-1].rsplit('/', 1)[0] for p in J)
print('archives: %d | folders: %d | conformers per folder: %d to %d' % (len(arch), len(fold), min(fold.values()), max(fold.values())))
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'jcheck_folders.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['folder', 'outputs'])
    for f, n in sorted(fold.items()):
        w.writerow([f, n])
for name, good, detail in res:
    print(('ok   ' if good else 'FAIL ') + name + ('' if good else ' | ' + detail[:300]))
print('%d checks, %d failed' % (len(res), sum(1 for r in res if not r[1])))
