"""2023 SI: MM-DP4+ at the 36 levels for the molecules of section 7 that come as a pair of isomers, each with its own
experimental data (each isomer is the alternative of the other). Recomputed from the printed shieldings with the
parameters shipped in DP4+App and compared with the whole-per-cent probabilities printed in section 8.
Stage 1: inputs by rule (carbons at or above 100 ppm are sp2, protons follow their carbon, no exchangeable pairs).
Stage 2: for the cases with records that are not reproduced, the unprinted inputs are searched (infer.py); the result
is kept only if fewer records are left unreproduced.
Writes jnp_s7_cases.csv, jnp_s7_records.csv, jnp_s7_candidates.csv.
Usage: jnp_verify_s7.py <lib dir> <scripts dir> <app_parameters.csv> <jnp_tables.json> <jnp_section8.csv> <out dir> [draws]"""
import sys, os, csv, json, collections
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[2])
import numpy as np
import dp4plus as D
import joc_data as J
import mmdp4 as M
import infer as I
import verify as V
par_app = M.load_parameters(sys.argv[3], 'MM'); T = json.load(open(sys.argv[4])); out_dir = sys.argv[6]; draws = int(sys.argv[7]) if len(sys.argv) > 7 else 200
pub = collections.OrderedDict()
for r in csv.DictReader(open(sys.argv[5])): pub.setdefault(r['molecule'], {})[r['level']] = float(r['probability'])
mols = list(pub); s7 = [t for t in T if t['line'] < 26358]; assert len(s7) == 3 * len(mols)
tab = {m: s7[3 * i:3 * i + 3] for i, m in enumerate(mols)}
def shield(m, code):
    for t in tab[m]:
        if code in t['levels']: return np.array([r[3 + t['levels'].index(code)] for r in t['rows']], float)
codes = [c for t in tab[mols[0]] for c in t['levels']]; pars = {c: par_app[M.app_level(c)] for c in codes}
KEY = ('DP4+', 'All')
pairs_of_isomers = [(a, b) for a in mols for b in mols if a[:-1] == b[:-1] and a != b]
singles = [m for m in mols if not any(m == a for a, _ in pairs_of_isomers)]
def assess(rows, sh, sp2, pairs, a):
    exp = [r['exp'] for r in rows]; nuc = [r['nuc'] for r in rows]; out = []
    for c in codes:
        s = [M.exchange(x, exp, pairs) for x in sh[c]] if pairs else sh[c]
        pr = M.probabilities(D, s, exp, nuc, sp2, pars[c]); p = 100 * pr[KEY]
        st, rg = V.status(float(p[0]), pub[a][c], lambda: V.rounding_range(D, sh[c], exp, nuc, sp2, pars[c], KEY, 0, 0.005, pairs=pairs, n=draws, seed=17))
        out.append((c, [float(x) for x in p], st, rg))
    return out
cases = []; rec = []; cand_rows = []
for a, b in pairs_of_isomers:
    ra, rb = tab[a][0]['rows'], tab[b][0]['rows']; assert [r[:2] for r in ra] == [r[:2] for r in rb]
    rows = [dict(nuc=r[1], label=r[0], exp=r[2]) for r in ra]; lab = lambda i: rows[i]['nuc'] + rows[i]['label']
    sh = {c: [shield(a, c), shield(b, c)] for c in codes}
    if all(np.array_equal(sh[c][0], sh[c][1]) for c in codes):
        cases.append(dict(case=a, alternative=b, status='excluded', note='the SI prints the same shieldings and experimental values for both isomers', reproduced=0, consistent=0, not_reproduced=0)); continue
    sp2 = J.sp2_flags(rows, 100.0, None)
    res = assess(rows, sh, sp2, [], a); nbad = sum(st == 'not reproduced' for _, _, st, _ in res); used = 'by rule'; flips = []; pairs = []
    if nbad:
        targets = [(c, KEY, 0, pub[a][c]) for c in codes]
        inf = I.infer(D, rows, sh, sp2, targets, pars)
        res2 = assess(rows, sh, inf['sp2'], inf['pairs'], a); nbad2 = sum(st == 'not reproduced' for _, _, st, _ in res2)
        if nbad2 < nbad: res, used, flips, pairs, sp2 = res2, 'recovered from the printed probabilities', inf['flips'], inf['pairs'], inf['sp2']
    n = collections.Counter(st for _, _, st, _ in res)
    status = 'verified' if n['not reproduced'] == 0 else ('partly verified' if n['not reproduced'] <= 3 else 'unverified')
    cases.append(dict(case=a, alternative=b, status=status, inputs=used, classes_changed=' '.join('%s->%s' % (lab(i), 'sp2' if sp2[i] else 'sp3') for i in flips), exchangeable_pairs=' '.join('%s/%s' % (lab(i), lab(j)) for i, j in pairs),
                      reproduced=n['reproduced'], consistent=n['consistent'], not_reproduced=n['not reproduced'], largest_difference=round(max(abs(p[0] - pub[a][c]) for c, p, _, _ in res), 1), note=''))
    for c, p, st, rg in res:
        rec.append(dict(case=a, alternative=b, level=c, published=pub[a][c], recomputed=round(p[0], 2), status=st, range_low='' if rg is None else round(rg[0], 1), range_high='' if rg is None else round(rg[1], 1)))
        for k, cc in enumerate((a, b)): cand_rows.append(dict(case=a, level=c, candidate=cc, is_correct=int(k == 0), probability=round(p[k], 4), record_status=st))
def write(name, rows_, keys=None):
    keys = keys or list(rows_[0].keys())
    w = csv.DictWriter(open(os.path.join(out_dir, name), 'w', newline=''), fieldnames=keys, extrasaction='ignore'); w.writeheader(); w.writerows(rows_)
write('jnp_s7_cases.csv', cases, ['case', 'alternative', 'status', 'inputs', 'classes_changed', 'exchangeable_pairs', 'reproduced', 'consistent', 'not_reproduced', 'largest_difference', 'note'])
write('jnp_s7_records.csv', rec); write('jnp_s7_candidates.csv', cand_rows)
print('molecules in section 7: %d | single molecules, alternative isomer not printed: %s | cases from pairs of isomers: %d' % (len(mols), ' '.join(singles), len(cases)))
print('cases:', dict(collections.Counter(c['status'] for c in cases)), '| inputs:', dict(collections.Counter(c.get('inputs', '-') for c in cases)))
n = collections.Counter(r['status'] for r in rec)
print('records (case, level): %d | reproduced within 1 point %d (%.1f%%) | consistent within the rounding of the printed shieldings %d | not reproduced %d (%.1f%%)' % (len(rec), n['reproduced'], 100 * n['reproduced'] / len(rec), n['consistent'], n['not reproduced'], 100 * n['not reproduced'] / len(rec)))
for c in cases:
    if c['status'] != 'verified' or c.get('inputs') != 'by rule': print('  %-5s %-16s reproduced %2s, consistent %2s, not %2s | largest difference %5s | %s%s%s' % (c['case'], c['status'], c['reproduced'], c['consistent'], c['not_reproduced'], c.get('largest_difference', ''), c.get('inputs', ''), (' | classes: ' + c['classes_changed']) if c.get('classes_changed') else '', (' | pairs: ' + c['exchangeable_pairs']) if c.get('exchangeable_pairs') else '') + (' | ' + c['note'] if c['note'] else ''))
