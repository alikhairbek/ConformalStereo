"""2023 SI: MM-DP4+ for compounds 101-121 (section 9: all candidates of a compound in one table, three levels)
against section 10, which prints the probability of the correct isomer in six formulations. Two of them can be
recomputed with the DP4+App parameters: 'ScTMS' (scaled shifts only) and 'ScTMS+Uns' (scaled and unscaled, the
MM-DP4+ of the program). The third level is SMD/wB97XD/6-31+G** although its columns are headed WB9_C_SMD.
For every case: fit of each candidate to the experimental column, then the six printed values (two formulations at
three levels) by rule, then with recovered inputs if needed.
Usage: jnp_verify_s9.py <lib dir> <scripts dir> <app_parameters.csv> <jnp_tables.json> <jnp_probabilities.csv> <out dir> [draws]"""
import sys, os, csv, json, math, collections
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[2])
import numpy as np
import dp4plus as D
import joc_data as J
import mmdp4 as M
import infer as I
import verify as V
par_app = M.load_parameters(sys.argv[3], 'MM'); T = json.load(open(sys.argv[4])); out_dir = sys.argv[6]; draws = int(sys.argv[7]) if len(sys.argv) > 7 else 200
pub = {r['case']: r for r in csv.DictReader(open(sys.argv[5]))}
s9 = [t for t in T if t['line'] > 26358]
LEV = [('B3L_A_GAS', 'MMFF_B3L', 'B3LYP.6-31G(d,p)'), ('MPW_B_PCM', 'MMFF_mPW', 'mPW1PW91.6-31+G(d,p).PCM'), ('WB9_C_SMD', 'MMFF_wB9', 'wB97XD.6-31+G(d,p).SMD')]
FORM = [('ScTMS', ('sDP4+', 'All')), ('ScTMS+Uns', ('DP4+', 'All'))]
def structure(t):
    ncol = len(t['rows'][0]) - 2
    if ncol == 13: return [0], [[1 + 3 * k + i for i in range(3)] for k in range(4)]
    if ncol == 7: return [0], [[1 + 3 * k + i for i in range(3)] for k in range(2)]
    if ncol == 8: return [0, 4], [[1, 2, 3], [5, 6, 7]]
    if ncol == 14: return [0, 4], [[1, 2, 3], [5, 6, 7], [8, 9, 10], [11, 12, 13]]
    return None, None
def scaled_rms(exp, sig):
    n = len(exp); mx, my = sum(sig) / n, sum(exp) / n; sxx = sum((s - mx) ** 2 for s in sig); sxy = sum((s - mx) * (e - my) for s, e in zip(sig, exp))
    b = sxy / sxx; a = my - b * mx
    return math.sqrt(sum((e - (a + b * s)) ** 2 for s, e in zip(sig, exp)) / n)
cases = []; rec = []; cand_rows = []; number = 100
for t in s9:
    ex, cands = structure(t)
    if not ex: continue                                  # the table of 121e-h: row labels shifted, not used
    number += 1
    for e_i, e in enumerate(ex):
        name = '%d%s' % (number, 'ab'[e_i]); true = e_i
        if name not in pub: continue
        rows = [dict(nuc=r[1], label=r[0], exp=r[2 + e]) for r in t['rows']]; lab = lambda i: rows[i]['nuc'] + rows[i]['label']
        exp = [r['exp'] for r in rows]; nuc = [r['nuc'] for r in rows]
        sh = {code: [np.array([r[2 + c[i]] for r in t['rows']], float) for c in cands] for i, (code, _, _) in enumerate(LEV)}
        pars = {code: par_app[app] for code, _, app in LEV}
        fit = []
        for k in range(len(cands)):
            v = sh['MPW_B_PCM'][k]
            fit.append((scaled_rms([x for x, n_ in zip(exp, nuc) if n_ == 'C'], [x for x, n_ in zip(v, nuc) if n_ == 'C']), scaled_rms([x for x, n_ in zip(exp, nuc) if n_ == 'H'], [x for x, n_ in zip(v, nuc) if n_ == 'H'])))
        bad_fit = [k for k, (c, h) in enumerate(fit) if c > 6 or h > 0.6]
        targets = [(code, key, true, float(pub[name]['%s|%s' % (plab, form)])) for code, plab, _ in LEV for form, key in FORM]
        def assess(sp2, pairs):
            out = []
            for code, key, k, val in targets:
                s = [M.exchange(x, exp, pairs) for x in sh[code]] if pairs else sh[code]
                p = 100 * float(M.probabilities(D, s, exp, nuc, sp2, pars[code])[key][k])
                st, rg = V.status(p, val, lambda: V.rounding_range(D, sh[code], exp, nuc, sp2, pars[code], key, k, 0.005, pairs=pairs, n=draws, seed=17), printed_half_unit=0.005)
                out.append((code, key, val, p, st, rg))
            return out
        sp2 = J.sp2_flags(rows, 100.0, None); res = assess(sp2, []); nbad = sum(x[4] == 'not reproduced' for x in res); used = 'by rule'; flips = []; pairs = []
        if nbad:
            inf = I.infer(D, rows, sh, sp2, targets, pars); res2 = assess(inf['sp2'], inf['pairs']); nbad2 = sum(x[4] == 'not reproduced' for x in res2)
            if nbad2 < nbad: res, used, flips, pairs, sp2 = res2, 'recovered from the printed probabilities', inf['flips'], inf['pairs'], inf['sp2']
        n = collections.Counter(x[4] for x in res)
        if bad_fit: status = 'excluded'; note = 'candidate %s does not fit the experimental column (13C/1H error %s ppm)' % (', '.join(str(k + 1) for k in bad_fit), '; '.join('%.1f/%.2f' % fit[k] for k in bad_fit))
        else: status = 'verified' if n['not reproduced'] == 0 else ('partly verified' if n['not reproduced'] <= 2 else 'unverified'); note = ''
        if number == 121: note = (note + '; ' if note else '') + 'four of the eight candidates (the table of 121e-h is not usable)'
        cases.append(dict(case=name, candidates=len(cands), status=status, inputs=used, classes_changed=' '.join('%s->%s' % (lab(i), 'sp2' if sp2[i] else 'sp3') for i in flips), exchangeable_pairs=' '.join('%s/%s' % (lab(i), lab(j)) for i, j in pairs),
                          reproduced=n['reproduced'], consistent=n['consistent'], not_reproduced=n['not reproduced'], largest_difference=round(max(abs(x[3] - x[2]) for x in res), 1), note=note))
        for (code, key, val, p, st, rg), (c2, plab, _) in zip(res, [l for l in LEV for _ in FORM]):
            rec.append(dict(case=name, level=plab, formulation='ScTMS' if key[0] == 'sDP4+' else 'ScTMS+Uns', published=val, recomputed=round(p, 2), status=st, case_status=status))
        for code, plab, _ in LEV:
            s = [M.exchange(x, exp, pairs) for x in sh[code]] if pairs else sh[code]
            p = 100 * M.probabilities(D, s, exp, nuc, sp2, pars[code])[('DP4+', 'All')]
            for k in range(len(cands)): cand_rows.append(dict(case=name, level=plab, candidate=k + 1, is_correct=int(k == true), probability=round(float(p[k]), 4), case_status=status))
def write(name, rows_, keys=None):
    keys = keys or list(rows_[0].keys())
    w = csv.DictWriter(open(os.path.join(out_dir, name), 'w', newline=''), fieldnames=keys, extrasaction='ignore'); w.writeheader(); w.writerows(rows_)
write('jnp_s9_cases.csv', cases); write('jnp_s9_records.csv', rec); write('jnp_s9_candidates.csv', cand_rows)
print('cases of compounds 101-121: %d | %s' % (len(cases), dict(collections.Counter(c['status'] for c in cases))))
n = collections.Counter(r['status'] for r in rec if r['case_status'] != 'excluded')
print('printed values of the cases that are not excluded (two formulations, three levels): %d | reproduced %d | consistent within rounding %d | not reproduced %d' % (sum(n.values()), n['reproduced'], n['consistent'], n['not reproduced']))
for c in cases: print('  %-5s K=%d %-16s reproduced %d, consistent %d, not %d | largest difference %5.1f | %s%s%s%s' % (c['case'], c['candidates'], c['status'], c['reproduced'], c['consistent'], c['not_reproduced'], c['largest_difference'], c['inputs'], (' | classes: ' + c['classes_changed']) if c['classes_changed'] else '', (' | pairs: ' + c['exchangeable_pairs']) if c['exchangeable_pairs'] else '', (' | ' + c['note']) if c['note'] else ''))
