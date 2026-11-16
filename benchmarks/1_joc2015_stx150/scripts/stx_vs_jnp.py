"""STX150 validation table against the group's 2023 paper: same cases, same levels (MMFF geometries with mPW1PW91 or
wB97XD shieldings, and DP4+ at the recommended DFT level). Which formulation does the STX150 table hold, and where do
the two publications disagree?
Usage: stx_vs_jnp.py <out_dir of this folder> <out dir of the SI extraction (stx110_pcorrect.csv, stx_validation_recomputed.csv)>"""
import sys, os, csv, collections
out, stx_out = sys.argv[1:3]
jnp = {r['case']: r for r in csv.DictReader(open(os.path.join(out, 'jnp_probabilities.csv')))}
stx = {r['case']: r for r in csv.DictReader(open(os.path.join(stx_out, 'stx110_pcorrect.csv')))}
common = [c for c in stx if c in jnp]
import re as _re
_n = lambda c: int(_re.match(r'\d+', c).group())
print('STX150 validation cases: %d | cases of the 2023 paper: %d | in both, by name: %d | STX150 cases numbered above 121: %d | numbered up to 121 under a name the 2023 SI does not use: %s' % (len(stx), len(jnp), len(common), sum(1 for c in stx if _n(c) > 121), ' '.join(c for c in stx if c not in jnp and _n(c) <= 121)))
FORMS = ['ScTMS', 'ScTMS+Uns', 'ScMSTD', 'ScMSTD+Uns', 'ScTMS+ScMSTD', 'ScTMS+ScMSTD+Uns']
best = {}
for col, lev in (('MMFF_mPW', 'MMFF_mPW'), ('MMFF_wB9', 'MMFF_wB9')):
    res = []
    for f in FORMS:
        d = [abs(float(stx[c][col]) - float(jnp[c]['%s|%s' % (lev, f)])) for c in common if jnp[c]['%s|%s' % (lev, f)] != '']
        res.append((sum(x <= 1.0 for x in d), len(d), f))
    print('STX150 column %-9s against the 2023 values, cases within 1 point: %s' % (col, ' | '.join('%s %d/%d' % (f, a, n) for a, n, f in res)))
    best[col] = max(res)[2]
d = [(c, float(stx[c]['B3L_mPW']), float(jnp[c]['DP4+'])) for c in common if jnp[c]['DP4+'] != '']
print('STX150 column B3L_mPW against DP4+ of the 2023 SI: within 1 point in %d of %d | differing: %s' % (sum(abs(a - b) <= 1 for _, a, b in d), len(d), ' '.join('%s %.0f/%.0f' % x for x in d if abs(x[1] - x[2]) > 1)))
print()
for col in ('MMFF_mPW', 'MMFF_wB9'):
    f = best[col]
    diff = [(c, float(stx[c][col]), float(jnp[c]['%s|%s' % (col, f)])) for c in common if abs(float(stx[c][col]) - float(jnp[c]['%s|%s' % (col, f)])) > 1.0]
    print('%s, formulation %s: %d of %d cases differ by more than 1 point (STX150 / 2023): %s' % (col, f, len(diff), len(common), ' '.join('%s %.0f/%.0f' % x for x in diff)))
rec = collections.defaultdict(dict)
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_validation_recomputed.csv'))):
    if r['published'] != '' and r['mine'] != '': rec[r['case']][r['level']] = (float(r['published']), float(r['mine']))
print('\ncontradicted cases (STX150 table says <= 1%, its printed shieldings say >= 95%) and what the 2023 paper gives for the same case:')
print('%-5s | STX150 table: MMFF_PBE MMFF_mPW MMFF_wB9 | xTB_PBE xTB_mPW xTB_wB9 | B3L_PBE B3L_mPW B3L_wB9 || 2023 paper: MMFF_B3L MMFF_mPW MMFF_wB9 | DP4+ (DFT)' % 'case')
for c in ('89a', '106a', '111a', '114i', '115a', '121a', '127a', '146a'):
    s = stx[c]; j = jnp.get(c)
    left = '%8s %8s %8s | %7s %7s %7s | %7s %7s %7s' % tuple(s[k] for k in ('MMFF_PBE', 'MMFF_mPW', 'MMFF_wB9', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9', 'B3L_PBE', 'B3L_mPW', 'B3L_wB9'))
    right = ('%8.0f %8.0f %8.0f | %s' % (float(j['MMFF_B3L|' + best['MMFF_mPW']]), float(j['MMFF_mPW|' + best['MMFF_mPW']]), float(j['MMFF_wB9|' + best['MMFF_wB9']]), j['DP4+'])) if j else 'not in the 2023 paper'
    mine = ' '.join('%s=%.0f' % (l, v[1]) for l, v in rec[c].items())
    print('%-5s | %s || %s' % (c, left, right)); print('      | recomputed from the printed STX150 shieldings: %s' % mine)
