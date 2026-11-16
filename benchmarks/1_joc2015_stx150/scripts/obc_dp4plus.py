"""DP4+ at Level B (PCM/mPW1PW91/6-31+G**//B3LYP/6-31G*) for the 39 compounds of the 2025 paper, recomputed from the
shieldings printed in its SI (obc_tables.py, version 2) and compared with the nine probabilities printed for the
correct isomer (SI section 4: H, C and all data, each scaled, unscaled and full).
The classes of the nuclei are not printed: carbons at or above <threshold> ppm are sp2, protons follow the carbon
with the same label and, when there is none, are sp2 from 5 ppm (the rule used for the other sets). Writes every candidate's probability.
Usage: obc_dp4plus.py <lib dir> <joc_dp4plus_parameters.csv> <obc_compounds.json> <d4ob01959h1_suppl.txt> <out dir> [threshold]"""
import sys, os, re, csv, json, collections
sys.path.insert(0, sys.argv[1])
import numpy as np
import dp4plus as D
import joc_data as J
par = {r['level']: {k: float(v) for k, v in r.items() if k != 'level'} for r in csv.DictReader(open(sys.argv[2]))}['mPW1PW91/6-31+G**/PCM']
C = json.load(open(sys.argv[3])); out_dir = sys.argv[5]; thr = float(sys.argv[6]) if len(sys.argv) > 6 else 100.0
lines = open(sys.argv[4], encoding='utf-8', errors='replace').read().split('\n')
ROW = re.compile(r'^\s*(\d{1,2})\s+((?:\d+\.\d+%\s+){8}\d+\.\d+%)\s*$')
sec = None; P = collections.defaultdict(dict)
for i, l in enumerate(lines):
    if i < 60: continue
    s = l.strip()
    if s.startswith('3. MM-DP4+ probabilities'): sec = 'A'
    elif s.startswith('4. DP4+ probabilities'): sec = 'B'
    elif s.startswith('5. MM-DP4+ probabilities'): sec = 'C'
    elif s.startswith('6.'): sec = None
    m = ROW.match(l)
    if m and sec: P[int(m.group(1))][sec] = [float(x.rstrip('%')) for x in m.group(2).split()]
ORDER = [('sDP4+', 'H'), ('sDP4+', 'C'), ('sDP4+', 'All'), ('uDP4+', 'H'), ('uDP4+', 'C'), ('uDP4+', 'All'), ('DP4+', 'H'), ('DP4+', 'C'), ('DP4+', 'All')]
NAMES = ['H sca', 'C sca', 'All sca', 'H uns', 'C uns', 'All uns', 'H full', 'C full', 'All full']
def classes(rows):
    """same rule as for the 2015 and STX150 sets (joc_data.sp2_flags)"""
    return J.sp2_flags([dict(nuc=r['nucleus'], label=r['label'], exp=r['exp']) for r in rows], thr, None)
def run(c, isos, trunc, sp2, keep):
    rows = [c['rows'][i] for i in keep]; nuc = np.array([r['nucleus'] for r in rows]); dexp = np.array([r['exp'] for r in rows])
    terms = [D.log_terms(np.array([c['isomers'][s]['B'][i] for i in keep]), dexp, nuc, np.array([sp2[i] for i in keep], bool), par, truncate_nu=trunc) for s in isos]
    pr = D.probabilities(terms, nuc)
    return pr, {s: [100 * pr[o][k] if o in pr else float('nan') for o in ORDER] for k, s in enumerate(isos)}
res = []; cand_rows = []; summ = collections.Counter()
for n, c in C.items():
    n = int(n); pub = P[n]['B']; isos = list(c['isomers']); sp2 = classes(c['rows'])
    badB = sorted({i for iso, l, i in c['bad_cells'] if l == 'B'})
    best = None
    for rows_used, keep in (('all rows', list(range(len(c['rows'])))), ('rows with a misprinted cell left out', [i for i in range(len(c['rows'])) if i not in badB])):
        if rows_used != 'all rows' and not badB: continue
        for trunc in (True, False):
            pr, nine = run(c, isos, trunc, sp2, keep)
            for s in isos:
                dev = max(abs(a - b) for a, b in zip(nine[s], pub) if a == a)
                key = (0, rows_used != 'all rows', s != 'a', dev) if dev <= 1 else (1, dev)   # among reproductions prefer all rows and isomer a; else the closest
                if best is None or key < best[0]: best = (key, dev, s, trunc, rows_used, pr, nine, keep)
    key, dev, s, trunc, rows_used, pr, nine, keep = best
    mine = nine[s]; worst = max(range(9), key=lambda i: abs(mine[i] - pub[i]) if mine[i] == mine[i] else -1)
    status = 'reproduced' if dev <= 1 else ('close' if dev <= 10 else 'not reproduced')
    res.append(dict(compound=n, isomers=len(isos), n_C=sum(c['rows'][i]['nucleus'] == 'C' for i in keep), n_H=sum(c['rows'][i]['nucleus'] == 'H' for i in keep), status=status, block_of_correct_isomer=s, rows_used=rows_used,
                    published_all_full=pub[8], recomputed_all_full=round(mine[8], 2), largest_difference=round(dev, 2), where=NAMES[worst], integer_nu=int(trunc),
                    isomer_a_all_full=round(nine['a'][8], 2), misprinted_cells_level_B=len([1 for iso, l, i in c['bad_cells'] if l == 'B'])))
    summ[status] += 1
    for k, s2 in enumerate(isos):
        cand_rows.append(dict(compound=n, isomer=s2, is_correct=int(s2 == s), status=status, **{('%s_%s' % o).replace('+', 'plus'): round(100 * pr[o][k], 4) for o in ORDER if o in pr}))
w = csv.DictWriter(open(os.path.join(out_dir, 'obc_dp4plus_levelB.csv'), 'w', newline=''), fieldnames=list(res[0].keys())); w.writeheader(); w.writerows(res)
keys = []
for r in cand_rows:
    for k in r:
        if k not in keys: keys.append(k)
w = csv.DictWriter(open(os.path.join(out_dir, 'obc_dp4plus_levelB_candidates.csv'), 'w', newline=''), fieldnames=keys); w.writeheader(); w.writerows(cand_rows)
print('compounds: %d | sp2 threshold %.0f ppm | nine printed probabilities: %s' % (len(res), thr, dict(summ)))
print('reproduced (all nine within 1 point): with the first printed block as the correct isomer %d | with another block %d %s | after leaving out rows with a misprinted cell %d' % (
    sum(r['status'] == 'reproduced' and r['block_of_correct_isomer'] == 'a' and r['rows_used'] == 'all rows' for r in res), sum(r['status'] == 'reproduced' and r['block_of_correct_isomer'] != 'a' for r in res),
    [(r['compound'], r['block_of_correct_isomer']) for r in res if r['status'] == 'reproduced' and r['block_of_correct_isomer'] != 'a'], sum(r['status'] == 'reproduced' and r['rows_used'] != 'all rows' for r in res)))
print('degrees of freedom as integers (the Excel sheet) reproduces better in %d compounds, as printed in %d' % (sum(r['integer_nu'] for r in res), sum(1 - r['integer_nu'] for r in res)))
for r in res:
    if r['status'] != 'reproduced': print('  compound %2d (%d isomers, %d C, %d H): %s | best block %s, %s | all data: published %.2f, best %.2f, first block %.2f | largest difference %.2f at %s | misprinted cells at Level B: %d' % (
        r['compound'], r['isomers'], r['n_C'], r['n_H'], r['status'], r['block_of_correct_isomer'], r['rows_used'], r['published_all_full'], r['recomputed_all_full'], r['isomer_a_all_full'], r['largest_difference'], r['where'], r['misprinted_cells_level_B']))
print('compounds with four isomers, probability of the correct isomer (all data) with all four and with the first two only (what the STX150 SI prints for 141a, 146a and 149a):')
for n, c in C.items():
    if len(c['isomers']) == 4:
        r = next(x for x in res if x['compound'] == int(n)); keep = list(range(len(c['rows']))); sp2 = classes(c['rows'])
        four = run(c, list('abcd'), bool(r['integer_nu']), sp2, keep)[1]['a'][8]; two = run(c, list('ab'), bool(r['integer_nu']), sp2, keep)[1]['a'][8]
        print('  compound %2s: published %.2f | four isomers %.2f | first two isomers only %.2f' % (n, P[int(n)]['B'][8], four, two))
