"""2025 SI (Org. Biomol. Chem. 2025, 23, 2480): the three levels for which it prints probabilities, recomputed from its
printed shieldings with the parameters shipped in DP4+App and compared with the nine probabilities printed per compound
and level (H, C and all data; scaled, unscaled, full).
  Level A  SMD/wB97XD/6-31+G**//MMFF            MM-DP4+ parameters
  Level B  PCM/mPW1PW91/6-31+G**//B3LYP/6-31G*  DP4+ parameters
  Level C  PCM/mPW1PW91/6-31+G**//MMFF          MM-DP4+ parameters
For each compound and level: which printed block reproduces the nine values as the correct isomer (the first one, as
the SI states, or another: then the blocks are exchanged at that level); inputs by rule, then recovered from the printed
values if needed (infer.py), jointly for the levels of the compound; rounding of the printed shieldings (three
decimals) taken into account (verify.py). A level whose columns hold cells that cannot be the shielding of their
nucleus is 'misprinted' unless it is reproduced all the same.
Writes obc_levels.csv (compound, level), obc_candidates.csv (probability of every isomer).
Usage: obc_verify.py <lib dir> <scripts dir> <app_parameters.csv> <obc_compounds.json> <d4ob01959h1_suppl.txt> <out dir> [draws]"""
import sys, os, re, csv, json, collections
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[2])
import numpy as np
import dp4plus as D
import joc_data as J
import mmdp4 as M
import infer as I
import verify as V
MMp = M.load_parameters(sys.argv[3], 'MM'); QMp = M.load_parameters(sys.argv[3], 'QM')
PAR = {'A': MMp['wB97XD.6-31+G(d,p).SMD'], 'B': QMp['mPW1PW91.6-31+G(d,p).PCM'], 'C': MMp['mPW1PW91.6-31+G(d,p).PCM']}
C = json.load(open(sys.argv[4])); out_dir = sys.argv[6]; draws = int(sys.argv[7]) if len(sys.argv) > 7 else 100
lines = open(sys.argv[5], encoding='utf-8', errors='replace').read().split('\n')
ROW = re.compile(r'^\s*(\d{1,2})\s+((?:\d+\.\d+%\s+){8}\d+\.\d+%)\s*$'); sec = None; P = collections.defaultdict(dict)
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
levels_out = []; cand_rows = []
for num, c in C.items():
    n = int(num); isos = list(c['isomers']); rows = [dict(nuc=r['nucleus'], label=r['label'], exp=r['exp']) for r in c['rows']]
    exp = [r['exp'] for r in rows]; nuc = [r['nuc'] for r in rows]; lab = lambda i: rows[i]['nuc'] + rows[i]['label']
    sh = {L: [np.array(c['isomers'][s][L], float) for s in isos] for L in 'ABC'}
    bad = {L: sorted({i for iso, l, i in c['bad_cells'] if l == L}) for L in 'ABC'}
    sp2 = J.sp2_flags(rows, 100.0, None)
    def nine(L, sp2_, pairs_):
        s = [M.exchange(x, exp, pairs_) for x in sh[L]] if pairs_ else sh[L]
        return M.probabilities(D, s, exp, nuc, sp2_, PAR[L])
    # block of the correct isomer at each level, with the inputs by rule
    block = {}
    for L in 'ABC':
        pr = nine(L, sp2, []); dev = [sum(abs(100 * float(pr[o][k]) - v) for o, v in zip(ORDER, P[n][L]) if o in pr) for k in range(len(isos))]
        block[L] = int(np.argmin(dev))
    def targets(Ls): return [(L, o, block[L], v) for L in Ls for o, v in zip(ORDER, P[n][L])]
    def assess(sp2_, pairs_):
        out = {}
        for L in 'ABC':
            pr = nine(L, sp2_, pairs_); st = []
            for o, v in zip(ORDER, P[n][L]):
                if o not in pr: continue
                p = 100 * float(pr[o][block[L]])
                st.append((o, v, p, V.status(p, v, lambda: V.rounding_range(D, sh[L], exp, nuc, sp2_, PAR[L], o, block[L], 0.0005, pairs=pairs_, n=draws, seed=17), printed_half_unit=0.005)[0]))
            out[L] = st
        return out
    res = assess(sp2, []); nbad = sum(x[3] == 'not reproduced' for L in 'ABC' for x in res[L]); used = 'by rule'; flips = []; pairs = []
    if nbad:
        clean = [L for L in 'ABC' if not bad[L]] or list('ABC')
        inf = I.infer(D, rows, sh, sp2, targets(clean), PAR)
        res2 = assess(inf['sp2'], inf['pairs']); nbad2 = sum(x[3] == 'not reproduced' for L in 'ABC' for x in res2[L])
        if nbad2 < nbad: res, used, flips, pairs, sp2 = res2, 'recovered', inf['flips'], inf['pairs'], inf['sp2']
    for L in 'ABC':
        cnt = collections.Counter(x[3] for x in res[L]); tot = sum(cnt.values()); notr = cnt['not reproduced']
        if notr == 0: status = 'verified'
        elif bad[L]: status = 'misprinted'
        elif notr <= 2: status = 'partly verified'
        else: status = 'unverified'
        allfull = [x for x in res[L] if x[0] == ('DP4+', 'All')][0]
        levels_out.append(dict(compound=n, level=L, isomers=len(isos), pair_with=c['pair_with'] or '', status=status, block_of_correct_isomer='abcd'[block[L]], values=tot, reproduced=cnt['reproduced'], consistent=cnt['consistent'], not_reproduced=notr,
                               printed_all_full=allfull[1], recomputed_all_full=round(allfull[2], 2), misprinted_cells=len([1 for iso, l, i in c['bad_cells'] if l == L]), inputs=used,
                               classes_changed=' '.join('%s->%s' % (lab(i), 'sp2' if sp2[i] else 'sp3') for i in flips), exchangeable_pairs=' '.join('%s/%s' % (lab(i), lab(j)) for i, j in pairs)))
        pr = nine(L, sp2, pairs)
        for k, s in enumerate(isos): cand_rows.append(dict(compound=n, level=L, isomer=s, is_correct=int(k == block[L]), probability=round(100 * float(pr[('DP4+', 'All')][k]), 4), level_status=status))
def write(name, rows_):
    w = csv.DictWriter(open(os.path.join(out_dir, name), 'w', newline=''), fieldnames=list(rows_[0].keys())); w.writeheader(); w.writerows(rows_)
write('obc_levels.csv', levels_out); write('obc_candidates.csv', cand_rows)
print('compounds: %d | records (compound, level): %d' % (len(C), len(levels_out)))
print('%-7s %s' % ('level', ' '.join('%16s' % s for s in ('verified', 'partly verified', 'unverified', 'misprinted'))))
for L in 'ABC':
    c_ = collections.Counter(r['status'] for r in levels_out if r['level'] == L); print('%-7s %s' % (L, ' '.join('%16d' % c_[s] for s in ('verified', 'partly verified', 'unverified', 'misprinted'))))
c_ = collections.Counter(r['status'] for r in levels_out); print('%-7s %s' % ('all', ' '.join('%16d' % c_[s] for s in ('verified', 'partly verified', 'unverified', 'misprinted'))))
tot = sum(r['values'] for r in levels_out); print('printed values: %d | reproduced %d | consistent within rounding %d | not reproduced %d' % (tot, sum(r['reproduced'] for r in levels_out), sum(r['consistent'] for r in levels_out), sum(r['not_reproduced'] for r in levels_out)))
print('compounds verified at all three levels: %d | at two: %d' % (sum(1 for n in C if all(r['status'] == 'verified' for r in levels_out if r['compound'] == int(n))), sum(1 for n in C if sum(r['status'] == 'verified' for r in levels_out if r['compound'] == int(n)) == 2)))
print('blocks exchanged (the correct isomer is not the first printed block):', [(r['compound'], r['level'], r['block_of_correct_isomer'], r['status']) for r in levels_out if r['block_of_correct_isomer'] != 'a'])
rc = collections.OrderedDict()
for r in levels_out:
    if r['inputs'] == 'recovered': rc[r['compound']] = (r['classes_changed'], r['exchangeable_pairs'])
print('inputs recovered for %d compounds:' % len(rc), '; '.join('%d: %s %s' % (k, v[0], v[1]) for k, v in rc.items()))
for r in levels_out:
    if r['status'] != 'verified': print('  compound %2d level %s %-16s values not reproduced %d of %d | all data: printed %.2f recomputed %.2f | misprinted cells %d' % (r['compound'], r['level'], r['status'], r['not_reproduced'], r['values'], r['printed_all_full'], r['recomputed_all_full'], r['misprinted_cells']))
