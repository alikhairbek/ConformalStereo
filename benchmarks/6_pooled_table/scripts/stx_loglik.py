"""STX150 validation set without the authors' files, version 2: which (case, level) records can be trusted?

A record is the experimental shifts of a case, the candidates printed in the SI and their shieldings at one level.
Its probabilities are recomputed for the printed candidates with the parameters of the paper.
  verified    the recomputed probability of the correct isomer reproduces the published one (within 1 point), or
              the published one lies inside the range left by the rounding of the printed shieldings (verify.py);
              either with the inputs by rule, or with unprinted inputs recovered from the six published values of the
              case (classes, exchangeable pairs; infer.py), accepted only if every usable level of the case is then
              reproduced with at most three changes
  corrected   127a and 137a: isomers exchanged between the (MMFF, xTB) and the B3LYP-geometry columns; the SI of the
              2025 paper fixes which block is the correct isomer; the probabilities are relabelled
              reproduced, or within 3 points, with at most three changes
  close       outside that range but within 10 points: same call, small numerical difference whose cause is not known
  unverified  more than 10 points away and outside that range (candidates missing from the SI, or no explanation)
  excluded    no experimental data printed; a copied shielding column at this level; two candidates with the same
              column; a candidate whose fit is not credible at this level (13C error above 6 ppm or 1H above
              0.6 ppm, and more than 3 times the error of the correct isomer); 154a (exchange not resolvable)
Writes stx_records.csv (one row per case and level) and stx_candidates.csv (one row per candidate).
Usage: stx_verified.py <lib dir> <scripts dir> <stx out dir> <stx_internal_consistency.csv> <stx_candidate_fit.csv> <out dir> [draws]"""
import sys, os, csv, collections
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[2])
import numpy as np
import dp4plus as D
import stx_data as S
import mmdp4 as M
import infer as I
stx_out, ic_csv, fit_csv, out_dir = sys.argv[3:7]; draws = int(sys.argv[7]) if len(sys.argv) > 7 else 200
blocks, cases, par = S.load(stx_out)
pub = {r['case']: r for r in csv.DictReader(open(os.path.join(stx_out, 'stx110_pcorrect.csv')))}
LEVELS = S.LEVELS; KEY = ('DP4+', 'All')
raw = collections.OrderedDict()
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_val_blocks.csv'))): raw.setdefault(r['block'], []).append(r)
copied = collections.defaultdict(set)
for name, rows in raw.items():
    for i, l1 in enumerate(LEVELS):
        for l2 in LEVELS[i + 1:]:
            a, b = [r[l1] for r in rows], [r[l2] for r in rows]
            if len(a) > 5 and sum(1 for x, y in zip(a, b) if x == y and x != '') >= 0.8 * len(a): copied[name].add(l1)
def shared_levels(cand):
    out = set()
    for i, a in enumerate(cand):
        for b in cand[i + 1:]:
            for l in LEVELS:
                va, vb = [r[l] for r in raw[a]], [r[l] for r in raw[b]]
                if len(va) == len(vb) and len(va) > 5 and sum(1 for x, y in zip(va, vb) if x == y and x != '') >= 0.8 * len(va): out.add(l)
    return out
exch = set()
for r in csv.DictReader(open(ic_csv)):
    if r['r_C'] != '' and r['r_H'] != '' and float(r['r_C']) < -0.8 and float(r['r_H']) < -0.8: exch.add(frozenset(r['pair'].split(' / ')))
fitv = collections.defaultdict(dict)
for r in csv.DictReader(open(fit_csv)):
    fitv[(r['case'], r['level'])][r['candidate']] = (float(r['rms_C']) if r['rms_C'] != '' else None, float(r['rms_H']) if r['rms_H'] != '' else None)
SECOND_IS_CORRECT = {'127a': ['MMFF_PBE', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9'], '137a': ['B3L_PBE', 'B3L_wB9']}
SOURCE_K = {'77a': 4, '77b': 4, '77d': 4, '106a': 4, '114i': 4, '116a': 4, '119c': 4, '119d': 4, '121a': 8, '141a': 4, '146a': 4, '149a': 4}

def probs(case, cand, lev, flip, pairs, noise=None):
    """probabilities (per cent) of the candidates; flip: set of row indices whose class is inverted; pairs: row pairs"""
    keep, nuc, dexp, sp2, sh = S.case_arrays(blocks, case, cand, lev)
    if len(keep) < 3: return None, 0
    sp2 = np.array([(not s) if k in flip else s for k, s in zip(keep, sp2)], bool); idx = {k: n for n, k in enumerate(keep)}
    pp = [(idx[i], idx[j]) for i, j in pairs if i in idx and j in idx]; terms = []
    for c in cand:
        v = np.array(sh[c], float)
        if noise is not None: v = v + noise.uniform(-0.005, 0.005, len(v))
        terms.append(D.log_terms(np.array(M.exchange(v, dexp, pp)), dexp, nuc, sp2, par[lev], truncate_nu=False))
    return 100 * D.probabilities(terms, nuc)[KEY], len(keep)
def state(case, cand, lev, flip, pairs, k, published):
    """-> probability vector, status, range"""
    p, n = probs(case, cand, lev, flip, pairs)
    if p is None: return None, 'too few nuclei', None, n
    if published is None: return p, 'no published value', None, n
    if abs(float(p[k]) - published) <= 1.0: return p, 'reproduced', None, n
    rng = np.random.default_rng(17); vals = [float(probs(case, cand, lev, flip, pairs, noise=rng)[0][k]) for _ in range(draws)]
    lo, hi = min(vals), max(vals)
    return p, ('consistent' if lo - 0.5 <= published <= hi + 0.5 else 'not reproduced'), (lo, hi), n
rec = []; cand_rows = []; recovered = collections.OrderedDict()
for case, info in cases.items():
    cand = info['cands']; rows = blocks[case]['rows']; has_exp = any(r['exp'] is not None for r in rows)
    sh_lev = shared_levels(cand); is_exch = any(e <= set(cand) for e in exch)
    P = {lev: (float(pub[case][lev]) if pub[case][lev] != '' else None) for lev in LEVELS}
    base = dict(case=case, candidate_set=info['set'], K_printed=len(cand), K_source=SOURCE_K.get(case, len(cand)))
    if not has_exp:
        for lev in LEVELS: rec.append(dict(base, level=lev, published=P[lev], status='excluded', reason='no experimental data printed'))
        continue
    true = {lev: (1 - cand.index(case) if lev in SECOND_IS_CORRECT.get(case, ()) else cand.index(case)) for lev in LEVELS}
    defects = {}
    for lev in LEVELS:
        f = fitv.get((case, lev), {}); fc = f.get(case, (None, None)); bad_fit = []
        for c2, (rc, rh) in f.items():
            if c2 == case: continue
            if (rc is not None and rc > 6 and fc[0] and rc > 3 * fc[0]) or (rh is not None and rh > 0.6 and fc[1] and rh > 3 * fc[1]): bad_fit.append(c2)
        d = []
        if any(lev in copied.get(c2, ()) for c2 in cand): d.append('copied column')
        if lev in sh_lev: d.append('two candidates share a column')
        if bad_fit: d.append('candidate without a credible fit: ' + ' '.join(bad_fit))
        if is_exch and case not in SECOND_IS_CORRECT: d.append('isomers exchanged between levels, not resolved')
        defects[lev] = d
    usable = [lev for lev in LEVELS if not defects[lev]]
    # as printed, inputs by rule; published value compared with the block the table took as correct (the first one)
    st0 = {lev: state(case, cand, lev, set(), [], cand.index(case), P[lev]) for lev in LEVELS}
    flip, pairs, how = set(), [], 'by rule'
    if case not in SECOND_IS_CORRECT and any(st0[lev][1] == 'not reproduced' for lev in usable):
        drows = [dict(nuc=r['nuc'], label=r['label'], exp=r['exp']) for r in rows]
        okrows = [i for i, r in enumerate(drows) if r['exp'] is not None]
        sub = [drows[i] for i in okrows]
        cpairs = [(okrows[i], okrows[j]) for i, j in I.candidate_pairs(sub)]; fmoves = [[okrows[i] for i in mv] for mv in I.flip_moves(sub)]
        tg = [lev for lev in usable if P[lev] is not None]
        def total(fl, pr): return sum(abs(float(probs(case, cand, lev, fl, pr)[0][cand.index(case)]) - P[lev]) for lev in tg)
        cur = total(flip, pairs); f2, p2 = set(), []
        for _ in range(3):
            best = None; used = {x for q in p2 for x in q}
            for mv in fmoves:
                v = total(f2 ^ set(mv), p2)
                if best is None or v < best[0]: best = (v, 'flip', mv)
            for q in cpairs:
                if q in p2 or q[0] in used or q[1] in used: continue
                v = total(f2, p2 + [q])
                if best is None or v < best[0]: best = (v, 'pair', q)
            if best is None or cur - best[0] < 0.5: break
            if best[1] == 'flip': f2 = f2 ^ set(best[2])
            else: p2 = p2 + [best[2]]
            cur = best[0]
            st2 = {lev: state(case, cand, lev, f2, p2, cand.index(case), P[lev]) for lev in usable}
            if all(st2[lev][1] in ('reproduced', 'consistent', 'no published value') or abs(float(st2[lev][0][cand.index(case)]) - P[lev]) <= 3.0 for lev in usable):
                flip, pairs, how = f2, p2, 'recovered'
                lab = lambda i: rows[i]['nuc'] + rows[i]['label']
                recovered[case] = ' '.join(['%s class changed' % lab(mv) for mv in sorted(f2) if rows[mv]['nuc'] == 'C'] + ['%s/%s exchangeable (%s, %s ppm)' % (lab(i), lab(j), rows[i]['exp'], rows[j]['exp']) for i, j in p2])
                break
    for lev in LEVELS:
        p0, s0, rg0, n0 = st0[lev]
        row = dict(base, level=lev, published=P[lev], n_nuclei=n0)
        if p0 is None: row.update(status='excluded', reason='too few nuclei'); rec.append(row); continue
        row['recomputed_as_printed'] = round(float(p0[cand.index(case)]), 2)
        if how == 'recovered': p, s, rg, _ = state(case, cand, lev, flip, pairs, cand.index(case), P[lev])
        else: p, s, rg = p0, s0, rg0
        k = true[lev]; row['p_correct'] = round(float(p[k]), 2); row['inputs'] = how if how == 'by rule' else 'recovered: ' + recovered[case]
        diff = None if P[lev] is None else float(p0[cand.index(case)]) - P[lev]; row['difference'] = None if diff is None else round(diff, 1)
        if rg: row['range_low'], row['range_high'] = round(rg[0], 1), round(rg[1], 1)
        reasons = list(defects[lev])
        if reasons: status = 'excluded'
        elif case in SECOND_IS_CORRECT: status = 'corrected'; reasons.append('isomers exchanged between levels; relabelled from the 2025 SI')
        elif s in ('reproduced', 'consistent'): status = 'verified'; reasons.append('reproduced' if s == 'reproduced' else 'inside the rounding range')
        elif P[lev] is not None and abs(float(p[cand.index(case)]) - P[lev]) <= 10.0: status = 'close'; reasons.append('within 10 points, outside the rounding range')
        else:
            status = 'unverified'
            reasons.append('fewer candidates printed than the source paper used' if (case in SOURCE_K and diff is not None and diff > 0) else 'recomputed and published differ, no explanation')
        row.update(status=status, reason='; '.join(reasons)); rec.append(row)
        for kk, c2 in enumerate(cand): cand_rows.append(dict(case=case, level=lev, candidate=c2, is_correct=int(kk == k), probability=round(float(p[kk]), 4), status=status,
                                                           candidate_set=base['candidate_set'], n_nuclei=n0,
                                                           log_ratio=round(float(np.log(max(float(p[kk]), 1e-300)) - np.log(max(float(max(p)), 1e-300))), 5)))
keys = ['case', 'level', 'candidate_set', 'K_printed', 'K_source', 'status', 'reason', 'inputs', 'published', 'recomputed_as_printed', 'p_correct', 'difference', 'range_low', 'range_high', 'n_nuclei']
w = csv.DictWriter(open(os.path.join(out_dir, 'stx_records.csv'), 'w', newline=''), fieldnames=keys); w.writeheader(); w.writerows(rec)
w = csv.DictWriter(open(os.path.join(out_dir, 'stx_candidates.csv'), 'w', newline=''), fieldnames=list(cand_rows[0].keys())); w.writeheader(); w.writerows(cand_rows)
ST = ['verified', 'close', 'corrected', 'unverified', 'excluded']
print('records (case, level): %d | rounding draws %d' % (len(rec), draws))
print('%-9s %s' % ('level', ' '.join('%10s' % s for s in ST)))
for lev in LEVELS:
    c = collections.Counter(r['status'] for r in rec if r['level'] == lev); print('%-9s %s' % (lev, ' '.join('%10d' % c[s] for s in ST)))
c = collections.Counter(r['status'] for r in rec); print('%-9s %s' % ('all', ' '.join('%10d' % c[s] for s in ST)))
v = [r for r in rec if r['status'] == 'verified']
print('verified: reproduced within 1 point %d | inside the rounding range %d | of these with recovered inputs %d records in %d cases' % (sum('reproduced' in r['reason'] for r in v), sum('rounding' in r['reason'] for r in v), sum(r['inputs'] != 'by rule' for r in v), len(recovered)))
for k, t in recovered.items(): print('   inputs recovered for %-5s %s' % (k, t))
full = collections.Counter(tuple(sorted({r['status'] for r in rec if r['case'] == case})) for case in cases)
good = sum(1 for case in cases if all(r['status'] in ('verified', 'close', 'corrected') for r in rec if r['case'] == case))
print('cases verified at all six levels: %d | verified, close or corrected at all six: %d | excluded at every level: %d' % (full[('verified',)], good, full[('excluded',)]))
print('reasons for exclusion:', dict(collections.Counter(r['reason'].split(':')[0] for r in rec if r['status'] == 'excluded')))
print('unverified:', dict(collections.Counter(r['reason'] for r in rec if r['status'] == 'unverified')))
u = [r for r in rec if r['status'] == 'unverified' and 'no explanation' in r['reason']]
print('unverified without explanation: 10 to 50 points %d | above 50 %d | cases: %s' % (sum(abs(r['p_correct'] - r['published']) <= 50 for r in u), sum(abs(r['p_correct'] - r['published']) > 50 for r in u), ' '.join(sorted({r['case'] for r in u}))))
