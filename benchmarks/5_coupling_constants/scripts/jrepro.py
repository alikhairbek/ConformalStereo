"""Reproduce the coupling part of the depositors' PyDP4 outputs from the raw coupling matrices.
Usage: jrepro.py <res/core run 2> <res/core run 3> <lib> <dep17_experiments.json> <jsource.json> <out folder>"""
import sys, os, json, math, collections, csv, itertools
core2, core3, lib, exp_json, src_json, out = sys.argv[1:7]
sys.path.insert(0, lib)
import numpy as np
import dep17, desc17, parse17, raw17, core17, jlib
from energy_source17 import renum_maps

S = dep17.Small(core2)
E = {(e['archive'], e['name']): e for e in json.load(open(exp_json))}
SRC = json.load(open(src_json))
data = raw17.load(core2, {r['archive'] for r in SRC})
M = jlib.load_matrices(core3)
print('conformers with both coupling matrices:', len(M))


def candidate(rec, desc, maps, stem, si, use_map, kind, last_labels=None):
    """Everything PyDP4 had for one candidate: assignment of the proton shifts and the pruned coupling matrix."""
    arch, base, folder = rec['archive'], rec['base'], rec['folder']
    lev = data[(arch, rec['jbase'], stem)]['J']
    ks = sorted(lev)
    pops = jlib.populations([lev[k][0] for k in ks])
    n = len(lev[ks[0]][1])
    sig = [sum(p * lev[k][1][i] for p, k in zip(pops, ks)) for i in range(n)]
    els = lev[ks[0]][2]
    if use_map:
        sig, els = raw17.renumber(sig, maps[si - 2]), raw17.renumber(els, maps[si - 2])
    at = core17.atoms_of(desc, els)
    labels = [lab for lab, idx in at['H']]
    shield = [sum(sig[i] for i in idx) / len(idx) for lab, idx in at['H']]
    pairs = core17.assign(labels, [-x for x in shield], [(s, l) for s, l, j in desc['H']])
    first_j = {}
    for s, l, j in desc['H']:
        first_j.setdefault(s, j if j else [0.0])              # PyDP4 takes the couplings of the first entry with that shift
    exp_j = [first_j[x] for lab, c, x in pairs]
    mats = [M.get((arch, folder, stem, k)) for k in ks]
    if any(m is None for m in mats):
        return None
    avg = sum(p * m[kind] for p, m in zip(pops, mats))
    # PyDP4 kept one list of hydrogen numbers for all candidates, that of the last candidate, and used it by position:
    # where the candidates are numbered differently (experiments with a renumbering file) rows end up under wrong labels
    hl = last_labels if last_labels is not None else [str(h) for h in mats[0]['h']]
    pruned, pl = jlib.zero_equiv(avg, hl, desc['groups'], desc['omit'])
    rows = [list(pruned[pl.index(lab[1:])]) if lab[1:] in pl else [0.0] for lab, c, x in pairs]
    exp_used, assigned = jlib.assign_j(exp_j, rows)
    return dict(labels=[p[0] for p in pairs], exp=[p[2] for p in pairs], exp_j=exp_used, calc_j=assigned, n_conf=len(ks), rows=rows, pops=pops,
                own_labels=[str(h) for h in mats[0]['h']])


fmt = lambda v: ['%4.1f' % x for x in v]
stat = collections.Counter()
rows_out = []
blocks_out, ties, deg = [], [], []
misplaced = set()
ratio = collections.defaultdict(list)
for rec in SRC:
    d = parse17.parse(S.read(rec['path']))
    desc = desc17.parse(S.read(rec['path'][:-4]).splitlines())
    maps = renum_maps(S, rec['archive'], rec['folder'])
    last = [f for f in rec['isomers'][-1]['fits'] if f[2] == 'J' and f[3] == 'J']
    last_labels = None
    if last:
        lv = data[(rec['archive'], rec['jbase'], last[0][0])]['J']
        m0 = M.get((rec['archive'], rec['folder'], last[0][0], sorted(lv)[0]))
        last_labels = [str(h) for h in m0['h']] if m0 else None
    for it in rec['isomers']:
        k = it['isomer']
        own = [f for f in it['fits'] if f[2] == 'J' and f[3] == 'J']
        if not own:
            stat['isomer blocks without identified structure'] += 1
            continue
        stem, si, sl, el, use_map, nconf = own[0]
        block = d['isomers'][k]
        res = {}
        for kind in ('total', 'FC'):
            c = candidate(rec, desc, maps, stem, si, use_map, kind, last_labels)
            if c is None:
                stat['isomer blocks with a conformer matrix missing'] += 1
                break
            res[kind] = c
        if len(res) < 2:
            continue
        stat['isomer blocks compared'] += 1
        # the same with the candidate's own hydrogen numbers: does the printed block need the numbers of the last candidate?
        own = candidate(rec, desc, maps, stem, si, use_map, 'total', None)
        pos_ok = all((['%4.1f' % x for x in a if abs(x) > jlib.J_THRESH] if e != [0.0] else []) == ['%4.1f' % x for x in r[6]]
                     for r, e, a in zip(block['H'], res['total']['exp_j'], res['total']['calc_j']))
        own_ok = all((['%4.1f' % x for x in a if abs(x) > jlib.J_THRESH] if e != [0.0] else []) == ['%4.1f' % x for x in r[6]]
                     for r, e, a in zip(block['H'], own['exp_j'], own['calc_j']))
        if pos_ok and not own_ok:
            stat['isomer blocks whose printed couplings follow only from the atom numbers of the last candidate'] += 1
            misplaced.add(rec['file'])
        printed = block['H']
        c = res['total']
        # matchings with the same (smallest) summed absolute error but different pairs: the choice between them is arbitrary
        par = jlib.PARAMS['total']
        for e, row in zip(c['exp_j'], c['rows']):
            if e == [0.0]:
                continue
            cj = [abs(j) for j in row if abs(j) > jlib.J_THRESH]
            best, group = None, {}
            for perm in itertools.permutations(cj, len(e)):
                err = sum(abs(perm[i] - e[i]) for i in range(len(e)))
                if best is None or err < best - 1e-9:
                    best, group = err, {}
                if abs(err - best) <= 1e-9:
                    key = tuple('%.6f' % x for x in perm)
                    group.setdefault(key, jlib.log_gauss(jlib.coupling_errors([e], [list(perm)], par['slope'], par['intercept']), par['mean'], par['sd']))
            stat['protons with experimental couplings matched'] += 1
            if len(group) > 1:
                stat['protons with several equally good matchings'] += 1
                deg.append(max(group.values()) - min(group.values()))
        lab_ok = [r[0] for r in printed] == c['labels']
        stat['proton assignment as printed' if lab_ok else 'proton assignment DIFFERENT'] += 1
        exp_ok = all(fmt(r[5]) == (fmt(e) if e != [0.0] else []) for r, e in zip(printed, c['exp_j'])) and len(printed) == len(c['exp_j'])
        stat['experimental couplings as printed' if exp_ok else 'experimental couplings DIFFERENT'] += 1
        for kind in ('total', 'FC'):
            cc = res[kind]
            same = n = 0
            for r, e, a, row in zip(printed, cc['exp_j'], cc['calc_j'], cc['rows']):
                if e != [0.0] or r[5]:
                    n += 1
                    mine = ['%4.1f' % x for x in a if abs(x) > jlib.J_THRESH] if e != [0.0] else []
                    same += mine == fmt(r[6])
                    if kind == 'total' and mine != fmt(r[6]):
                        # the same couplings in another order, with the same summed error: a tie in the matching
                        tie = sorted(mine) == sorted(fmt(r[6])) and abs(sum(abs(float(x) - y) for x, y in zip(mine, e)) - sum(abs(x - y) for x, y in zip(r[6], e))) < 0.11
                        ties.append((rec['file'], k, r[0], fmt(e), fmt(r[6]), mine, tie))
            stat['protons with couplings (%s)' % kind] += n
            stat['calculated couplings as printed (%s)' % kind] += same
            if same == n:
                stat['isomer blocks with every calculated coupling as printed (%s)' % kind] += 1
            par = jlib.PARAMS[kind]
            errs = jlib.coupling_errors(cc['exp_j'], cc['calc_j'], par['slope'], par['intercept'])
            lg = jlib.log_gauss(errs, par['mean'], par['sd'])
            pj = block['p']['J']
            if pj > 0:
                ratio[kind].append((lg - math.log(pj), rec['file'], k, len(errs)))
            rows_out.append([rec['file'], k, stem, kind, n, same, len(errs), '%.4g' % math.exp(lg) if lg > -700 else '0', '%.3g' % pj])
            if kind == 'total':
                blocks_out.append(dict(file=rec['file'], isomer=k, structure=stem, all_as_printed=(same == n), printed=pj, log_gauss=lg,
                                       n_without=sum(1 for e in cc['exp_j'] if e == [0.0]),
                                       pairs=[[e_, a_] for e, a in zip(cc['exp_j'], cc['calc_j']) if e != [0.0] for e_, a_ in zip(e, a)]))
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'jrepro_blocks.csv'), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['output', 'isomer', 'structure', 'matrix', 'protons_with_couplings', 'protons_with_calculated_couplings_as_printed', 'couplings', 'my_J_term', 'printed_J_term'])
    w.writerows(rows_out)
json.dump(blocks_out, open(os.path.join(out, 'jrepro_couplings.json'), 'w'))
for k, v in sorted(stat.items()):
    print('%6d  %s' % (v, k))
print('protons whose calculated couplings differ from the printed ones: %d, of which ties in the matching (same values, other order, same summed error): %d' % (
    len(ties), sum(1 for t in ties if t[6])))
for f, k, lab, e, pr, mine, tie in ties:
    print('   %-20s isomer %2d %-4s exp %s | printed %s | here %s%s' % (f, k, lab, ','.join(e), ','.join(pr), ','.join(mine), '' if tie else '   NOT A TIE'))
print('outputs with such blocks: %d (%s)' % (len(misplaced), ', '.join(sorted(misplaced))))
deg.sort()
print('where several matchings are equally good, the Gaussian coupling term of the source changes with the choice by a factor of e^x: median x %.1f, largest %.1f, above 2.3 (a factor of ten) in %d protons' % (
    deg[len(deg) // 2], deg[-1], sum(1 for x in deg if x > 2.3)))
for kind in ('total', 'FC'):
    v = ratio[kind]
    close = sum(1 for x in v if abs(x[0]) < 0.006)
    files = collections.Counter(x[1] for x in v if abs(x[0]) < 0.006)
    print('coupling term as the Gaussian of the PyDP4 source with the %s matrix: equal to the printed value (three digits) in %d of %d blocks%s' % (
        kind, close, len(v), (' (' + ', '.join('%s %d' % kv for kv in sorted(files.items())) + ')') if files else ''))
