"""Shift and coupling likelihoods of every candidate, recomputed from the raw outputs of the coupling folders.
Usage: jbench.py <res/core run 2> <res/core run 3> <lib> <jsource.json> <dep17_truth.csv> <jmodel_recovered.json> <out folder> [draws] [seed]

One record per compound with experimental couplings (33) and candidate. Shieldings, energies and coupling matrices are
those of the coupling jobs (B3LYP/6-31G(d,p) on MMFF geometries, second job step). The coupling term is defined so that
it does not depend on the order or numbering of the atoms, which the term of the deposit outputs does:
  * candidates whose atoms are numbered differently are brought to the numbering of the first candidate for the couplings
    as well (PyDP4 did this for the shieldings only);
  * every experimental signal keeps its own couplings (PyDP4 gave signals with equal shifts the couplings of the first);
  * calculated couplings are matched to the experimental ones so that the likelihood under the scoring model is largest
    (PyDP4 minimised the summed absolute error, which is the same for every matching whenever all calculated values lie
    on one side of the experimental ones, and then took the first matching in the order of the atoms);
  * an experimental coupling without a calculated partner is compared with zero (PyDP4 skipped a proton without
    calculated couplings, and set all calculated values of a proton to zero if they were fewer than the experimental).
Kept from PyDP4: couplings below 0.2 Hz are not partners; equivalent protons are averaged and their mutual couplings set
to zero; the linear correction of the calculated couplings and the absolute corrected error as the scored quantity.
Coupling models:  gauss = linear correction and Gaussian of the PyDP4 source (full coupling);  fc = the same for the
Fermi-contact part with its own parameters;  rec = scaling and density recovered from the deposit outputs (jmodel.py).
Energy-averaged likelihoods: conformer energies get independent Gaussian noise (s kJ/mol), the pairing of signals with
atoms and of couplings with partners stays that of the unperturbed calculation, and exp(L) is averaged over the draws."""
import sys, os, json, math, csv, collections
core2, core3, lib, src_json, truth_csv, model_json, out = sys.argv[1:8]
draws = int(sys.argv[8]) if len(sys.argv) > 8 else 500
seed = int(sys.argv[9]) if len(sys.argv) > 9 else 20261114
sys.path.insert(0, lib)
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import linear_sum_assignment
from scipy.special import log_ndtr, logsumexp
import dep17, desc17, parse17, raw17, core17, jlib
from energy_source17 import renum_maps

RT = 8.3145 * 298.15 / 1000.0
SCALES = [0.0, 2.0, 4.0, 8.0]
LOG2 = math.log(2.0)
S = dep17.Small(core2)
SRC = json.load(open(src_json))
T = {r['experiment']: r for r in csv.DictReader(open(truth_csv, encoding='utf-8'))}
data = raw17.load(core2, {r['archive'] for r in SRC})
M = jlib.load_matrices(core3)
REC = json.load(open(model_json))
_spl = BSpline(np.array(REC['knots']), np.array(REC['coef']), REC['degree'])
X_SUP = 10.0                                              # above this the recovered density rests on a handful of couplings
_k_tail = float(_spl(8.0) - _spl(X_SUP)) / (X_SUP - 8.0)


def log_rec(err):
    """log of the recovered density; beyond X_SUP it falls off at the rate found between 8 Hz and X_SUP"""
    err = np.asarray(err, dtype=float)
    return np.where(err <= X_SUP, _spl(np.minimum(err, X_SUP)), float(_spl(X_SUP)) - _k_tail * (err - X_SUP))


def log_gauss(err, mean, sd):
    return -0.5 * ((err - mean) / sd) ** 2 - math.log(sd * math.sqrt(2 * math.pi))


MODELS = {'gauss': ('total', jlib.PARAMS['total']['slope'], jlib.PARAMS['total']['intercept'], lambda e: log_gauss(e, jlib.PARAMS['total']['mean'], jlib.PARAMS['total']['sd'])),
          'fc': ('FC', jlib.PARAMS['FC']['slope'], jlib.PARAMS['FC']['intercept'], lambda e: log_gauss(e, jlib.PARAMS['FC']['mean'], jlib.PARAMS['FC']['sd'])),
          'rec': ('total', REC['slope'], REC['intercept'], log_rec)}


def assign_idx(labels, calc, exp):
    """core17.assign with the position of each experimental signal in the description kept: [(label, position)]"""
    order = sorted(range(len(exp)), key=lambda i: exp[i][0])
    pairs = core17.assign(labels, calc, exp)
    # core17.assign walks the signals in this order and drops none in these descriptions
    assert len(pairs) == len(order), 'signals without an atom'
    assert all(abs(pairs[v][2] - exp[order[v]][0]) < 1e-9 for v in range(len(order)))
    return [(pairs[v][0], order[v]) for v in range(len(order))]


def shift_loglik(u, exp, sigma):
    mx = exp.mean()
    my = u.mean(axis=1, keepdims=True)
    sxx = ((exp - mx) ** 2).sum()
    slope = ((u - my) * (exp - mx)).sum(axis=1, keepdims=True) / sxx
    icpt = my - slope * mx
    err = (u - icpt) / slope - exp
    return (LOG2 + log_ndtr(-np.abs(err) / sigma)).sum(axis=1)


primary = collections.OrderedDict()
pref = lambda r: 0 if r['file'].endswith('NMRj') else 1 if r['file'].endswith('jreal') else 2
for rec in SRC:
    cur = primary.get(rec['experiment'])
    if cur is None or pref(rec) < pref(cur):
        primary[rec['experiment']] = rec
rng = np.random.default_rng(seed)
os.makedirs(out, exist_ok=True)
fh = open(os.path.join(out, 'jbench_candidates.csv'), 'w', newline='')
w = csv.writer(fh)
head = ['archive', 'experiment', 'compound', 'file', 'candidates', 'candidate', 'structure', 'is_actual', 'renumbered', 'conformers', 'n_C', 'n_H',
        'protons_with_couplings', 'couplings', 'couplings_without_partner', 'logL_C', 'logL_H', 'logL_J_gauss', 'logL_J_fc', 'logL_J_rec']
for s in SCALES[1:]:
    head += ['E%g_CH' % s, 'E%g_CHJ_gauss' % s, 'E%g_CHJ_fc' % s, 'E%g_CHJ_rec' % s]
w.writerow(head)
wc = csv.writer(open(os.path.join(out, 'jbench_couplings.csv'), 'w', newline=''))
wc.writerow(['experiment', 'candidate', 'structure', 'is_actual', 'proton', 'partner', 'experimental_Hz', 'calculated_total_Hz', 'corrected_Hz'])   # matching of the Gaussian model
stat = collections.Counter()
check = []
dup_j = 0
for name, rec in primary.items():
    if name not in T:
        stat['experiments without a known actual structure'] += 1
        continue
    arch, jbase, folder, stems = rec['archive'], rec['jbase'], rec['folder'], rec['stems']
    actual = int(T[name]['actual_isomer'])
    desc = desc17.parse(S.read(rec['path'][:-4]).splitlines())
    printed = parse17.parse(S.read(rec['path']))['isomers']
    block_of = {}
    for it in rec['isomers']:                              # candidate (main order) -> isomer block of the deposit output
        for f in it['fits']:
            if f[2] == 'J' and f[3] == 'J':
                block_of.setdefault(f[1], it['isomer'])
    maps = renum_maps(S, arch, folder)
    seen = {}
    for s_, l_, j_ in desc['H']:
        if s_ in seen and seen[s_] != j_:
            dup_j += 1
        seen.setdefault(s_, j_)
    rows = []
    for si, stem in enumerate(stems, 1):
        lev = data.get((arch, jbase, stem), {}).get('J')
        if not lev:
            stat['candidates without coupling jobs'] += 1
            rows = None
            break
        ks = sorted(lev)
        mats = [M.get((arch, folder, stem, k)) for k in ks]
        if any(m is None for m in mats):
            stat['candidates with a conformer matrix missing'] += 1
            rows = None
            break
        en = np.array([lev[k][0] for k in ks]) * jlib.HARTREE
        en -= en.min()
        Sh = np.array([lev[k][1] for k in ks])
        els = list(lev[ks[0]][2])
        own_h = mats[0]['h']
        mapping = maps[si - 2] if (maps and si >= 2 and len(maps[si - 2]) == Sh.shape[1]) else None
        if mapping is not None:
            Sh = Sh[:, mapping]
            els = [els[i] for i in mapping]
        hpos = {a: i for i, a in enumerate(own_h)}
        new_h = [i for i, e_ in enumerate(els) if e_ == 'H']                       # positions in the numbering of the first candidate
        idx = [hpos[(mapping[i] if mapping is not None else i) + 1] for i in new_h]
        hl = [str(i + 1) for i in new_h]
        at = core17.atoms_of(desc, els)
        lw0 = -en / RT
        w0 = np.exp(lw0 - logsumexp(lw0))
        part, pairing = {}, {}
        for el in 'CH':
            labels = [lab for lab, ix in at[el]]
            G = np.zeros((len(labels), Sh.shape[1]))
            for j, (lab, ix) in enumerate(at[el]):
                for a in ix:
                    G[j, a] = 1.0 / len(ix)
            SL = Sh @ G.T                                                       # conformers x labels
            exp = [(s_, l_) for s_, l_, j_ in desc[el]]
            pr = assign_idx(labels, list(-(w0 @ SL)), exp)
            pairing[el] = pr
            col = [labels.index(lab) for lab, pos in pr]
            part[el] = (SL[:, col], np.array([desc[el][pos][0] for lab, pos in pr]))
        # couplings: pruned matrix of every conformer (pruning is linear in the matrix)
        pruned = {}
        for kind in ('total', 'FC'):
            Pc = []
            for m in mats:
                p_, pl = jlib.zero_equiv(m[kind][np.ix_(idx, idx)], hl, desc['groups'], desc['omit'])
                Pc.append(p_)
            Pc = np.array(Pc)
            pruned[kind] = (Pc, np.tensordot(w0, Pc, axes=1), pl)
        jparts = {}
        for mname, (kind, slope, icpt, logf) in MODELS.items():
            Pc, P0, pl = pruned[kind]
            ev, vals, info, n_prot, n_pad = [], [], [], 0, 0
            for lab, pos in pairing['H']:
                jexp = desc['H'][pos][2]
                if not jexp or lab[1:] not in pl:
                    continue
                r = pl.index(lab[1:])
                cols = [q for q in range(len(pl)) if abs(P0[r, q]) > jlib.J_THRESH]
                cvals = np.array([abs(P0[r, q]) for q in cols] + [0.0] * max(0, len(jexp) - len(cols)))
                cost = -logf(np.abs(np.array(jexp)[:, None] - (cvals[None, :] - icpt) / slope))
                ri, ci = linear_sum_assignment(cost)                         # the matching with the largest likelihood
                n_prot += 1
                for i_, c_ in zip(ri, ci):
                    ev.append(jexp[i_])
                    if c_ < len(cols):
                        vals.append(Pc[:, r, cols[c_]]); info.append((lab, 'H' + pl[cols[c_]], jexp[i_]))
                    else:
                        vals.append(np.zeros(len(ks))); info.append((lab, '', jexp[i_])); n_pad += 1
            jparts[mname] = (np.array(ev), np.array(vals).T if vals else np.zeros((len(ks), 0)), info, n_prot, n_pad)
        res = dict(si=si, stem=stem, n_conf=len(ks), ren=int(mapping is not None), nC=part['C'][0].shape[1], nH=part['H'][0].shape[1],
                   nprot=jparts['gauss'][3], nJ=len(jparts['gauss'][0]), npad=jparts['gauss'][4])
        per_scale = {}
        for s in SCALES:
            m = 1 if s == 0 else draws
            noise = rng.normal(0.0, s, size=(m, len(ks))) if s > 0 else np.zeros((1, len(ks)))
            lw = -(en[None, :] + noise) / RT
            W = np.exp(lw - logsumexp(lw, axis=1, keepdims=True))
            L = {el: shift_loglik(-(W @ part[el][0]), part[el][1], core17.SIGMA[el]) for el in 'CH'}
            LJ = {}
            for mname, (kind, slope, icpt, logf) in MODELS.items():
                ev, V, info, n_prot, n_pad = jparts[mname]
                if len(ev) == 0:
                    LJ[mname] = np.zeros(m)
                    continue
                calc = np.abs(W @ V)                                        # draws x couplings
                LJ[mname] = logf(np.abs(ev[None, :] - (calc - icpt) / slope)).sum(axis=1)
            per_scale[s] = (L, LJ, m)
        L, LJ, _ = per_scale[0.0]
        res.update(logC=float(L['C'][0]), logH=float(L['H'][0]), J={k_: float(v[0]) for k_, v in LJ.items()}, E={})
        for s in SCALES[1:]:
            L, LJ, m = per_scale[s]
            ch = L['C'] + L['H']
            res['E'][s] = [float(logsumexp(ch) - math.log(m))] + [float(logsumexp(ch + LJ[k_]) - math.log(m)) for k_ in ('gauss', 'fc', 'rec')]
        ev, V, info, n_prot, n_pad = jparts['gauss']
        c0 = np.abs(w0 @ V) if len(ev) else []
        for (lab, partner, e_), cval in zip(info, c0):
            wc.writerow([name, si, stem, int(si == actual), lab, partner, e_, '%.3f' % cval, '%.3f' % ((cval - MODELS['gauss'][2]) / MODELS['gauss'][1])])
        # shift terms against the deposit output (same shieldings, energies and model)
        b = block_of.get(si)
        if b is not None and printed[b]['p'].get('C', 0) > 0 and printed[b]['p'].get('H', 0) > 0:
            check.append((res['logC'] - math.log(printed[b]['p']['C']), res['logH'] - math.log(printed[b]['p']['H'])))
        rows.append(res)
    if rows is None:
        stat['experiments left out'] += 1
        continue
    stat['experiments'] += 1
    for r in rows:
        line = [arch, name, T[name]['si_compound'], rec['file'], len(rows), r['si'], r['stem'], int(r['si'] == actual), r['ren'], r['n_conf'], r['nC'], r['nH'],
                r['nprot'], r['nJ'], r['npad'], '%.5f' % r['logC'], '%.5f' % r['logH'], '%.5f' % r['J']['gauss'], '%.5f' % r['J']['fc'], '%.5f' % r['J']['rec']]
        for s in SCALES[1:]:
            line += ['%.5f' % v for v in r['E'][s]]
        w.writerow(line)
        stat['candidates'] += 1
fh.close()
for k, v in sorted(stat.items()):
    print('%5d  %s' % (v, k))
c = np.array(check)
print('shift terms against the deposit outputs (%d candidate blocks): largest difference of log(term) %.4f for carbon, %.4f for hydrogen' % (
    len(c), np.abs(c[:, 0]).max(), np.abs(c[:, 1]).max()))
print('signals with a shift equal to an earlier one but other couplings (where PyDP4 used the couplings of the first): %d' % dup_j)
print('energy noise: %s kJ/mol, %d draws, seed %d' % (', '.join('%g' % s for s in SCALES[1:]), draws, seed))
