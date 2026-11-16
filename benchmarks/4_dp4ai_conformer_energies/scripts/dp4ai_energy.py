"""DP4 for the DP4-AI deposit with a pairing rule for incomplete shift lists and with the energy-averaged likelihood.

For every molecule (46) and protocol (4: MMFF geometries; B3LYP geometries with B3LYP, mPW1PW91 or M06-2X energies):
  - the log-likelihood of every candidate as in the 2020 code (single Gaussian, internal scaling, sorted pairing)
  - where the description lists fewer or more values than there are atoms: the order-preserving pairing with the
    highest likelihood (pairing.py) instead of the 2020 rule of pairing from the top
  - the likelihood averaged over Gaussian noise on the conformer energies (0, 2, 4, 8, 16 kJ/mol; 200 draws); in each
    draw the shifts are sorted and paired again, within the atoms and signals kept by the pairing
Usage: dp4ai_energy.py <res/core of Kaggle run 1> <tables dir of the DP4-AI archive> <out dir> [draws] [seed]
(needs numpy and scipy)
"""
import collections
import csv
import json
import math
import os
import re
import sys

import numpy as np
from scipy.special import log_ndtr, logsumexp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dp4ai_data as A
import dp4core as D
import dp4lib as L
import pairing

PROTOCOLS = [('MM', 'MM', 'nmr', 'MMFF'), ('DFT_EB3LYP', 'Opt', 'opt', 'B3LYP/6-31G(d,p)'),
             ('DFT_EmPW', 'Opt', 'nmr', 'B3LYP/6-31G(d,p)'), ('DFT_EM062X', 'Opt_E', 'e', 'B3LYP/6-31G(d,p)')]
ENERGIES = {'MM': 'NMR level', 'DFT_EB3LYP': 'B3LYP', 'DFT_EmPW': 'NMR level', 'DFT_EM062X': 'M06-2X/def2-TZVP'}
TMS = (188.452125, 32.1243166667)
SIGMA = {'C': D.SIGMA_C, 'H': D.SIGMA_H}
SCALES = [0.0, 2.0, 4.0, 8.0, 16.0]
RT = A.R * A.T / 1000.0
LOG2 = math.log(2.0)


def loglik(u, exp, sigma):
    """u: draws x pairs, each row ascending; exp ascending -> log-likelihood per draw (internal scaling)"""
    mx = exp.mean()
    my = u.mean(axis=1, keepdims=True)
    slope = ((u - my) * (exp - mx)).sum(axis=1, keepdims=True) / ((exp - mx) ** 2).sum()
    icpt = my - slope * mx
    err = (u - icpt) / slope - exp
    return (LOG2 + log_ndtr(-np.abs(err) / sigma)).sum(axis=1)


def main(core, work, out, draws=200, seed=20261114):
    draws, seed = int(draws), int(seed)
    dep = A.Deposit(core)
    S = L.Small(core, 'dp4ai')
    ps = S.paths()
    truth = json.load(open(os.path.join(work, 'dp4ai_truth.json')))
    # geometry check of the DP4-AI archive: is every conformer that carries population on the intended geometry?
    CLEAN = {(r['case'], r['protocol']): int(r['qc_clean']) for r in csv.DictReader(open(os.path.join(work, 'dp4ai_cases.csv')))}
    rng = np.random.default_rng(seed)

    def description(case):
        seen = collections.OrderedDict()
        for folder in ('MM', 'Opt', 'Opt_E'):
            for p in ps:
                if re.match(r'^%s\.zip!/[^/]+/%s/%sNMR$' % (re.escape(case), folder, re.escape(case)), p):
                    try:
                        d = L.parse_desc(S.read(p))
                    except Exception:
                        continue
                    if d['bad'] or len(d['H']) < 2:
                        continue
                    seen.setdefault((tuple(d['C']), tuple(d['H'])), []).append((folder, d))
        if not seen:
            return None
        key, lst = max(seen.items(), key=lambda kv: (len(kv[1]), kv[1][0][0] == 'MM'))
        return lst[0][1]

    os.makedirs(out, exist_ok=True)
    cols = ['archive', 'experiment', 'shielding_level', 'energy_level', 'geometry', 'nmr', 'candidates', 'isomer', 'structure',
            'n_C', 'n_H', 'conformers', 'exact_counts', 'qc_clean', 'logL_2020_rule'] + ['logL_noise_%g' % s for s in SCALES]
    rows = []
    tr = {}
    worst = 0.0
    changed = []
    for case in dep.cases():
        desc = description(case)
        if desc is None:
            continue
        fam = truth[case]['_summary']['family']
        TRUE = truth[case]['_summary']['true_isomer']
        for prot, folder, how, geo in PROTOCOLS:
            isos = dep.isomers(case, folder)
            if TRUE not in isos:
                continue
            cand = []
            ok = True
            for k in isos:
                av = dep.averaged(case, folder, k, how)
                if av is None:
                    ok = False
                    break
                els, sig, pop = av
                src = dep.E[(case, folder, how)]
                use = sorted(pop)
                M = np.array([dep.SIG[(case, folder)][(k, c)][1] for c in use])
                en = np.array([src[(k, c)] for c in use]) * A.HARTREE
                cand.append(dict(iso=k, els=els, sig=sig, M=M, en=en - en.min()))
            if not ok:
                continue
            built = []
            for c in cand:
                one = {}
                for nuc in 'CH':
                    idx = [a for a, e in enumerate(c['els']) if e == nuc]
                    one[nuc] = (['%s%d' % (nuc, a + 1) for a in idx], [L.to_shift(c['sig'][a], nuc, TMS) for a in idx])
                built.append(one)
            old = D.dp4(built, desc['C'], desc['H'], rule='2020')
            exact = True
            use_n = {}
            for c, b in zip(cand, built):
                c['keep'] = {}
                for nuc in 'CH':
                    exp = desc[nuc]
                    idx = [a for a, e in enumerate(c['els']) if e == nuc]
                    if not exp or len(exp) < 2 or not idx:
                        continue
                    if len(exp) == len(idx):
                        ka, ke = list(range(len(idx))), list(range(len(exp)))
                    else:
                        exact = False
                        ka, ke = pairing.optimal(b[nuc][1], exp, SIGMA[nuc])
                    c['keep'][nuc] = (np.array([idx[i] for i in ka]), np.sort(np.array([exp[j] for j in ke])))
                    use_n[nuc] = len(ka)
            L0 = []
            for c in cand:
                vals = []
                for s in SCALES:
                    m = 1 if s == 0 else draws
                    noise = rng.normal(0.0, s, size=(m, len(c['en']))) if s > 0 else np.zeros((1, len(c['en'])))
                    lw = -(c['en'][None, :] + noise) / RT
                    lw -= logsumexp(lw, axis=1, keepdims=True)
                    W = np.exp(lw)
                    tot = np.zeros(m)
                    for nuc, (atoms, exp) in c['keep'].items():
                        u = np.sort(-(W @ c['M'][:, atoms]), axis=1)
                        tot += loglik(u, exp, SIGMA[nuc])
                    vals.append(float(logsumexp(tot) - math.log(m)))
                c['vals'] = vals
                L0.append(vals[0])
            old_log = [sum(old['logs'][n][i] for n in old['logs']) for i in range(len(cand))] if 'logs' in old else [float('nan')] * len(cand)
            if exact:
                worst = max(worst, max(abs(a - b) for a, b in zip(L0, old_log)))
            else:
                t = isos.index(TRUE)
                p_new = D.normalise(L0)
                p_old = old['all']
                changed.append((case, prot, 1 + sum(1 for x in p_old if x > p_old[t]), 100 * p_old[t], 1 + sum(1 for x in p_new if x > p_new[t]), 100 * p_new[t]))
            clean = CLEAN.get((case, prot), 0)
            for pos, c in enumerate(cand):
                rows.append([fam, case, prot, ENERGIES[prot], geo, 'mPW1PW91', len(cand), pos + 1, '%s_%d' % (case, c['iso']),
                             use_n.get('C', 0), use_n.get('H', 0), len(c['en']), int(exact), clean, '%.5f' % old_log[pos]]
                            + ['%.5f' % v for v in c['vals']])
            tr[case] = isos.index(TRUE) + 1
    for name, keep in (('dp4ai_energy_marginal.csv', lambda r: True), ('dp4ai_energy_marginal_clean.csv', lambda r: r[12] == 1 and r[13] == 1)):
        w = csv.writer(open(os.path.join(out, name), 'w', newline=''))
        w.writerow(cols)
        w.writerows(r for r in rows if keep(r))
    w = csv.writer(open(os.path.join(out, 'dp4ai_truth.csv'), 'w', newline=''))
    w.writerow(['experiment', 'actual_isomer'])
    for case in sorted(tr):
        w.writerow([case, tr[case]])
    recs = {(r[1], r[2]) for r in rows}
    clean = {(r[1], r[2]) for r in rows if r[12] == 1 and r[13] == 1}
    print('molecules %d, records (molecule x protocol) %d, of them clean with one value per atom %d' % (len(tr), len(recs), len(clean)))
    print('descriptions with one value per atom: without noise the log-likelihood differs from the 2020 rule by at most %.1e' % worst)
    print()
    print('descriptions with fewer or more values than atoms: the reported structure under the 2020 rule and under the optimal pairing')
    print('%-6s %-11s  2020 rule: rank, probability   optimal pairing: rank, probability' % ('case', 'protocol'))
    for case, prot, r0, p0, r1, p1 in changed:
        print('%-6s %-11s  %10d %12.1f%%  %21d %12.1f%%' % (case, prot, r0, p0, r1, p1))
    print('reported structure first: %d of %d records under the 2020 rule, %d under the optimal pairing'
          % (sum(1 for c in changed if c[2] == 1), len(changed), sum(1 for c in changed if c[4] == 1)))


if __name__ == '__main__':
    main(*sys.argv[1:6])
