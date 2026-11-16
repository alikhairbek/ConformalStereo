"""Coupling constants of the 2017 deposit: conformer matrices, and the coupling term of PyDP4 (2015-2017) step by step.

The routines zero_equiv, assign_shifts_j, assign_j and coupling_errors follow the PyDP4 source of the time
(repository KristapsE/PyDP4, files NMRAnalysis.py and DP4.py as of 18 April 2017; copyright 2015 K. Ermanis and
J. M. Goodman), including the points where that code is not what one would write today; they are marked."""
import gzip
import glob
import itertools
import json
import math
import os
import re

import numpy as np

R, T, HARTREE = 8.3145, 298.15, 2625.499629554010
J_THRESH = 0.2
# linear correction (calculated = slope * experimental + intercept) and Gaussian of the corrected absolute errors
PARAMS = {'total': dict(slope=1.1598043367, intercept=-0.365587317, mean=-0.13133138905769429, sd=0.79315485469419067),
          'FC': dict(slope=1.1701625846, intercept=-0.438405757, mean=-0.15436128540661589, sd=0.92117647579294348)}
NAME = re.compile(r'^(?P<stem>.+?)ginp(?P<conf>\d+)\.out$')


def load_matrices(core3, key='dp4x2'):
    """-> {(archive, folder, stem, conformer): dict(h=[atom numbers of the hydrogens], total=array, FC=array)}
    total: step 2 of the job; FC: step 1 (the first Fermi-contact matrix of the file, the one PyDP4 read)."""
    out = {}
    tmp = {}
    for fn in sorted(glob.glob(os.path.join(core3, 'jcoup_%s_*.jsonl.gz' % key))):
        with gzip.open(fn, 'rt', encoding='utf-8') as fh:
            for line in fh:
                r = json.loads(line)
                p = r['path'].split('!/')
                if len(p) != 3:
                    continue
                arch = p[1].split('/')[-1].replace('.zip', '')
                folder, _, name = p[2].rpartition('/')
                m = NAME.match(name)
                if not m or '/' in folder:
                    continue
                k = (arch, folder, m.group('stem'), int(m.group('conf')))
                el = r['el'].split()
                hs = [i for i, e in enumerate(el) if e == 'H']
                rec = tmp.setdefault(k, dict(h=[i + 1 for i in hs]))
                for name_, want_job in (('FC', 0), ('total', 1)):
                    if r['job'] == want_job and name_ in r['J']:
                        b = r['J'][name_]
                        assert b['lab'] == list(range(1, len(el) + 1))
                        v = b['v']
                        a = np.zeros((len(hs), len(hs)))
                        for x, i in enumerate(hs):
                            base = i * (i + 1) // 2
                            for y in range(x):
                                a[x, y] = a[y, x] = v[base + hs[y]]
                        rec[name_] = a
    for k, rec in tmp.items():
        if 'FC' in rec and 'total' in rec:
            out[k] = rec
    return out


def populations(energies):
    e0 = min(energies)
    w = [math.exp(-(x - e0) * HARTREE * 1000 / (R * T)) for x in energies]
    z = sum(w)
    return [x / z for x in w]


def zero_equiv(mat, labels, equivs, omits):
    """PyDP4 ZeroEquivJ for one candidate. mat: square array over the hydrogens, labels: their atom numbers as strings,
    equivs: groups of labels such as ['H13', 'H14', 'H15'], omits: labels. -> (pruned array, pruned labels)"""
    m = np.array(mat, dtype=float)
    remove = []
    for g in equivs:
        if g[0][0] != 'H' or not all(x[1:] in labels for x in g):
            continue
        idx = [labels.index(x[1:]) for x in g]
        for r in idx:
            for c in idx:
                m[r, c] = 0.0
        first = idx[0]
        for r in range(len(m)):                               # columns first, then rows, as in the source
            m[r, first] = np.mean([m[r, x] for x in idx])
        for c in range(len(m)):
            m[first, c] = np.mean([m[x, c] for x in idx])
        remove.extend(idx[1:])
    remove.extend(labels.index(x[1:]) for x in omits if x[0] == 'H' and x[1:] in labels)
    keep = [i for i in range(len(m)) if i not in remove]
    return m[np.ix_(keep, keep)], [labels[i] for i in keep]


def assign_j(exp_j, calc_rows):
    """PyDP4 AssignJvals. exp_j: list per proton of experimental couplings ([0.0] = none); calc_rows: list per proton
    of the calculated couplings to the other protons. -> (experimental lists as used, assigned calculated lists)"""
    pruned = [[abs(j) for j in row if abs(j) > J_THRESH] for row in calc_rows]
    exp_used = [list(e) for e in exp_j]
    for i in range(len(exp_j)):
        if len(pruned[i]) == 0:
            exp_used[i] = [0.0]
    assigned = []
    for e, c in zip(exp_used, pruned):
        if e != [0.0]:
            best, best_err = [0.0] * len(e), 10000
            for perm in itertools.permutations(c, len(e)):     # none if fewer calculated than experimental couplings
                err = sum(abs(perm[i] - e[i]) for i in range(len(e)))
                if err < best_err:
                    best_err, best = err, list(perm)
            assigned.append(best)
        else:
            assigned.append([0.0])
    return exp_used, assigned


def coupling_errors(exp_used, assigned, slope, intercept):
    """Absolute errors after the linear correction, as PyDP4 CalculateJerrors after ScaleJvals."""
    out = []
    for e, c in zip(exp_used, assigned):
        if e != [0.0]:
            out.extend(abs(e[i] - (c[i] - intercept) / slope) for i in range(len(e)))
    return out


def log_gauss(errors, mean, sd):
    """Sum of log densities of a Gaussian (PyDP4 CalculatePDP4 multiplies the densities)."""
    return sum(-0.5 * ((e - mean) / sd) ** 2 - math.log(sd * math.sqrt(2 * math.pi)) for e in errors)
