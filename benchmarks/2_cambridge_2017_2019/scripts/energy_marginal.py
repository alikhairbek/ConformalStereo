"""A score that allows for errors in the conformer energies.

DP4 weights the conformers of a candidate with Boltzmann factors of computed energies and treats the weights as exact.
Here the energies of the conformers get independent Gaussian noise (standard deviation s kJ/mol), the shifts and the
DP4 log-likelihood are recomputed for each draw, and the candidate's likelihood is averaged over the draws:
    L_s(candidate) = log mean_draws exp(L(candidate | energies + noise))
A candidate whose good fit hangs on one fragile set of populations loses; one that fits under any plausible populations
keeps its likelihood. The pairing of shifts and atoms is the one of the unperturbed calculation.
Usage: energy_marginal.py <res/core> <experiments.json> <energy_source.json> <structures.jsonl.gz> <truth.csv> <out.csv> [draws] [seed]
(needs numpy and scipy; CS_DEPOSIT selects the deposit)
"""
import csv
import json
import math
import os
import re
import sys

import numpy as np
from scipy.special import log_ndtr, logsumexp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core17
import dep17
import raw17
from energy_source17 import renum_maps

RT = 8.3145 * 298.15 / 1000.0          # kJ/mol
HARTREE = 2625.499629554010
SCALES = [0.0, 2.0, 4.0, 8.0, 16.0]
LOG2 = math.log(2.0)


def loglik(u, exp, sigma):
    """u: draws x pairs (unscaled values, any affine convention); exp: pairs -> log-likelihood per draw"""
    mx = exp.mean()
    my = u.mean(axis=1, keepdims=True)
    sxx = ((exp - mx) ** 2).sum()
    slope = ((u - my) * (exp - mx)).sum(axis=1, keepdims=True) / sxx
    icpt = my - slope * mx
    err = (u - icpt) / slope - exp
    return (LOG2 + log_ndtr(-np.abs(err) / sigma)).sum(axis=1)


def main(core, experiments, source, structures, truth, out, draws=200, seed=20261114):
    draws, seed = int(draws), int(seed)
    S = dep17.Small(core)
    T = {r['experiment']: int(r['actual_isomer']) for r in csv.DictReader(open(truth))}
    E = [e for e in json.load(open(experiments)) if e['name'] in T]
    const = core17.reference_constants(json.load(open(source)), core17.solvents_of(structures))
    data = raw17.load(core, {e['archive'] for e in E})
    rng = np.random.default_rng(seed)
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['archive', 'experiment', 'shielding_level', 'energy_level', 'geometry', 'nmr', 'candidates', 'isomer', 'structure',
                'n_C', 'n_H', 'conformers'] + ['logL_noise_%g' % s for s in SCALES])
    worst = 0.0
    n = 0
    for e in E:
        for level in dep17.LEVELS:
            if level not in e['levels']:
                continue
            desc, rows = core17.candidates(S, data, e, level, level, const)
            if rows is None:
                continue
            maps = renum_maps(S, e['archive'], dep17.split(e['levels'][level])[1])
            for row in rows:
                k = row['isomer']
                conf = data[(e['archive'], e['base'], row['stem'])][level]
                keys = sorted(conf)
                M = np.array([conf[q][1] for q in keys])                 # conformers x atoms
                els = conf[keys[0]][2]
                if maps and k >= 2 and len(maps[k - 2]) == M.shape[1]:
                    mp = maps[k - 2]
                    M = M[:, mp]
                    els = [els[i] for i in mp]
                en = np.array([conf[q][0] for q in keys]) * HARTREE
                en = en - en.min()
                at = core17.atoms_of(desc, els)
                parts = {}
                for el in 'CH':
                    idx = dict(at[el])
                    G = np.zeros((len(row[el]), M.shape[1]))
                    for j, pair in enumerate(row[el]):
                        for a in idx[pair[0]]:
                            G[j, a] = 1.0 / len(idx[pair[0]])
                    parts[el] = (G, np.array([pair[4] for pair in row[el]]))
                vals = []
                for s in SCALES:
                    m = 1 if s == 0 else draws
                    noise = rng.normal(0.0, s, size=(m, len(keys))) if s > 0 else np.zeros((1, len(keys)))
                    lw = -(en[None, :] + noise) / RT
                    lw -= logsumexp(lw, axis=1, keepdims=True)
                    sig = np.exp(lw) @ M                                    # draws x atoms
                    L = np.zeros(m)
                    for el in 'CH':
                        G, exp = parts[el]
                        L += loglik(-(sig @ G.T), exp, core17.SIGMA[el])
                    vals.append(float(logsumexp(L) - math.log(m)))
                worst = max(worst, abs(vals[0] - (row['log']['C'] + row['log']['H'])))
                w.writerow([e['archive'], e['name'], level, level, dep17.GEOMETRY[level], dep17.NMR[level], len(rows), k, row['stem'],
                            len(row['C']), len(row['H']), len(keys)] + ['%.5f' % v for v in vals])
            n += 1
    print('records %d; without noise the log-likelihood differs from the main calculation by at most %.1e' % (n, worst))


if __name__ == '__main__':
    main(*sys.argv[1:9])
