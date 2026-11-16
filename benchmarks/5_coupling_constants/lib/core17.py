"""Shared routines for the 2017 deposit: atoms taken into account, pairing with the experimental shifts, DP4."""
import collections
import csv
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dep17
import desc17
import parse17
import raw17
from energy_source17 import LEVELS, renum_maps  # noqa: F401

SIGMA = {'C': 2.269372270818724, 'H': 0.18731058105269952}
LEVEL_NAMES = dep17.LEVEL_NAMES


def num(label):
    return int(re.sub(r'\D', '', label)) - 1


def atoms_of(desc, elements):
    """-> {'C': [(label, [atom indices])], 'H': [...]} in atom order"""
    omit = set(desc['omit'])
    group = {}
    for g in desc['groups']:
        for lab in g:
            group[lab] = g
    out = {'C': [], 'H': []}
    for i, el in enumerate(elements):
        if el not in out:
            continue
        lab = '%s%d' % (el, i + 1)
        if lab in omit:
            continue
        g = group.get(lab)
        if g is None:
            out[el].append((lab, [i]))
        elif g[0] == lab:
            out[el].append((lab, [num(x) for x in g]))
    return out


def assign(labels, calc, exp):
    """PyDP4 (2015-2017) AssignExpNMR. exp = [(shift, label text)] -> [(label, calc, exp)] in ascending exp order"""
    order = sorted(range(len(exp)), key=lambda i: exp[i][0])
    exp_sorted = [exp[i][0] for i in order]
    exp_labels = []
    for i in order:
        t = [x.strip() for x in re.split(r',|\bor\b', exp[i][1]) if x.strip()]
        exp_labels.append([] if (not t or 'any' in t) else t)
    calc_order = sorted(range(len(calc)), key=lambda i: calc[i])
    calc_labels = [labels[i] for i in calc_order]
    given = [''] * len(exp_sorted)
    for v, t in enumerate(exp_labels):
        if len(t) == 1:
            if t[0] not in labels:
                raise KeyError(t[0])
            given[v] = t[0]
    for v, t in enumerate(exp_labels):
        if len(t) > 1:
            for lab in calc_labels:
                if lab in t and lab not in given:
                    given[v] = lab
                    break
    for v in range(len(exp_sorted)):
        if given[v] == '':
            for lab in calc_labels:
                if lab not in given:
                    given[v] = lab
                    break
    pos = {lab: i for i, lab in enumerate(labels)}
    return [(lab, calc[pos[lab]], e) for lab, e in zip(given, exp_sorted) if lab]


def scale(calc, exp):
    n = len(calc)
    mx, my = sum(exp) / n, sum(calc) / n
    sxx = sum((a - mx) ** 2 for a in exp)
    sxy = sum((a - mx) * (b - my) for a, b in zip(exp, calc))
    slope = sxy / sxx
    icpt = my - slope * mx
    return [(c - icpt) / slope for c in calc]


def log_prob(errors, sigma):
    return sum(math.log(max(math.erfc(abs(e) / sigma / math.sqrt(2)), 1e-300)) for e in errors)


def normalise(logs):
    m = max(logs)
    w = [math.exp(x - m) for x in logs]
    z = sum(w)
    return [x / z for x in w]


def reference_constants(source, solvents=None):
    """Reference constants (shift = (constant - shielding) / (1 - constant * 1e-6)) as found in the deposit outputs:
    the median over the isomer blocks whose shifts are reproduced, per calculation folder (molecule archive and base
    folder; one solvent each), level and element. With solvents = {(archive, base): solvent} a missing value is taken
    from the folders calculated in the same solvent."""
    acc = collections.defaultdict(list)
    for r in source:
        if r['reproduced_with'] and r['const']:
            for el, v in r['const'].items():
                acc[(r['archive'], r['base'], r['level'], el)].append(v)
    const = {k: sorted(v)[len(v) // 2] for k, v in acc.items()}
    if solvents:
        pool = collections.defaultdict(list)
        for (a, b, l, el), v in const.items():
            pool[(solvents.get((a, b)), l, el)].append(v)
        for (a, b), solv in solvents.items():
            for l in LEVELS:
                for el in 'CH':
                    if (a, b, l, el) not in const and pool.get((solv, l, el)):
                        v = pool[(solv, l, el)]
                        const[(a, b, l, el)] = sorted(v)[len(v) // 2]
    return const


def solvents_of(structures):
    """{(archive, base folder): solvent named by Gaussian} from the table of structures"""
    import gzip
    out = {}
    with gzip.open(structures, 'rt') as f:
        for line in f:
            r = json.loads(line)
            if r.get('solvent') and not r['sub']:
                out.setdefault((r['archive'], r['base']), r['solvent'])
    return out


def candidates(S, data, e, shield_level, energy_level, const):
    """All candidates of experiment e with the shieldings of one level and the Boltzmann weights of another.
    -> (description, [dict(isomer, stem, n_conf, C=[(label, shielding, unscaled, scaled, exp)], H=[...], log={'C': .., 'H': ..})])"""
    path = e['levels'][shield_level]
    desc = desc17.parse(S.read(path[:-4]).splitlines())
    maps = renum_maps(S, e['archive'], dep17.split(path)[1])
    exp = {'C': [(s, l) for s, l, j in desc['C']], 'H': [(s, l) for s, l, j in desc['H']]}
    out = []
    for k, stem in enumerate(e['stems'], 1):
        a = raw17.average(data.get((e['archive'], e['base'], stem), {}), shield_level, energy_level)
        if a is None:
            return desc, None
        sig, els, ks = a
        if maps and k >= 2 and len(maps[k - 2]) == len(sig):
            sig, els = raw17.renumber(sig, maps[k - 2]), raw17.renumber(els, maps[k - 2])
        at = atoms_of(desc, els)
        row = dict(isomer=k, stem=stem, n_conf=len(ks), log={})
        for el in 'CH':
            labels = [lab for lab, idx in at[el]]
            shield = [sum(sig[i] for i in idx) / len(idx) for lab, idx in at[el]]
            c = const.get((e['archive'], e['base'], shield_level, el))
            unscaled = [(c - x) / (1 - c * 1e-6) for x in shield] if c else [-x for x in shield]
            pairs = assign(labels, unscaled, exp[el])
            sc = scale([p[1] for p in pairs], [p[2] for p in pairs])
            sh = dict(zip(labels, shield))
            row[el] = [(lab, sh[lab], u if c else None, s2, x) for (lab, u, x), s2 in zip(pairs, sc)]
            row['log'][el] = log_prob([s2 - x for (lab, u, x), s2 in zip(pairs, sc)], SIGMA[el])
        out.append(row)
    return desc, out
