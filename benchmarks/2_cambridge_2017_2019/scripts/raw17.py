"""Conformer-level data of the 2017 deposit and Boltzmann averaging with a chosen energy source."""
import collections
import json
import math
import os
import re

import dep17
from runs17 import level_of

R, T, HARTREE = 8.3145, 298.15, 2625.499629554010
NAME = re.compile(r'^(?P<stem>.+?)ginp(?P<conf>\d+)\.out$')


def load(core, archives=None):
    """-> {(archive, base, stem): {level: {conformer: (scf energy, [shieldings], [elements], ended normally)}}}
    only files that sit directly in a level folder"""
    data = collections.defaultdict(lambda: collections.defaultdict(dict))
    for fn in dep17.gaussian_files(core):
        for r in dep17.jsonl(fn):
            if 'iso' not in r or 'scf' not in r:
                continue
            s = dep17.split(r['path'])
            if s is None or '/' in s[1]:
                continue
            arch, folder, name = s
            if archives and arch not in archives:
                continue
            m = NAME.match(name)
            if not m:
                continue
            base, level = level_of(folder)
            iso = json.loads(r['iso']) if isinstance(r['iso'], str) else r['iso']
            d = data[(arch, base, m.group('stem'))][level]
            d.setdefault(int(m.group('conf')), (float(r['scf']), iso, r['el'].split(), str(r.get('norm')) == '1'))
    return data


def average(levels, shield_level, energy_level):
    """Boltzmann average of the shieldings of one level with the energies of another (or the same) level.
    -> (shieldings, elements, conformers used) or None"""
    a = levels.get(shield_level)
    b = levels.get(energy_level)
    if not a or not b:
        return None
    ks = sorted(k for k in a if k in b)
    if not ks or len(ks) != len(a):
        return None
    e = [b[k][0] for k in ks]
    e0 = min(e)
    w = [math.exp(-(x - e0) * HARTREE * 1000 / (R * T)) for x in e]
    z = sum(w)
    n = len(a[ks[0]][1])
    if any(len(a[k][1]) != n for k in ks):
        return None
    return [sum(wi / z * a[k][1][i] for wi, k in zip(w, ks)) for i in range(n)], a[ks[0]][2], ks


def renumber(values, mapping):
    """PyDP4 atom renumbering: position i takes the value of atom mapping[i] (0-based)."""
    return [values[mapping[i]] for i in range(len(values))]
