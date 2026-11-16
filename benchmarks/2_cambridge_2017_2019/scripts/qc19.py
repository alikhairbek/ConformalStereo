"""Quality checks on the raw files of the 2019 deposit, per candidate structure and kind of geometry.
  - conformers with shieldings at each of the three NMR levels, and whether the three levels share one geometry
    (principal moments of inertia printed by Gaussian)
  - optimisation jobs (files ...temp.out): found, ended normally with a converged geometry
  - conformers that fell into the same minimum (same energy to 1e-6 hartree at the first NMR level)
Usage (with CS_DEPOSIT=conf2018): qc19.py <res/core> <experiments.json> <out.csv>
"""
import collections
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dep17

NMR = re.compile(r'^(?P<stem>.+?)ginp(?P<conf>\d+)\.out$')
OPT = re.compile(r'^(?P<stem>.+?)ginp(?P<conf>\d+)temp\.out$')
GEO = [('', 'MMFF'), ('s', 'PM6'), ('pm7', 'PM7'), ('h', 'HF/6-31G(d,p)'), ('d', 'B3LYP/6-31G(d,p)'), ('m', 'M06-2X/6-31G(d,p)')]


def main(core, experiments, out):
    E = json.load(open(experiments))
    want = {(e['archive'], e['base'], s) for e in E for s in e['stems']}
    nmr = collections.defaultdict(lambda: collections.defaultdict(dict))
    opt = collections.defaultdict(lambda: collections.defaultdict(dict))
    for fn in dep17.gaussian_files(core):
        for r in dep17.jsonl(fn):
            s = dep17.split(r['path'])
            if s is None or '/' in s[1]:
                continue
            base, level = dep17.level_of(s[1])
            if level not in dep17.LEVELS:
                continue
            m = NMR.match(s[2])
            if m and (s[0], base, m.group('stem')) in want and 'iso' in r and 'scf' in r:
                pm = json.loads(r['pm']) if isinstance(r['pm'], str) else r['pm']
                nmr[(s[0], base, m.group('stem'))][level].setdefault(int(m.group('conf')), (float(r['scf']), pm, str(r.get('norm')) == '1'))
                continue
            m = OPT.match(s[2])
            if m and (s[0], base, m.group('stem')) in want:
                ok = str(r.get('opt_done')) == '1' and str(r.get('norm')) == '1'
                d = opt[(s[0], base, m.group('stem'))][level]
                d[int(m.group('conf'))] = d.get(int(m.group('conf')), False) or ok
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['archive', 'folder', 'structure', 'geometry', 'conformers', 'nmr_levels_complete', 'same_geometry_at_three_levels',
                'optimisation_jobs', 'optimisations_not_converged', 'distinct_minima'])
    tot = collections.defaultdict(collections.Counter)
    for key in sorted(want):
        for g, gname in GEO:
            levels = ['BF1' + g, 'BF4' + g, 'BF7' + g]
            sets = [sorted(nmr[key].get(l, {})) for l in levels]
            if not any(sets):
                continue
            ref = nmr[key].get(levels[0]) or nmr[key].get(levels[2]) or {}
            confs = sorted(ref)
            complete = sum(1 for s in sets if s == confs and s)
            same = True
            for l in levels[1:]:
                for k in confs:
                    v = nmr[key].get(l, {}).get(k)
                    if v and any(b and abs(a - b) / b > 1e-4 for a, b in zip(v[1], ref[k][1])):
                        same = False
            o = {}
            for l in levels:
                for k, ok in opt[key].get(l, {}).items():
                    o[k] = o.get(k, False) or ok
            bad = sum(1 for k in confs if k in o and not o[k])
            minima = len({round(ref[k][0], 6) for k in confs})
            w.writerow([key[0], key[1], key[2], gname, len(confs), complete, int(same), len(o), bad, minima])
            t = tot[gname]
            t['structures'] += 1
            t['conformers'] += len(confs)
            t['complete'] += complete == 3
            t['same'] += same
            t['opt'] += len([k for k in confs if k in o])
            t['bad'] += bad
            t['bad_structures'] += bad > 0
            t['minima'] += minima
    print('%-20s %10s %10s %14s %14s %12s %16s %22s %16s' % ('geometry', 'structures', 'conformers', 'three levels', 'same geometry', 'opt jobs',
                                                             'not converged', 'structures affected', 'distinct minima'))
    for g, gname in GEO:
        t = tot[gname]
        print('%-20s %10d %10d %14d %14d %12d %16d %22d %16d' % (gname, t['structures'], t['conformers'], t['complete'], t['same'], t['opt'],
                                                               t['bad'], t['bad_structures'], t['minima']))


if __name__ == '__main__':
    main(*sys.argv[1:4])
