"""Boltzmann-averaged shieldings of every structure of the 2017 deposit, from the parsed Gaussian outputs.

One record per (molecule archive, folder, structure): the conformers found (files <structure>ginp<nnn>.out in the
folder itself), their SCF energies at the NMR level, populations at 298.15 K and the averaged isotropic shieldings.
Usage: build17.py <res/core> <out.jsonl.gz>
"""
import collections
import gzip
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dep17
from runs17 import level_of

R, T, HARTREE = 8.3145, 298.15, 2625.499629554010
NAME = re.compile(r'^(?P<stem>.+?)ginp(?P<conf>\d+)\.out$')


def main(core, out):
    groups = collections.defaultdict(dict)
    skipped = collections.Counter()
    for fn in dep17.gaussian_files(core):
        for r in dep17.jsonl(fn):
            s = dep17.split(r['path'])
            if s is None:
                skipped['not in a molecule archive'] += 1
                continue
            arch, folder, name = s
            m = NAME.match(name)
            if not m:
                skipped['other file name'] += 1
                continue
            key = (arch, folder, m.group('stem'))
            conf = int(m.group('conf'))
            c = groups[key].setdefault(conf, {'steps': 0})
            c['steps'] += 1
            if 'iso' in r and 'scf' in r and 'iso' not in c:
                c['iso'] = json.loads(r['iso']) if isinstance(r['iso'], str) else r['iso']
                c['el'] = r['el'].split()
                c['scf'] = float(r['scf'])
                c['route'] = r['route']
                c['solv'] = r.get('solv', '')
                c['norm'] = str(r.get('norm')) == '1'
                c['nbf'] = int(r['nbf'])
    n = 0
    stats = collections.Counter()
    with gzip.open(out, 'wt') as f:
        for (arch, folder, stem) in sorted(groups):
            confs = groups[(arch, folder, stem)]
            use = sorted(k for k, c in confs.items() if 'iso' in c)
            base, level = level_of(folder)
            rec = dict(archive=arch, folder=folder, base=base, level=level, sub='/'.join(folder.split('/')[1:]),
                       stem=stem, n_files=len(confs), n_conf=len(use))
            if not use:
                stats['structures without shieldings'] += 1
                rec['sigma'] = None
                f.write(json.dumps(rec) + '\n')
                continue
            first = confs[use[0]]
            same = [k for k in use if confs[k]['el'] == first['el'] and len(confs[k]['iso']) == len(first['iso'])]
            if len(same) != len(use):
                stats['structures with conformers of another atom list'] += 1
            use = same
            e = [confs[k]['scf'] for k in use]
            e0 = min(e)
            rel = [(x - e0) * HARTREE for x in e]
            w = [math.exp(-x * 1000 / (R * T)) for x in rel]
            z = sum(w)
            pop = [x / z for x in w]
            sigma = [sum(p * confs[k]['iso'][a] for p, k in zip(pop, use)) for a in range(len(first['iso']))]
            # the same shieldings weighted with the energies of the first level (B3LYP/6-31G(d,p)) of the molecule
            other = groups.get((arch, '/'.join([base] + folder.split('/')[1:]), stem), {})
            sigma1 = None
            if level != 'BF1' and all(k in other and 'scf' in other[k] for k in use):
                e1 = [other[k]['scf'] for k in use]
                m1 = min(e1)
                w1 = [math.exp(-(x - m1) * HARTREE * 1000 / (R * T)) for x in e1]
                z1 = sum(w1)
                sigma1 = [round(sum(p / z1 * confs[k]['iso'][a] for p, k in zip(w1, use)), 5) for a in range(len(first['iso']))]
            rec.update(elements=first['el'], sigma=[round(x, 5) for x in sigma], sigma_e1=sigma1, conformers=use,
                       energies=e, populations=[round(p, 6) for p in pop], route=first['route'], solvent=first['solv'],
                       nbf=first['nbf'], routes=len({confs[k]['route'].lower().replace(' ', '') for k in use}),
                       not_normal=sum(1 for k in use if not confs[k]['norm']))
            if rec['routes'] > 1:
                stats['structures with more than one route'] += 1
            if rec['not_normal']:
                stats['structures with a job that did not end normally'] += 1
            f.write(json.dumps(rec) + '\n')
            n += 1
    print('structures with averaged shieldings:', n)
    for k, v in sorted(skipped.items()):
        print('   job steps left out (%s): %d' % (k, v))
    for k, v in sorted(stats.items()):
        print('   %s: %d' % (k, v))


if __name__ == '__main__':
    main(*sys.argv[1:3])
