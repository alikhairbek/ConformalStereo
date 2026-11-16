"""For every isomer of every primary PyDP4 output: which energies reproduce the printed shifts?
The shieldings are those of the output's own level; the Boltzmann weights are tried with the energies of each of the
twelve levels. Usage: energy_source17.py <res/core> <experiments.json> <out.json>
"""
import collections
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dep17
import desc17
import parse17
import raw17

LEVELS = dep17.LEVELS
TOL = 0.0111


def num(label):
    return int(re.sub(r'\D', '', label)) - 1


def renum_maps(S, archive, folder):
    for q in sorted(S.index):
        t = dep17.split(q)
        if t and t[0] == archive and t[1] == folder and 'Renum' in t[2] and not t[2].endswith('~'):
            return [[int(x) - 1 for x in line.split(',')] for line in S.read(q).splitlines() if line.strip()]
    return []


def spreads(tabs, group, sig, els):
    out = {}
    for el in 'CH':
        off = []
        for row in tabs[el]:
            atoms = [num(x) for x in group.get(row[0], [row[0]])]
            if any(a >= len(sig) or els[a] != el for a in atoms):
                return None
            off.append((row[1] + sum(sig[a] for a in atoms) / len(atoms)) / (1 + row[1] * 1e-6))
        if off:
            out[el] = (max(off) - min(off), sorted(off)[len(off) // 2])
    return out


def main(core, experiments, out):
    S = dep17.Small(core)
    E = json.load(open(experiments))
    data = raw17.load(core, {e['archive'] for e in E})
    res = []
    for e in E:
        for level in LEVELS:
            path = e['levels'].get(level)
            if not path:
                continue
            d = parse17.parse(S.read(path))
            desc = desc17.parse(S.read(path[:-4]).splitlines())
            group = {lab: g for g in desc['groups'] for lab in g}
            folder = dep17.split(path)[1]
            maps = renum_maps(S, e['archive'], folder)
            for k in sorted(d['isomers']):
                stem = e['stems'][k - 1]
                levels = data.get((e['archive'], e['base'], stem), {})
                fits = {}
                for el_level in LEVELS:
                    a = raw17.average(levels, level, el_level)
                    if a is None:
                        continue
                    sig, els, ks = a
                    if maps and k >= 2 and len(maps[k - 2]) == len(sig):
                        sig, els = raw17.renumber(sig, maps[k - 2]), raw17.renumber(els, maps[k - 2])
                    sp = spreads(d['isomers'][k], group, sig, els)
                    if sp:
                        fits[el_level] = sp
                ok = [l for l in LEVELS if l in fits and all(v[0] <= TOL for v in fits[l].values())]
                copy_of = []
                if not ok:
                    # is the output a copy of the output of another level (shieldings and energies of that level)?
                    for other in LEVELS:
                        if other == level:
                            continue
                        a = raw17.average(levels, other, other)
                        if a is None:
                            continue
                        sig, els, ks = a
                        if maps and k >= 2 and len(maps[k - 2]) == len(sig):
                            sig, els = raw17.renumber(sig, maps[k - 2]), raw17.renumber(els, maps[k - 2])
                        sp = spreads(d['isomers'][k], group, sig, els)
                        if sp and all(v[0] <= TOL for v in sp.values()):
                            copy_of.append(other)
                own = fits.get(level)
                res.append(dict(archive=e['archive'], base=e['base'], name=e['name'], level=level, isomer=k, stem=stem,
                                n_conf=len(levels.get(level, {})), reproduced_with=ok, copy_of=copy_of,
                                own_spread={el: round(v[0], 4) for el, v in own.items()} if own else None,
                                const={el: round(v[1], 4) for el, v in fits[ok[0]].items()} if ok else None))
    json.dump(res, open(out, 'w'))
    tot = len(res)
    own = sum(1 for r in res if r['level'] in r['reproduced_with'])
    other = sum(1 for r in res if r['reproduced_with'] and r['level'] not in r['reproduced_with'])
    print('isomer blocks %d: reproduced with the energies of the same level %d, only with another level %d, not reproduced %d'
          % (tot, own, other, tot - own - other))
    ncp = sum(1 for r in res if r.get('copy_of'))
    if ncp:
        print('of the blocks not reproduced, %d hold the shifts of another level (a copied output file)' % ncp)
    print()
    print('isomer blocks with more than one conformer, by level of the output (shieldings):')
    tab = collections.defaultdict(collections.Counter)
    for r in res:
        if r['n_conf'] < 2:
            continue
        rw = r['reproduced_with']
        if r['level'] in rw:
            k = 'same level'
        elif rw:
            k = 'energies of ' + (rw[0] if len(rw) == 1 else 'another level (several fit)')
        else:
            k = 'not reproduced'
        tab[r['level']][k] += 1
    for l in LEVELS:
        print('   %-5s %s' % (l, ', '.join('%s %d' % kv for kv in sorted(tab[l].items(), key=lambda x: -x[1]))))
    print()
    print('outputs behind which other energies stand (experiment: level of the output <- level of the energies, isomer blocks):')
    oth = collections.Counter((r['name'], r['level'], r['reproduced_with'][0] if len(r['reproduced_with']) == 1 else 'several')
                              for r in res if r['reproduced_with'] and r['level'] not in r['reproduced_with'])
    for (n, l, w), c in sorted(oth.items(), key=lambda x: (x[0][0], dep17.level_index(x[0][1]))):
        print('   %-24s %-5s <- %-8s %d' % (n, l, w, c))
    cp = collections.Counter((r['name'], r['level'], r['copy_of'][0]) for r in res if r.get('copy_of'))
    if cp:
        print()
        print('outputs that hold the shifts of another level (experiment: level of the folder <- level of the shifts, isomer blocks):')
        for (n, l, w), c in sorted(cp.items(), key=lambda x: (x[0][0], dep17.level_index(x[0][1]))):
            print('   %-24s %-6s <- %-6s %d' % (n, l, w, c))
    print()
    print('not reproduced (experiment, level, isomer, structure, conformers, spread of the offsets C / H in ppm):')
    for r in res:
        if not r['reproduced_with'] and not r.get('copy_of'):
            sp = r['own_spread'] or {}
            print('   %-24s %-5s %2d %-16s %3d  %.2f / %.3f' % (r['name'], r['level'], r['isomer'], r['stem'], r['n_conf'], sp.get('C', 0), sp.get('H', 0)))


if __name__ == '__main__':
    main(*sys.argv[1:4])
