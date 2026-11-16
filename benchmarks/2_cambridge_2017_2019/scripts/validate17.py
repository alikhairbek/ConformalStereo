"""Which structure of the deposit is behind each isomer of each PyDP4 output, and do my Boltzmann-averaged
shieldings reproduce the 'calc' column?  shift = (reference - shielding) / (1 - reference * 1e-6), so (calc + shielding) / (1 + calc * 1e-6)
must be one constant per element; its spread over the listed atoms is the test (the outputs print two decimals).
Usage: validate17.py <res/core> <structures.jsonl.gz> <runs.json> <out.json>
"""
import collections
import gzip
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dep17
import desc17
import parse17


def num(label):
    return int(re.sub(r'\D', '', label)) - 1


def main(core, structures, runs, out):
    S = dep17.Small(core)
    by_folder = collections.defaultdict(dict)
    with gzip.open(structures, 'rt') as f:
        for line in f:
            r = json.loads(line)
            if r['sigma']:
                by_folder[(r['archive'], r['folder'])][r['stem']] = r
    R = json.load(open(runs))
    cache = {}
    res = []
    for run in R:
        key = (run['md5'], run['archive'], run['folder'])
        if key in cache:
            res.append(dict(cache[key], path=run['path']))
            continue
        d = parse17.parse(S.read(run['path']))
        src = run['path'][:-4]
        desc = desc17.parse(S.read(src).splitlines() if src in S.index else d['desc'])
        group = {}
        for g in desc['groups']:
            for lab in g:
                group[lab] = g
        cands = by_folder.get((run['archive'], run['folder']), {})
        maps = []
        for q in sorted(S.index):
            t = dep17.split(q)
            if t and t[0] == run['archive'] and t[1] == run['folder'] and 'Renum' in t[2] and not t[2].endswith('~'):
                maps = [[int(x) - 1 for x in line.split(',')] for line in S.read(q).splitlines() if line.strip()]
        isomers = []
        for k in sorted(d['isomers']):
            tabs = d['isomers'][k]
            best = None
            for stem, s in cands.items():
                variants = [('own energies', s['sigma'])]
                if s.get('sigma_e1'):
                    variants.append(('energies of level 1', s['sigma_e1']))
                if maps and 2 <= k <= len(maps) + 1 and len(maps[k - 2]) == len(s['sigma']):
                    mp = maps[k - 2]
                    extra = []
                    for name, sig in variants:
                        fwd = [sig[mp[i]] for i in range(len(sig))]
                        inv = [0.0] * len(sig)
                        for i in range(len(sig)):
                            inv[mp[i]] = sig[i]
                        extra.append((name + ', renumbered', fwd, [s['elements'][mp[i]] for i in range(len(sig))]))
                        el_inv = [''] * len(sig)
                        for i in range(len(sig)):
                            el_inv[mp[i]] = s['elements'][i]
                        extra.append((name + ', renumbered (inverse)', inv, el_inv))
                    variants = [(n, sg, s['elements']) for n, sg in variants] + extra
                else:
                    variants = [(n, sg, s['elements']) for n, sg in variants]
                for vname, sig, els in variants:
                    spread, const, ok = {}, {}, True
                    for el in 'CH':
                        off = []
                        for row in tabs[el]:
                            atoms = [num(x) for x in group.get(row[0], [row[0]])]
                            if any(a >= len(sig) or els[a] != el for a in atoms):
                                ok = False
                                break
                            off.append((row[1] + sum(sig[a] for a in atoms) / len(atoms)) / (1 + row[1] * 1e-6))
                        if not ok:
                            break
                        if off:
                            off.sort()
                            spread[el] = off[-1] - off[0]
                            const[el] = off[len(off) // 2]
                    if not ok or not spread:
                        continue
                    score = max(spread.get('C', 0.0) / 10.0, spread.get('H', 0.0))
                    if best is None or score < best[0] - 1e-9:
                        best = (score, stem, spread, const, vname)
            if best is None:
                isomers.append(dict(isomer=k, stem=None))
            else:
                isomers.append(dict(isomer=k, stem=best[1], spread=best[2], const=best[3], variant=best[4]))
        rec = dict(isomers=isomers, n_candidates=len(cands), description=src in S.index)
        cache[key] = rec
        res.append(dict(rec, path=run['path']))
    json.dump(res, open(out, 'w'))
    tot = sum(len(r['isomers']) for r in res)
    good = sum(1 for r in res for i in r['isomers'] if i['stem'] and i['spread'].get('C', 0) <= 0.0111 and i['spread'].get('H', 0) <= 0.0111)
    none = sum(1 for r in res for i in r['isomers'] if not i['stem'])
    full = sum(1 for r in res if r['isomers'] and all(i['stem'] and i['spread'].get('C', 0) <= 0.0111 and i['spread'].get('H', 0) <= 0.0111 for i in r['isomers']))
    print('outputs %d, isomer blocks %d: reproduced to the printed precision %d, no structure found %d' % (len(res), tot, good, none))
    print('outputs with every isomer reproduced: %d' % full)


if __name__ == '__main__':
    main(*sys.argv[1:5])
