"""The experiments of the 2017 deposit: one experimental description with shift data only, run at all twelve levels.
Side runs are left out (J couplings, test folders, outputs written under a one-letter name). The structure behind
each isomer number is the one whose Boltzmann-averaged shieldings reproduce the output (validate17.py); it must be
the same at every level.
Usage: experiments17.py <runs.json> <validate.json> <out.json>
"""
import collections
import json
import re
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dep17


def good(i):
    s = i.get('spread', {})
    return bool(i.get('stem')) and s.get('C', 0) <= 0.0111 and s.get('H', 0) <= 0.0111


def main(runs, validate, out):
    V = {r['path']: r for r in json.load(open(validate))}
    R = json.load(open(runs))
    exps = collections.defaultdict(dict)
    left = collections.Counter()
    for r in R:
        if r['sub']:
            left['in a sub-folder'] += 1
        elif r['level'] == 'J' or r['hasJ'] or re.search(r'NMRj\w*$', r['name']):
            left['with J couplings'] += 1
        elif len(r['name']) < 4 or r['name'] == 'Solifenacin2' or r['name'].endswith('NMRb') or ',' in r['name']:
            left['second output of the same run under another name'] += 1
        elif r['base'].endswith(('_test', 'mmo')) or r['level'] not in dep17.LEVELS:
            left['test folder'] += 1
        else:
            exps[(r['archive'], r['base'], r['name'])][r['level']] = r
    res = []
    amb = 0
    for k in sorted(exps):
        lv = exps[k]
        n = max(r['n'] for r in lv.values())
        stems = collections.defaultdict(collections.Counter)
        for r in lv.values():
            for i in V[r['path']]['isomers']:
                if good(i):
                    stems[i['isomer']][i['stem']] += 1
        amb += sum(1 for i in stems if len(stems[i]) > 1)
        cons = [stems[i].most_common(1)[0][0] if stems[i] else None for i in range(1, n + 1)]
        res.append(dict(archive=k[0], base=k[1], name=k[2], stems=cons, levels={l: r['path'] for l, r in lv.items()}))
    json.dump(res, open(out, 'w'))
    print('PyDP4 outputs %d; left out: %s' % (len(R), ', '.join('%s %d' % kv for kv in sorted(left.items()))))
    print('experiments %d, records %d (experiment x level); experiments with all %d levels %d'
          % (len(res), sum(len(e['levels']) for e in res), len(dep17.LEVELS), sum(1 for e in res if len(e['levels']) == len(dep17.LEVELS))))
    print('isomer numbers without a structure: %d; isomer numbers matched by different structures at different levels: %d'
          % (sum(1 for e in res for s in e['stems'] if s is None), amb))
    print('candidates per experiment: ' + ', '.join('%d x %d' % (c, n) for n, c in sorted(collections.Counter(len(e['stems']) for e in res).items())))


if __name__ == '__main__':
    main(*sys.argv[1:4])
