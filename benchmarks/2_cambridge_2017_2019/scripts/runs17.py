"""Table of all PyDP4 outputs (*.dp4) of the 2017 deposit.
Usage: runs17.py <res/core> <out.json>
"""
import json
import re
import sys

import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dep17
import parse17


level_of = dep17.level_of


def main(core, out):
    S = dep17.Small(core)
    runs = []
    for p in S.paths():
        if not p.endswith('.dp4'):
            continue
        s = dep17.split(p)
        if s is None:
            arch, folder, name = '(top)', '', p.split('/')[-1]
        else:
            arch, folder, name = s
        d = parse17.parse(S.read(p))
        base, level = level_of(folder) if folder else ('', '')
        iso = d['isomers']
        n = len(iso)
        first = iso[min(iso)] if iso else {'C': [], 'H': []}
        runs.append(dict(path=p, md5=S.index[p][1], archive=arch, folder=folder, base=base, level=level,
                         sub='/'.join(folder.split('/')[1:]), name=name[:-4], n=n,
                         nC=len(first['C']), nH=len(first['H']),
                         hasJ=any(r[5] for v in iso.values() for r in v['H']),
                         results=d['results']))
    json.dump(runs, open(out, 'w'))
    print('outputs', len(runs), 'distinct contents', len({r['md5'] for r in runs}))


if __name__ == '__main__':
    main(*sys.argv[1:3])
