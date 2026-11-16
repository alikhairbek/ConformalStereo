"""Which structures, shieldings and energies stand behind every candidate of the PyDP4 outputs with couplings?
Usage: jsource.py <res/core of run 2> <lib> <dep17_experiments.json> <out.json>
For each isomer block of each output in a _J folder, the printed shifts are compared with Boltzmann averages built from
the raw outputs: shieldings of the coupling jobs (level 'J') or of one of the twelve levels, weights from the energies of
the coupling jobs or of one of the twelve levels. A block is reproduced when the offsets of all its carbon and of all
its proton shifts agree within 0.011 ppm."""
import sys, os, json, collections
core, lib, exp_json, out = sys.argv[1:5]
sys.path.insert(0, lib)
import dep17, desc17, parse17, raw17
from energy_source17 import renum_maps, spreads, TOL

S = dep17.Small(core)
E = json.load(open(exp_json))
main = {}
for e in E:
    p = e['levels'].get('BF1')
    iso = parse17.parse(S.read(p))['isomers']
    main[(e['archive'], e['name'])] = dict(e=e, exp=sorted(r[3] for r in iso[min(iso)]['C']))
data = raw17.load(core, {e['archive'] for e in E})
ALL = ['J'] + dep17.LEVELS
res = []
for path in S.paths():
    if '_J/' not in path or not path.endswith('.dp4'):
        continue
    arch, folder, name = dep17.split(path)
    if '/' in folder:
        continue
    d = parse17.parse(S.read(path))
    iso = d['isomers']
    if not iso or not all('J' in iso[k]['p'] for k in iso) or max(sum(len(r[5]) for r in iso[k]['H']) for k in iso) == 0:
        continue
    mine = sorted(r[3] for r in iso[min(iso)]['C'])
    hits = [k for k, m in main.items() if k[0] == arch and len(m['exp']) == len(mine) and max(abs(x - y) for x, y in zip(m['exp'], mine)) < 0.006]
    assert len(hits) == 1, (path, hits)
    e = main[hits[0]]['e']
    desc = desc17.parse(S.read(path[:-4]).splitlines())
    group = {lab: g for g in desc['groups'] for lab in g}
    maps = renum_maps(S, arch, folder)
    jbase = folder[:-2]                                   # the coupling folder may hang off another base folder than the experiment's
    rec = dict(path=path, archive=arch, base=e['base'], jbase=jbase, folder=folder, experiment=e['name'], file=name[:-4], stems=e['stems'], isomers=[])
    for k in sorted(iso):
        fits = []
        for si, stem in enumerate(e['stems'], 1):
            levels = dict(data.get((arch, e['base'], stem), {}))
            if 'J' in data.get((arch, jbase, stem), {}):
                levels['J'] = data[(arch, jbase, stem)]['J']
            for sl in ALL:
                for el in ALL:
                    a = raw17.average(levels, sl, el)
                    if a is None:
                        continue
                    sig, els, ks = a
                    for use_map in ((False, True) if maps else (False,)):
                        s2, e2 = sig, els
                        if use_map:
                            if not (si >= 2 and len(maps[si - 2]) == len(sig)):
                                continue
                            s2, e2 = raw17.renumber(sig, maps[si - 2]), raw17.renumber(els, maps[si - 2])
                        sp = spreads(iso[k], group, s2, e2)
                        if sp and all(v[0] <= TOL for v in sp.values()):
                            fits.append((stem, si, sl, el, use_map, len(ks)))
        rec['isomers'].append(dict(isomer=k, fits=fits))
    res.append(rec)
json.dump(res, open(out, 'w'))
tot = sum(len(r['isomers']) for r in res)
rep = sum(1 for r in res for i in r['isomers'] if i['fits'])
print('outputs with couplings: %d | isomer blocks: %d | reproduced from raw data: %d' % (len(res), tot, rep))
how = collections.Counter()
for r in res:
    for i in r['isomers']:
        if not i['fits']:
            how['not reproduced'] += 1
            continue
        stems = sorted(set(f[0] for f in i['fits']))
        own = [f for f in i['fits'] if f[2] == 'J' and f[3] == 'J']
        sl = sorted(set(f[2] for f in i['fits']))
        how[('one structure' if len(stems) == 1 else '%d structures fit' % len(stems), 'own shieldings and energies' if own else 'shieldings %s, energies %s' % ('/'.join(sl)[:20], '/'.join(sorted(set(f[3] for f in i['fits'])))[:30]))] += 1
for k, v in how.most_common():
    print('   %5d  %s' % (v, k))
order = collections.Counter()
for r in res:
    same = all(any(f[1] == i['isomer'] for f in i['fits']) for i in r['isomers'] if i['fits'])
    order[same] += 1
    bad = [i['isomer'] for i in r['isomers'] if not i['fits']]
    if bad or not same:
        print('   %-26s isomers %2d | not reproduced: %s | order as in the shift-only output: %s | %s' % (r['file'], len(r['isomers']), bad, same,
              [(i['isomer'], sorted(set(f[1] for f in i['fits']))) for i in r['isomers']][:8]))
print('outputs whose candidates stand in the order of the shift-only output:', dict(order))
