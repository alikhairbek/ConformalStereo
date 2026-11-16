"""Quality checks on the raw files behind the benchmark of the 2017 deposit.
  - the same conformers at all twelve levels, and the same geometry (principal moments of inertia printed by Gaussian)
  - jobs that did not end normally, conformers without shieldings
  - conformers that occur twice within a structure (same moments and energy)
Usage: qc17.py <res/core> <experiments.json> <dep17_truth.csv> <out.csv>
"""
import collections
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dep17
from runs17 import level_of

LEVELS = ['BF%d' % i for i in range(1, 13)]
NAME = re.compile(r'^(?P<stem>.+?)ginp(?P<conf>\d+)\.out$')
HARTREE = 2625.499629554010


def main(core, experiments, truth, out):
    keep = {r['experiment'] for r in csv.DictReader(open(truth))}
    E = [e for e in json.load(open(experiments)) if e['name'] in keep]
    want = {(e['archive'], e['base'], s) for e in E for s in e['stems']}
    data = collections.defaultdict(lambda: collections.defaultdict(dict))
    for fn in sorted(os.listdir(core)):
        if not re.match(r'gaussian_dp4x2_\d+\.jsonl\.gz$', fn):
            continue
        for r in dep17.jsonl(os.path.join(core, fn)):
            s = dep17.split(r['path'])
            if s is None or '/' in s[1]:
                continue
            m = NAME.match(s[2])
            if not m:
                continue
            base, level = level_of(s[1])
            key = (s[0], base, m.group('stem'))
            if key not in want or level not in LEVELS:
                continue
            d = data[key][level].setdefault(int(m.group('conf')), {'steps': 0})
            d['steps'] += 1
            if 'iso' in r:
                pm = json.loads(r['pm']) if isinstance(r['pm'], str) else r['pm']
                d.update(scf=float(r['scf']), pm=pm, norm=str(r.get('norm')) == '1')
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['archive', 'folder', 'structure', 'conformers', 'levels_with_all_conformers', 'max_moment_deviation_percent',
                'jobs_not_ended_normally', 'conformers_without_shieldings', 'repeated_conformers', 'max_relative_energy_kJ_level1',
                'conformers_above_1_percent_level1'])
    tot = collections.Counter()
    nconf = []
    for key in sorted(want):
        lv = data.get(key, {})
        ref = lv.get('BF1', {})
        confs = sorted(k for k, v in ref.items() if 'scf' in v)
        full = sum(1 for l in LEVELS if sorted(k for k, v in lv.get(l, {}).items() if 'scf' in v) == confs)
        dev = 0.0
        for l in LEVELS:
            for k in confs:
                v = lv.get(l, {}).get(k)
                if v and 'pm' in v:
                    for a, b in zip(v['pm'], ref[k]['pm']):
                        if b:
                            dev = max(dev, abs(a - b) / b * 100)
        bad = sum(1 for l in LEVELS for v in lv.get(l, {}).values() if 'scf' in v and not v['norm'])
        nosh = sum(1 for l in LEVELS for v in lv.get(l, {}).values() if 'scf' not in v)
        seen = collections.Counter((round(ref[k]['scf'], 6), tuple(round(x, 2) for x in ref[k]['pm'])) for k in confs)
        rep = sum(c - 1 for c in seen.values())
        e = [ref[k]['scf'] for k in confs]
        rel = [(x - min(e)) * HARTREE for x in e] if e else []
        import math
        wts = [math.exp(-x * 1000 / (8.3145 * 298.15)) for x in rel]
        z = sum(wts) or 1.0
        w.writerow([key[0], key[1], key[2], len(confs), full, '%.4f' % dev, bad, nosh, rep, '%.1f' % (max(rel) if rel else 0),
                    sum(1 for x in wts if x / z > 0.01)])
        nconf.append(len(confs))
        tot['structures'] += 1
        tot['all twelve levels complete'] += full == 12
        tot['same geometry at all levels (moments within 0.01%)'] += dev <= 0.01
        tot['with a job that did not end normally'] += bad > 0
        tot['with a conformer without shieldings'] += nosh > 0
        tot['with repeated conformers'] += rep > 0
        tot['one conformer only'] += len(confs) == 1
    nconf.sort()
    print('candidate structures of the 42 compounds: %d; conformers per structure: median %d, range %d-%d, total %d (x 12 levels = %d NMR jobs)'
          % (tot['structures'], nconf[len(nconf) // 2], nconf[0], nconf[-1], sum(nconf), 12 * sum(nconf)))
    for k in ('all twelve levels complete', 'same geometry at all levels (moments within 0.01%)', 'with a job that did not end normally',
              'with a conformer without shieldings', 'with repeated conformers', 'one conformer only'):
        print('   %s: %d' % (k, tot[k]))


if __name__ == '__main__':
    main(*sys.argv[1:5])
