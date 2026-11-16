"""DP4 for every workflow that the 2019 deposit allows.

A workflow = kind of geometry (MMFF, PM6, PM7, HF, B3LYP, M06-2X) x NMR level (B3LYP/6-31G(d,p), mPW1PW91/6-311G(d),
M06-2X/6-31G(d,p)) x energies for the Boltzmann weights: those of the NMR level itself, or M06-2X/6-31G(d,p) on the
same geometries (the combination the paper recommends, and the one behind many output files of the deposit).
Usage (with CS_DEPOSIT=conf2018): workflows19.py <res/core> <experiments.json> <energy_source.json> <structures.jsonl.gz> <out.csv>
"""
import csv
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core17
import dep17
import raw17


def main(core, experiments, source, structures, out):
    S = dep17.Small(core)
    E = json.load(open(experiments))
    const = core17.reference_constants(json.load(open(source)), core17.solvents_of(structures))
    data = raw17.load(core, {e['archive'] for e in E})
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['archive', 'experiment', 'geometry', 'nmr', 'energies', 'shielding_level', 'energy_level', 'candidates', 'isomer',
                'structure', 'conformers', 'p_C', 'p_H', 'p_all', 'rms_C', 'rms_H'])
    done = missing = 0
    for e in E:
        for level in dep17.LEVELS:
            if level not in e['levels']:
                continue
            geo = re.sub(r'^BF\d', '', level)
            pairs = [(level, 'same level')]
            if not level.startswith('BF7') and ('BF7' + geo) in e['levels']:
                pairs.append(('BF7' + geo, 'M06-2X/6-31G(d,p)'))
            for elevel, ename in pairs:
                desc, rows = core17.candidates(S, data, e, level, elevel, const)
                if rows is None:
                    missing += 1
                    continue
                logs = {el: [r['log'][el] for r in rows] for el in 'CH'}
                pc, ph = core17.normalise(logs['C']), core17.normalise(logs['H'])
                pa = core17.normalise([a + b for a, b in zip(logs['C'], logs['H'])])
                for k, r in enumerate(rows):
                    err = {el: [p[3] - p[4] for p in r[el]] for el in 'CH'}
                    w.writerow([e['archive'], e['name'], dep17.GEOMETRY[level], dep17.NMR[level], ename, level, elevel, len(rows), k + 1,
                                r['stem'], r['n_conf'], '%.6g' % pc[k], '%.6g' % ph[k], '%.6g' % pa[k],
                                '%.3f' % math.sqrt(sum(x * x for x in err['C']) / len(err['C'])),
                                '%.4f' % math.sqrt(sum(x * x for x in err['H']) / len(err['H']))])
                done += 1
    print('workflow records written: %d; left out because raw outputs are missing for a candidate: %d' % (done, missing))


if __name__ == '__main__':
    main(*sys.argv[1:6])
