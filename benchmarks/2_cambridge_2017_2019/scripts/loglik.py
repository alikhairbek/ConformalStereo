"""Log-likelihoods of every candidate, the raw material for conformal scores.

For the active deposit (CS_DEPOSIT) and every experiment and level: the DP4 log-likelihood of each candidate for 13C and
for 1H (sum over the nuclei of log(2 * Phi(-|error| / sigma)), single-Gaussian model, internally scaled shifts).
Probabilities are their softmax; nothing is lost to rounding or underflow here. In the 2019 deposit each level is
written twice where possible: with the energies of the NMR level and with the M06-2X/6-31G(d,p) energies.
Usage: loglik.py <res/core> <experiments.json> <energy_source.json> <structures.jsonl.gz> <out.csv>
"""
import csv
import json
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
    w.writerow(['archive', 'experiment', 'shielding_level', 'energy_level', 'geometry', 'nmr', 'candidates', 'isomer', 'structure',
                'n_C', 'n_H', 'logL_C', 'logL_H'])
    n = 0
    for e in E:
        for level in dep17.LEVELS:
            if level not in e['levels']:
                continue
            pairs = [level]
            if dep17.KEY == 'conf2018' and not level.startswith('BF7'):
                other = 'BF7' + re.sub(r'^BF\d', '', level)
                if other in e['levels']:
                    pairs.append(other)
            for elevel in pairs:
                desc, rows = core17.candidates(S, data, e, level, elevel, const)
                if rows is None:
                    continue
                for r in rows:
                    w.writerow([e['archive'], e['name'], level, elevel, dep17.GEOMETRY[level], dep17.NMR[level], len(rows), r['isomer'],
                                r['stem'], len(r['C']), len(r['H']), '%.6f' % r['log']['C'], '%.6f' % r['log']['H']])
                n += 1
    print('records written: %d' % n)


if __name__ == '__main__':
    main(*sys.argv[1:6])
