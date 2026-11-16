"""DP4 for the 2017 deposit, recomputed from the raw Gaussian outputs.

For every experiment (one experimental description, 2-32 candidate structures) and each of the twelve NMR levels:
  shieldings of each candidate = Boltzmann average over its conformers, weights from the SCF energies of the same level
  atoms as in PyDP4: every C and H, atoms after OMIT left out, each group of equivalent atoms averaged (first label kept)
  pairing with the experimental shifts as in PyDP4 (2015-2017): assigned labels first, then limited sets, the rest in
  ascending order
  DP4 = single Gaussian error model on internally scaled shifts (sigma 2.269 ppm for 13C, 0.187 ppm for 1H)
and the comparison with the output file kept in the deposit for the same experiment and level.
Usage: dp4_17.py <res/core> <experiments.json> <energy_source.json> <structures.jsonl.gz> <out dir>
"""
import collections
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core17
import dep17
import parse17
import raw17
from core17 import LEVELS, LEVEL_NAMES, normalise


def main(core, experiments, source, structures, out):
    S = dep17.Small(core)
    E = json.load(open(experiments))
    SRC = json.load(open(source))
    src = {(r['archive'], r['name'], r['level'], r['isomer']): r for r in SRC}
    const = core17.reference_constants(SRC, core17.solvents_of(structures))
    data = raw17.load(core, {e['archive'] for e in E})
    os.makedirs(out, exist_ok=True)
    fc = open(os.path.join(out, dep17.PREFIX + '_candidates.csv'), 'w', newline='')
    wc = csv.writer(fc)
    wc.writerow(['archive', 'experiment', 'level', 'isomer', 'structure', 'conformers', 'n_C', 'n_H',
                 'p_C', 'p_H', 'p_all', 'rms_C', 'max_C', 'rms_H', 'max_H',
                 'deposit_p_C', 'deposit_p_H', 'deposit_p_all', 'deposit_shifts'])
    fs = open(os.path.join(out, dep17.PREFIX + '_shifts.csv'), 'w', newline='')
    ws = csv.writer(fs)
    ws.writerow(['archive', 'experiment', 'level', 'isomer', 'nucleus', 'label', 'shielding', 'shift_unscaled', 'shift_scaled', 'exp'])
    fr = open(os.path.join(out, dep17.PREFIX + '_records.csv'), 'w', newline='')
    wr = csv.writer(fr)
    wr.writerow(['archive', 'experiment', 'level', 'level_name', 'candidates', 'n_C', 'n_H', 'top_all', 'p_top_all', 'top_C', 'top_H',
                 'deposit_top_all', 'deposit_model', 'deposit_energies', 'max_abs_diff_all_points', 'pairing_same_as_deposit'])
    fk = open(os.path.join(out, dep17.PREFIX + '_reference_constants.csv'), 'w', newline='')
    wk = csv.writer(fk)
    solv = core17.solvents_of(structures)
    wk.writerow(['archive', 'folder', 'solvent', 'level', 'level_name', 'reference_C', 'reference_H'])
    for (a, b, l) in sorted({k[:3] for k in const}, key=lambda x: (x[0], x[1], dep17.level_index(x[2]))):
        wk.writerow([a, b, solv.get((a, b), ''), l, LEVEL_NAMES[l], const.get((a, b, l, 'C'), ''), const.get((a, b, l, 'H'), '')])
    fk.close()
    stats = collections.Counter()
    worst = []
    incomplete = []
    for e in E:
        for level in LEVELS:
            path = e['levels'].get(level)
            if not path:
                continue
            dep = parse17.parse(S.read(path))
            desc, rows = core17.candidates(S, data, e, level, level, const)
            if rows is None:
                incomplete.append((e['name'], level))
                continue
            n = len(rows)
            same_pairing = True
            for row in rows:
                k = row['isomer']
                for el in 'CH':
                    for lab, sh, u, sc, x in row[el]:
                        ws.writerow([e['archive'], e['name'], level, k, el, lab, '%.4f' % sh, '%.3f' % u if u is not None else '', '%.3f' % sc, x])
                    if sorted((r[0], r[3]) for r in dep['isomers'][k][el]) != sorted((p[0], p[4]) for p in row[el]):
                        same_pairing = False
            logs = {el: [row['log'][el] for row in rows] for el in 'CH'}
            pc, ph = normalise(logs['C']), normalise(logs['H'])
            pa = normalise([a + b for a, b in zip(logs['C'], logs['H'])])
            dres = dep['results']
            d_all = dres.get('CH') or dres.get('all')
            model_default = all(abs(math.log10(max(dep['isomers'][k]['p'].get('C', 1e-300), 1e-300)) - logs['C'][k - 1] / math.log(10)) < 0.2
                                for k in range(1, n + 1)) if all('C' in dep['isomers'][k]['p'] for k in range(1, n + 1)) else False
            srcs = [src[(e['archive'], e['name'], level, k)] for k in range(1, n + 1)]
            if all(level in s['reproduced_with'] for s in srcs):
                energies = 'same level'
            elif all(s['reproduced_with'] for s in srcs):
                other = collections.Counter(l for s in srcs if level not in s['reproduced_with'] for l in s['reproduced_with'])
                energies = 'other level (%s)' % other.most_common(1)[0][0]
            else:
                energies = 'not reproduced'
            diff = max(abs(100 * a - b) for a, b in zip(pa, d_all)) if d_all and len(d_all) == n else float('nan')
            top = max(range(n), key=lambda i: pa[i])
            wr.writerow([e['archive'], e['name'], level, LEVEL_NAMES[level], n, len(rows[0]['C']), len(rows[0]['H']),
                         top + 1, '%.4f' % pa[top], max(range(n), key=lambda i: pc[i]) + 1, max(range(n), key=lambda i: ph[i]) + 1,
                         (max(range(n), key=lambda i: d_all[i]) + 1) if d_all else '', 'default' if model_default else 'other',
                         energies, '%.1f' % diff, int(same_pairing)])
            for k, row in enumerate(rows):
                s = srcs[k]
                err = {el: [p[3] - p[4] for p in row[el]] for el in 'CH'}
                wc.writerow([e['archive'], e['name'], level, k + 1, row['stem'], row['n_conf'], len(row['C']), len(row['H']),
                             '%.6g' % pc[k], '%.6g' % ph[k], '%.6g' % pa[k],
                             '%.3f' % math.sqrt(sum(x * x for x in err['C']) / len(err['C'])), '%.2f' % max(abs(x) for x in err['C']),
                             '%.4f' % math.sqrt(sum(x * x for x in err['H']) / len(err['H'])), '%.3f' % max(abs(x) for x in err['H']),
                             dres['C'][k] if 'C' in dres else '', dres['H'][k] if 'H' in dres else '', d_all[k] if d_all else '',
                             'same level' if level in s['reproduced_with'] else (','.join(s['reproduced_with']) or 'not reproduced')])
            stats['records'] += 1
            stats['same first isomer'] += bool(d_all) and (max(range(n), key=lambda i: d_all[i]) == top)
            stats['pairing same as deposit'] += same_pairing
            key = (energies == 'same level', model_default)
            stats[key] += 1
            if key == (True, True):
                stats['comparable'] += 1
                if diff <= 0.1001:
                    stats['comparable, within 0.1 point'] += 1
                else:
                    worst.append((diff, e['name'], level))
    for f in (fc, fs, fr):
        f.close()
    print('records (experiment x level): %d' % stats['records'])
    if incomplete:
        print('left out, Gaussian outputs missing for at least one candidate: %d (%s)'
              % (len(incomplete), ', '.join('%s %s' % x for x in incomplete)))
    print('pairing of shifts and atoms identical to the deposit output: %d' % stats['pairing same as deposit'])
    print('deposit output with the same energies and the default error model: %d; my probabilities within 0.1 point of every printed value: %d'
          % (stats['comparable'], stats['comparable, within 0.1 point']))
    print('deposit output with the same energies, another error model: %d' % stats[(True, False)])
    print('first-ranked isomer the same as in the deposit output: %d of %d' % (stats['same first isomer'], stats['records']))
    print('deposit output with other energies: default model %d, other model %d' % (stats[(False, True)], stats[(False, False)]))
    for d, nme, l in sorted(worst, reverse=True)[:12]:
        print('   differs by %.1f points: %s %s' % (d, nme, l))


if __name__ == '__main__':
    main(*sys.argv[1:6])
