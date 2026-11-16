"""Actual structures for the experiments of the 2019 deposit (doi 10.17863/CAM.32463).

The deposit holds 25 of the experiments of the 2017 deposit again, now with six kinds of geometry. The actual
structure of each is taken from the 2017 benchmark (identified there from the paper's SI), after two checks that the
candidates are the same structures in the same order: the names of the structure files, and the Boltzmann-averaged
shieldings at B3LYP/6-31G(d,p) on MMFF geometries, the one level both deposits have.
Usage: truth19.py <dep17 tables dir> <dep19 tables dir> <out.csv>
"""
import collections
import csv
import json
import os
import sys


def shieldings(path, level):
    out = collections.defaultdict(dict)
    for r in csv.DictReader(open(path)):
        if r['level'] == level:
            out[(r['experiment'], int(r['isomer']))][(r['nucleus'], r['label'])] = float(r['shielding'])
    return out


def main(d17, d19, out):
    T = {r['experiment']: r for r in csv.DictReader(open(os.path.join(d17, 'dep17_truth.csv')))}
    E17 = {e['name']: e for e in json.load(open(os.path.join(d17, 'dep17_experiments.json')))}
    E19 = json.load(open(os.path.join(d19, 'dep19_experiments.json')))
    s17 = shieldings(os.path.join(d17, 'dep17_shifts.csv'), 'BF1')
    s19 = shieldings(os.path.join(d19, 'dep19_shifts.csv'), 'BF1')
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['si_compound', 'experiment', 'actual_isomer', 'structure', 'candidates', 'levels', 'same_structure_names',
                'mmff_level_in_both', 'atoms_compared', 'largest_shielding_difference_ppm'])
    n_ok = n_names = n_ident = 0
    for e in E19:
        name = e['name']
        if name not in T:
            print('no actual structure known for', name)
            continue
        t = T[name]
        same = e['stems'] == E17[name]['stems']
        diffs = []
        for k in range(1, len(e['stems']) + 1):
            a, b = s17.get((name, k)), s19.get((name, k))
            if a and b:
                diffs += [abs(a[q] - b[q]) for q in a if q in b]
        worst = max(diffs) if diffs else None
        w.writerow([t['si_compound'], name, t['actual_isomer'], t['structure'], len(e['stems']), len(e['levels']), int(same),
                    int(bool(diffs)), len(diffs), '' if worst is None else '%.4f' % worst])
        n_ok += 1
        n_names += same
        n_ident += bool(diffs) and worst <= 0.0002
        print('%-26s isomer %2s of %2d (%-14s) names %s | MMFF level: %s' % (
            name, t['actual_isomer'], len(e['stems']), t['structure'], 'same' if same else 'DIFFERENT',
            ('%d atoms, largest difference %.4f ppm' % (len(diffs), worst)) if diffs else 'not in the 2019 deposit'))
    print()
    print('experiments %d: actual structure taken over for %d; same structure names in the same order %d; '
          'averaged MMFF-level shieldings identical in both deposits %d' % (len(E19), n_ok, n_names, n_ident))


if __name__ == '__main__':
    main(*sys.argv[1:4])
