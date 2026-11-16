"""ConformalStereo: the numbers behind the exhibits of the paper that the final-method package does not hold.
Usage: exhibits.py <pooled_candidates.csv> <folder with the final tables> <jbench_candidates.csv> <out folder> [bootstrap=4000] [seed=20261114]

 1. worked example of Figure 1 (case 87a of JOC2015 at B3LYP/6-31G*/gas)
 2. temperature of the stated probabilities by number of shifts (Figure 3)
 3. coverage of every workflow by origin of the cut-off (Figure 5)
 4. misses per family (Figure 6)
 5. proton-proton couplings with the final score (Table S11)
"""
import sys, os, csv, math, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cslib as cs

table, final, jtable, out = sys.argv[1:5]
BOOT = int(sys.argv[5]) if len(sys.argv) > 5 else 4000
SEED = int(sys.argv[6]) if len(sys.argv) > 6 else 20261114
FINAL = 'likelihood ratio per shift'
os.makedirs(out, exist_ok=True)
recs = cs.load(table)
BENCH = list(collections.OrderedDict((r['bench'], 1) for r in recs))
by_wf = collections.OrderedDict()
for r in recs:
    by_wf.setdefault(r['wf'], []).append(r)
frec = {tuple(r[c] for c in ('benchmark', 'method', 'geometry', 'level', 'case')): r for r in csv.DictReader(open(os.path.join(final, 'final_records.csv')))}


def writer(name, header):
    w = csv.writer(open(os.path.join(out, name), 'w', newline=''))
    w.writerow(header)
    return w


# ------------------------------------------------------------------ 1. worked example
print('1. WORKED EXAMPLE: JOC2015, case 87a, DP4+ at B3LYP/6-31G*/gas')
key = ('JOC2015', 'DP4+', 'B3LYP', 'B3LYP/6-31G*/gas', '87a')
r = next(x for x in recs if (x['bench'], x['method'], x['geo'], x['level'], x['case']) == key)
row = frec[key]
w = writer('exhibit_example.csv', ['candidate', 'is_actual', 'stated_probability', 'log_likelihood_below_best', 'distance_per_shift', 'in_90_set', 'in_80_set', 'shifts', 'cutoff_90', 'cutoff_80'])
c90, c80 = float(row['cutoff_90']), float(row['cutoff_80'])
print('   shifts %d | cut-off at 90%% %.4f, at 80%% %.4f (calibrated on the other 17 families of the level)' % (r['N'], c90, c80))
for i in np.argsort(-r['ell']):
    s = (r['ell'].max() - r['ell'][i]) / r['N']
    w.writerow([r['names'][i], int(i == r['k']), '%.6g' % r['p'][i], '%.8f' % (r['ell'].max() - r['ell'][i]), '%.8f' % s, int(s <= c90), int(s <= c80), r['N'], '%.8f' % c90, '%.8f' % c80])
    print('   %-4s %s stated %.4f%% | log-likelihood below best %7.2f | per shift %.4f | in the 90%% set: %s, in the 80%% set: %s' % (
        r['names'][i], 'actual' if i == r['k'] else '      ', 100 * r['p'][i], r['ell'].max() - r['ell'][i], s, 'yes' if s <= c90 else 'no', 'yes' if s <= c80 else 'no'))

# ------------------------------------------------------------------ 2. temperature by number of shifts
print()
print('2. TEMPERATURE BY NUMBER OF SHIFTS: families of a benchmark in three groups of equal size by their number of shifts')
w = writer('exhibit_temperature_bins.csv', ['benchmark', 'group', 'families', 'records', 'median_shifts', 'lowest_shifts', 'highest_shifts', 'temperature', 'temperature_of_benchmark', 'ratio'])
for b in BENCH:
    rr = [r for r in recs if r['bench'] == b]
    famN = collections.defaultdict(list)
    for r in rr:
        famN[r['family']].append(r['N'])
    famN = {f: float(np.median(v)) for f, v in famN.items()}
    order = sorted(famN, key=lambda f: (famN[f], f))
    T_all = cs.Fit(rr).temperature()
    cells = []
    for g in range(3):
        fams = set(order[g * len(order) // 3:(g + 1) * len(order) // 3])
        sub = [r for r in rr if r['family'] in fams]
        T = cs.Fit(sub).temperature()
        ns = [famN[f] for f in fams]
        w.writerow([b, g + 1, len(fams), len(sub), '%.1f' % float(np.median(ns)), '%.0f' % min(ns), '%.0f' % max(ns), '%.8f' % T, '%.8f' % T_all, '%.8f' % (T / T_all)])
        cells.append('%2d families, %4.0f shifts (%2.0f to %2.0f): T %6.2f' % (len(fams), float(np.median(ns)), min(ns), max(ns), T))
    print('   %-8s T of all %5.2f | %s' % (b, T_all, ' | '.join(cells)))

# ------------------------------------------------------------------ 3. coverage of every workflow by origin of the cut-off
print()
print('3. COVERAGE OF EVERY WORKFLOW BY ORIGIN OF THE CUT-OFF (the family of the record always left out)')
w = writer('exhibit_workflow_coverage.csv', ['target', 'cutoff_from', 'benchmark', 'method', 'geometry', 'level', 'families', 'records', 'coverage', 'mean_set_size'])
ORIGINS = {'own workflow': lambda r, x: x['wf'] == r['wf'], 'same kind of calculation': lambda r, x: x['group'] == r['group'], 'other benchmarks': lambda r, x: x['bench'] != r['bench']}
for target in (0.90, 0.80):
    for origin, same in ORIGINS.items():
        cache, cov, size = {}, [], []
        for r in recs:
            k = (r['wf'] if origin == 'own workflow' else r['group'] if origin == 'same kind of calculation' else r['bench'], r['family'])
            if k not in cache:
                cache[k] = cs.threshold([(x['t'][FINAL], x['family']) for x in recs if same(r, x) and x['family'] != r['family']], target)
            c, z = cs.apply_set(r, FINAL, cache[k])
            cov.append(c)
            size.append(z)
        pos = {id(r): i for i, r in enumerate(recs)}
        per = []
        for wf, rr in by_wf.items():
            idx = [pos[id(r)] for r in rr]
            x = cs.summarize(rr, [cov[i] for i in idx], [size[i] for i in idx])
            per.append(x['coverage'])
            w.writerow([target, origin] + list(wf) + [x['families'], x['records'], '%.8f' % x['coverage'], '%.8f' % x['size']])
        per = np.array(per)
        print('   %d%% %-26s workflows %d | lowest %.1f%%, median %.1f%%, highest %.1f%% | more than 5 points below target: %d' % (
            100 * target, origin, len(per), 100 * per.min(), 100 * np.median(per), 100 * per.max(), int((per < target - 0.05).sum())))

# ------------------------------------------------------------------ 4. misses per family
print()
print('4. MISSES PER FAMILY (sets with one family left out)')
w = writer('exhibit_family_misses.csv', ['benchmark', 'family', 'candidates', 'median_shifts', 'records', 'top_candidate_right', 'missed_90', 'missed_80', 'share_missed_90', 'share_missed_80'])
fam = collections.OrderedDict()
for k, row in frec.items():
    fam.setdefault((row['benchmark'], row['family']), []).append(row)
for (b, f), rr in fam.items():
    n = len(rr)
    m90, m80 = sum(1 - int(x['covered_90']) for x in rr), sum(1 - int(x['covered_80']) for x in rr)
    w.writerow([b, f, rr[0]['candidates'], '%.0f' % float(np.median([int(x['shifts']) for x in rr])), n, sum(int(x['top_right']) for x in rr), m90, m80, '%.8f' % (m90 / n), '%.8f' % (m80 / n)])
for b in BENCH:
    rows = [(f, rr) for (bb, f), rr in fam.items() if bb == b]
    never = sum(1 for f, rr in rows if all(int(x['covered_80']) for x in rr))
    half = sum(1 for f, rr in rows if sum(1 - int(x['covered_80']) for x in rr) >= len(rr) / 2)
    print('   %-8s families %2d | never missed at 80%%: %2d | missed in half their records or more at 80%%: %2d' % (b, len(rows), never, half))

# ------------------------------------------------------------------ 5. couplings with the final score
print()
print('5. PROTON-PROTON COUPLINGS WITH THE FINAL SCORE (33 compounds of the 2017 deposit, 22 families, one level)')
ex = collections.OrderedDict()
for row in csv.DictReader(open(jtable)):
    ex.setdefault(row['experiment'], []).append(row)
MODELS = collections.OrderedDict([
    ('shifts', (lambda r: float(r['logL_C']) + float(r['logL_H']), False)),
    ('shifts and couplings, Gaussian of the PyDP4 source', (lambda r: float(r['logL_C']) + float(r['logL_H']) + float(r['logL_J_gauss']), True)),
    ('shifts and couplings, Fermi-contact part', (lambda r: float(r['logL_C']) + float(r['logL_H']) + float(r['logL_J_fc']), True)),
    ('shifts and couplings, recovered density', (lambda r: float(r['logL_C']) + float(r['logL_H']) + float(r['logL_J_rec']), True)),
    ('shifts, energies averaged', (lambda r: float(r['E4_CH']), False)),
    ('shifts and couplings (Gaussian), energies averaged', (lambda r: float(r['E4_CHJ_gauss']), True)),
    ('shifts and couplings (recovered), energies averaged', (lambda r: float(r['E4_CHJ_rec']), True))])
w = writer('exhibit_couplings.csv', ['data', 'actual_first', 'compounds', 'target', 'coverage', 'mean_set_size', 'size_lo', 'size_hi', 'single_structure_sets', 'compared_with', 'difference_of_mean_set_size', 'diff_lo', 'diff_hi'])
base = {}
for name, (fn, with_j) in MODELS.items():
    rr = []
    for e, rows in ex.items():
        ell = np.array([fn(r) for r in rows])
        ell = ell - ell.max()
        n = int(rows[0]['n_C']) + int(rows[0]['n_H']) + (int(rows[0]['couplings']) if with_j else 0)
        k = [int(r['is_actual']) for r in rows].index(1)
        s = (ell.max() - ell) / n
        top = int(np.argmax(ell))
        rr.append(dict(wf=('J',), family=rows[0]['archive'], K=len(rows), N=n, first=(k == top), t={FINAL: float(s[k])}, s={FINAL: np.sort(s)}, tops={FINAL: float(s[top])}))
    line = []
    for target in (0.90, 0.80):
        cov, size = [], []
        for r in rr:
            thr = cs.threshold([(x['t'][FINAL], x['family']) for x in rr if x['family'] != r['family']], target)
            c, z = cs.apply_set(r, FINAL, thr)
            cov.append(c)
            size.append(z)
        x = cs.summarize(rr, cov, size, boot=BOOT, seed=SEED + 11)
        ref = 'shifts, energies averaged' if name.endswith('energies averaged') and name != 'shifts, energies averaged' else 'shifts'
        if name in ('shifts', 'shifts, energies averaged'):
            base[(name, target)] = size
        if name == 'shifts':
            d = (0.0, 0.0, 0.0)
        else:
            um = cs.unit_means(rr, [(float(a) - float(b),) for a, b in zip(size, base[(ref, target)])])
            A = np.array([v[0] for v in um.values()])
            rg = np.random.default_rng(SEED + 12)
            bt = A[rg.integers(0, len(A), size=(BOOT, len(A)))].mean(axis=1)
            d = (float(A.mean()), float(np.quantile(bt, 0.025)), float(np.quantile(bt, 0.975)))
        w.writerow([name, sum(r['first'] for r in rr), len(rr), target, '%.8f' % x['coverage'], '%.8f' % x['size'], '%.8f' % x['size_lo'], '%.8f' % x['size_hi'], '%.8f' % x['single'],
                    ref, '%.8f' % d[0], '%.8f' % d[1], '%.8f' % d[2]])
        line.append('%d%%: %.2f (%.2f to %.2f), single %2.0f%%, against %s %+.2f (%+.2f to %+.2f)' % (100 * target, x['size'], x['size_lo'], x['size_hi'], 100 * x['single'],
                                                                                              'the same without couplings' if ref != 'shifts' or name != 'shifts, energies averaged' else 'plain shifts', d[0], d[1], d[2]))
    print('   %-52s actual first %2d of %d | %s' % (name, sum(r['first'] for r in rr), len(rr), ' | '.join(line)))

# ------------------------------------------------------------------ 6. sources and reproduction (Table S1; facts from the logs of the earlier packages)
SOURCES = [
    ['JOC2015', 'Grimblat, Zanardi and Sarotti, J. Org. Chem. 2015, 80, 12526', 'Supporting Information of the paper: shieldings of every candidate of 48 cases at 24 NMR levels on B3LYP/6-31G* geometries',
     'DP4+ with the parameters of the paper for each level', 'every candidate within 0.5 point of the published probability in 1,040 of 1,152 records', '1,152 of 1,152'],
    ['STX150', 'Luciano et al., J. Chem. Inf. Model. 2026, 66, 11062', 'Supporting Information of the paper: shieldings of the 110 validation cases at six levels (MMFF, xTB and B3LYP geometries)',
     'MM-DP4+, xTB-DP4+ and DP4+ with the parameters of the paper', '505 records verified against the published table (481 within 1 point, 24 within the rounding of the printed shieldings), 42 close, 12 corrected',
     '559 of 660; not used: 63 unverified, 38 excluded (no experimental data printed, a copied shielding column, or too few nuclei)'],
    ['DP4-AI', 'Howarth, Ermanis and Goodman, Chem. Sci. 2020, 11, 4351', 'Data deposit doi 10.17863/CAM.47721: Gaussian outputs of every conformer of 46 molecules, 344 candidates, four protocols',
     'DP4 as in the 2020 code of PyDP4: one Gaussian per nucleus (2.269 ppm for 13C, 0.187 ppm for 1H), internal scaling, unassigned shift lists paired by order',
     '109 of the 131 stored DP4 outputs reproduced to 0.01 ppm in the shifts, 107 of them to 0.05 point in the probabilities', '182 of 184 (the MMFF protocol exists for 44 molecules); 141 without a data flag'],
    ['Cam2017', 'Ermanis, Parkes, Agback and Goodman, Org. Biomol. Chem. 2017, 15, 8998', 'Data deposit doi 10.17863/CAM.13222: 69,540 NMR jobs of 42 compounds at twelve NMR levels on MMFF geometries',
     'DP4 with one Gaussian per nucleus (2.269 ppm for 13C, 0.187 ppm for 1H) on internally scaled shifts, atoms and pairing as in PyDP4 of 2015 to 2017',
     'probabilities within 0.1 point of every printed value in all 455 outputs that are comparable in full; the four per-compound tables of section S-4 of the SI reproduced', '504 of 504'],
    ['Cam2019', 'data deposit of Ermanis and co-workers, doi 10.17863/CAM.32463', '25 of the 42 experiments on MMFF, PM6, PM7, HF, B3LYP and M06-2X geometries, three NMR levels, two choices of conformer energies: 30 workflows',
     'as for Cam2017', 'probabilities within 0.1 point of every printed value in all 207 outputs that are comparable in full; 202 of 2,412 isomer blocks of the deposit not reproduced',
     '703 of 750 (three experiments exist on HF and DFT geometries only; two records on B3LYP geometries are missing)']]
w = writer('exhibit_sources.csv', ['benchmark', 'paper', 'data', 'recomputed_as', 'agreement_with_the_published_numbers', 'records_used'])
for row in SOURCES:
    w.writerow(row)
print()
print('6. SOURCES AND REPRODUCTION: %d benchmarks written to exhibit_sources.csv' % len(SOURCES))

# ------------------------------------------------------------------ 7. kinds of sets per benchmark (Figure 4), from the sets of every record
print()
print('7. KINDS OF SETS PER BENCHMARK (shares averaged per family within a workflow, then over these units)')
w = writer('exhibit_set_kinds.csv', ['benchmark', 'target', 'one_structure', 'several_not_all', 'every_candidate', 'mean_set_size', 'mean_candidates'])
for b in BENCH:
    rows = [row for row in frec.values() if row['benchmark'] == b]
    for t in ('90', '80'):
        unit = collections.defaultdict(list)
        for row in rows:
            z, K = int(row['set_size_' + t]), int(row['candidates'])
            unit[(row['method'], row['geometry'], row['level'], row['family'])].append((float(z == 1), float(1 < z < K), float(z == K and K > 1), float(z), float(K)))
        m = [sum(sum(v[i] for v in vs) / len(vs) for vs in unit.values()) / len(unit) for i in range(5)]
        w.writerow([b, t] + ['%.6f' % x for x in m[:3]] + ['%.8f' % m[3], '%.8f' % m[4]])
        print('   %-8s %s%%: one structure %5.1f%% | several, not all %5.1f%% | every candidate %5.1f%% | mean size %.2f of %.2f' % (b, t, 100 * m[0], 100 * m[1], 100 * m[2], m[3], m[4]))
