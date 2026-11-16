"""Writes the test files of this folder. Only needed when the fixtures are to be made again (needs pandas and openpyxl;
the tests themselves need nothing but Python). All numbers are made up.

    python make_fixtures.py

dp4plus_app_results.xlsx   written with the same calls as the results file of DP4+App 2.2 (output_module.gen_out_dp4)
dp4plus_sheet.xlsx         a workbook with the cell layout of the DP4+ Excel sheet (row "DP4+", header "Experimental")
pydp4_new.dp4, pydp4_old.dp4   texts in the two layouts of PyDP4 / DP4-AI output
"""
import os, math, random
import numpy as np
import pandas as pd
from openpyxl import Workbook, load_workbook

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')
random.seed(11); np.random.seed(11)

# ------------------------------------------------------------------ DP4+App results
isoms = ['01', '02', '03', '04']
nuc = ['C'] * 14 + ['H'] * 11
exp = pd.DataFrame({'nuclei': nuc, 'sp2': [None] * len(nuc), 'exp_data': np.round(np.r_[np.linspace(15, 170, 14), np.linspace(0.9, 5.8, 11)], 2),
                    'exchange': [None] * len(nuc)})
spread = np.array([1.0, 1.6, 2.6, 5.0])
p_sca = np.clip(1 - np.abs(np.random.standard_t(5, size=(len(nuc), 4))) * spread / 9.0, 2e-9, 0.9999)
p_uns = np.clip(1 - np.abs(np.random.standard_t(4, size=(len(nuc), 4))) * spread / 8.0, 2e-9, 0.9999)
p_sca[20, 3] = 3.0e-9                                           # one very unlikely nucleus: its digits matter
isC = np.array([n == 'C' for n in nuc])
p_abs = pd.DataFrame(columns=isoms)
for n_, mask in (('C', isC), ('H', ~isC)):
    p_abs.loc['%s_uns' % n_] = np.prod(p_uns[mask], axis=0)
    p_abs.loc['%s_sca' % n_] = np.prod(p_sca[mask], axis=0)
    p_abs.loc['%s_full' % n_] = p_abs.loc['%s_uns' % n_] * p_abs.loc['%s_sca' % n_]
p_abs.loc['Uns'] = p_abs.loc['H_uns'] * p_abs.loc['C_uns']
p_abs.loc['Sca'] = p_abs.loc['H_sca'] * p_abs.loc['C_sca']
p_abs.loc['Full'] = p_abs.loc['C_full'] * p_abs.loc['H_full']
p_rel = pd.DataFrame(columns=isoms)
for kind, vec in p_abs.iterrows():
    p_rel.loc[kind] = vec / np.sum(vec)
results = pd.DataFrame(columns=isoms)
for row in ['H_sca', 'C_sca', 'Sca', 'H_uns', 'C_uns', 'Uns', 'H_full', 'C_full', 'Full']:
    results.loc[row] = p_rel.loc[row]
outputs = {'exp': exp, 'results': results, 'p_sca': p_sca, 'p_uns': p_uns}
path = os.path.join(HERE, 'dp4plus_app_results.xlsx')
with pd.ExcelWriter(path) as writer:
    for sheet, matrix in outputs.items():
        if 'exp' in sheet:
            continue
        try:
            matrix.to_excel(writer, sheet_name=sheet, index=True, float_format="%.4f", startrow=2)
        except Exception:
            matrix = pd.DataFrame(matrix, columns=isoms, index=outputs['exp'].index)
            matrix = pd.concat([outputs['exp'], matrix], axis=1)
            matrix.to_excel(writer, sheet_name=sheet, index=False, float_format="%.10f")
wb = load_workbook(path)
ws = wb['results']
ws.cell(1, 1).value = 'DP4+ QM'
ws.cell(14, 1).value = 'Cite this: J. Nat. Prod. 2023, 86, 10, 2360-2367 - '
ws.cell(18, 1).value = 'Theory Level selected : B3LYP/6-31G(d) (GAS)'
wb.save(path)
exact = (p_abs.loc['Full'] / p_abs.loc['Full'].sum()).astype(float).values
with open(os.path.join(HERE, 'dp4plus_app_results.expected.txt'), 'w') as f:
    f.write('\n'.join('%s %.15e' % (i, x) for i, x in zip(isoms, exact)) + '\n')

# ------------------------------------------------------------------ DP4+ Excel sheet (layout only)
wb = Workbook()
ws = wb.active; ws.title = 'Main'
for c, t in ((1, 'Functional'), (3, 'Solvent?'), (5, 'Basis Set'), (7, 'Type of Data')):
    ws.cell(1, c).value = t
for c, t in ((1, 'B3LYP'), (3, 'GAS'), (5, '6-31G(d)'), (7, 'Shielding Tensors')):
    ws.cell(2, c).value = t
labels = ['H (sDP4+)', 'C (sDP4+)', 'sDP4+', 'H (uDP4+)', 'C (uDP4+)', 'uDP4+', 'H (DP4+)', 'C (DP4+)', 'DP4+']
probs = [0.862, 0.1379999, 1.0e-7, 0.0]
for i, lab in enumerate(labels):
    ws.cell(4 + i, 3).value = lab
    for k in range(16):
        ws.cell(4 + i, 4 + k).value = (probs[k] if lab == 'DP4+' else round(1.0 / 4, 6)) if k < 4 else '-'
ws.cell(14, 1).value = 'Nuclei'; ws.cell(14, 2).value = 'sp2?'; ws.cell(14, 3).value = 'Experimental'
for k in range(16):
    ws.cell(14, 4 + k).value = 'Isomer %d' % (k + 1)
for r in range(18):
    ws.cell(15 + r, 1).value = 'C' if r < 10 else 'H'
    ws.cell(15 + r, 3).value = round(20 + 9.1 * r, 2) if r < 10 else round(0.8 + 0.55 * (r - 10), 2)
    for k in range(4):
        ws.cell(15 + r, 4 + k).value = round(100 + random.random(), 4)
wb.create_sheet('Calculo DP4')
wb.save(os.path.join(HERE, 'dp4plus_sheet.xlsx'))

# ------------------------------------------------------------------ PyDP4 output, both layouts
def table(nuc, rows):
    return '\n'.join('%-5s %8.2f %6.2f %6.2f %6.2f' % (('%s%d' % (nuc, i + 1),) + r) for i, r in enumerate(rows))
expC = [9.63, 13.65, 20.21, 25.19, 27.61, 47.67, 49.50, 61.75, 170.00]
expH = [0.92, 1.28, 1.45, 1.85, 2.74, 3.12, 4.10]
errs = {1: ([-0.07, 0.36, 2.18, -2.15, -2.02, 1.45, -1.82, -0.03, 0.28], [-0.07, 0.17, -0.04, 0.13, 0.01, 0.18, -0.22]),
        2: ([-0.64, 1.26, 0.41, 0.65, -3.52, 0.13, -0.27, 2.41, -0.48], [0.11, -0.25, 0.31, -0.08, 0.27, -0.33, 0.19]),
        3: ([3.10, -4.26, 2.75, -1.95, 5.02, -3.88, 2.14, -6.31, 1.07], [0.41, -0.52, 0.38, 0.66, -0.47, 0.29, -0.71])}
def rows(expv, ev):
    return [(x + e + 2.0, x + e, x, e) for x, e in zip(expv, ev)]
new = ['TEST_1', '', 'Solvent = chloroform', 'Force Field = mmff', '', 'DFT NMR Functional = mPW1PW91', 'DFT NMR Basis = 6-311g(d)', '', 'Number of isomers = 3', '']
for nuc_, key in (('C', 0), ('H', 1)):
    for iso in (1, 2, 3):
        new += ['Assigned %s shifts for isomer %d: ' % (nuc_, iso), 'label, calc, corrected, exp, error', table(nuc_, rows(expC if nuc_ == 'C' else expH, errs[iso][key])), '']
SIG = {'C': 2.269372270818724, 'H': 0.18731058105269952}
def ll(iso):
    return sum(math.log(math.erfc(abs(e) / SIG['C'] / math.sqrt(2))) for e in errs[iso][0]) + sum(math.log(math.erfc(abs(e) / SIG['H'] / math.sqrt(2))) for e in errs[iso][1])
L = [ll(i) for i in (1, 2, 3)]
P = [math.exp(x - max(L)) for x in L]; P = [100 * x / sum(P) for x in P]
new += ['Results of DP4 using Proton: ', 'Isomer 1: 50.0%', 'Isomer 2: 49.9%', 'Isomer 3:  0.1%', '', 'Results of DP4 using Carbon: ', 'Isomer 1: 60.0%', 'Isomer 2: 40.0%',
        'Isomer 3:  0.0%', '', 'Results of DP4: '] + ['Isomer %d: %4.1f%%' % (i + 1, p) for i, p in enumerate(P)] + ['']
open(os.path.join(HERE, 'pydp4_new.dp4'), 'w').write('\n'.join(new))
old = ['170(any),61.75(any)', '', '']
for iso in (1, 2, 3):
    old += ['Assigned shifts for isomer %d: ' % iso, '', 'Assigned C shifts: (label, calc, corrected, exp, error)', table('C', rows(expC, errs[iso][0])), 'Max C error:  %5.2f' % max(errs[iso][0], key=abs), '',
            'Assigned H shifts: (label, calc, corrected, exp, error)', table('H', rows(expH, errs[iso][1])), 'Max H error:  %5.2f' % max(errs[iso][1], key=abs), '']
old += ['Results of DP4 using all available data:'] + ['Isomer %d: %4.1f%%' % (i + 1, p) for i, p in enumerate(P)] + ['', 'Results of DP4 using carbon data only:', 'Isomer 1: 60.0%', 'Isomer 2: 40.0%', 'Isomer 3:  0.0%', '']
open(os.path.join(HERE, 'pydp4_old.dp4'), 'w').write('\n'.join(old))
with open(os.path.join(HERE, 'pydp4.expected.txt'), 'w') as f:
    f.write('\n'.join('%.15e' % x for x in L) + '\n')
print('fixtures written to', HERE)
