"""Parameter sets shipped with DP4+App 2.2.3 (PyPI wheel, folder parameters_bank): TMS standards and the six Student-t
parameter sets (n, m, s for Csp3, Csp2, Csca, Hsp3, Hsp2, Hsca) of every level, for the three geometry families
(QM: B3LYP/6-31G* geometries, DP4+; MM: MMFF geometries, MM-DP4+; xTB: GFN-xTB geometries).
Writes one CSV in the column layout used by dp4plus.py. The scaled sets have m = 0 in every sheet (checked).
Usage: app_parameters.py <parameters_bank dir> <out csv>     (needs openpyxl)"""
import sys, os, csv, openpyxl
bank, out_csv = sys.argv[1:3]
rows = []; nonzero = 0
for fam, fname, std_sheet in (('QM', 'data_base_QM.xlsx', 'standard'), ('MM', 'data_base_MM.xlsx', 'Standard'), ('xTB', 'data_base_xTB.xlsx', 'CHCl3')):
    wb = openpyxl.load_workbook(os.path.join(bank, fname), read_only=True, data_only=True)
    std = {str(r[0]).strip(): (float(r[1]), float(r[2])) for r in wb[std_sheet].iter_rows(values_only=True) if r and r[0] and isinstance(r[1], (int, float))}
    for name in wb.sheetnames:
        if name not in std: continue
        tab = {str(r[0]).strip(): r[1:4] for r in wb[name].iter_rows(values_only=True) if r and r[0] and isinstance(r[1], (int, float))}
        if not all(k in tab for k in ('Csp3', 'Csp2', 'Csca', 'Hsp3', 'Hsp2', 'Hsca')): continue
        nonzero += (tab['Csca'][1] != 0) + (tab['Hsca'][1] != 0)
        rec = dict(family=fam, level=name, TMS_C=std[name][0], TMS_H=std[name][1])
        for nuc in 'CH':
            for cls in ('sp2', 'sp3'):
                n, m, s = tab[nuc + cls]; rec.update({'%s_%s_mu' % (nuc, cls): m, '%s_%s_sigma' % (nuc, cls): s, '%s_%s_nu' % (nuc, cls): n})
            n, m, s = tab[nuc + 'sca']; rec.update({'%s_sc_sigma' % nuc: s, '%s_sc_nu' % nuc: n})
        rows.append(rec)
keys = ['family', 'level', 'TMS_C', 'TMS_H'] + ['%s_%s_%s' % (n, c, p) for n in 'CH' for c in ('sp2', 'sp3') for p in ('mu', 'sigma', 'nu')][:6] + ['C_sc_sigma', 'C_sc_nu'] + ['H_%s_%s' % (c, p) for c in ('sp2', 'sp3') for p in ('mu', 'sigma', 'nu')] + ['H_sc_sigma', 'H_sc_nu']
w = csv.DictWriter(open(out_csv, 'w', newline=''), fieldnames=keys); w.writeheader(); w.writerows(rows)
import collections
print('levels written:', dict(collections.Counter(r['family'] for r in rows)), '| scaled sets with m different from 0:', nonzero)
for r in rows:
    if r['level'] in ('mPW1PW91.6-31+G(d,p).PCM', 'wB97XD.6-31+G(d,p).SMD', 'PBE1PBE.6-31G(d,p).SMD'): print('  %-3s %-26s TMS %.4f %.4f | Csp2 %.3f %.3f %.3f | Csp3 %.3f %.3f %.3f | Csca %.3f %.3f' % (r['family'], r['level'], r['TMS_C'], r['TMS_H'], r['C_sp2_mu'], r['C_sp2_sigma'], r['C_sp2_nu'], r['C_sp3_mu'], r['C_sp3_sigma'], r['C_sp3_nu'], r['C_sc_sigma'], r['C_sc_nu']))
