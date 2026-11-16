"""Loader and DP4+ evaluation for the STX150 validation set (110 cases, six levels)."""
import csv, collections, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dp4plus as D
import joc_data as J

LEVELS = ['MMFF_PBE', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9', 'B3L_PBE', 'B3L_wB9']


def load(out):
    blocks = collections.OrderedDict()      # name -> dict(set, table, rows=[...])
    for r in csv.DictReader(open(os.path.join(out, 'stx_val_blocks.csv'))):
        b = blocks.setdefault(r['block'], dict(set=int(r['set']), table=int(r['table']), has_exp=r['has_exp'] == '1', rows=[]))
        b['rows'].append(dict(row=int(r['row']), nuc=r['nucleus'], label=r['label'], exp=(float(r['exp']) if r['exp'] != '' else None),
                              sh={l: (float(r[l]) if r[l] != '' else None) for l in LEVELS}))
    cases = collections.OrderedDict()
    for r in csv.DictReader(open(os.path.join(out, 'stx_val_cases.csv'))):
        cases[r['case']] = dict(set=int(r['set']), cands=r['candidates'].split())
    par = {r['level']: {k: float(v) for k, v in r.items() if k != 'level'} for r in csv.DictReader(open(os.path.join(out, 'stx_dp4plus_parameters.csv')))}
    return blocks, cases, par


def pair_rows(case_block, cand_block):
    """Map each row of the case block to the candidate's row: same table -> same printed row; otherwise by atom label."""
    if case_block['table'] == cand_block['table']:
        return {i: i for i in range(len(case_block['rows']))}
    return {i: j for i, j in J.pair(case_block['rows'], cand_block['rows'])}


def case_arrays(blocks, case, cands, level, sp2_overrides=None, c_threshold=100.0):
    """Rows usable at this level: experimental shift present and a shielding printed for every candidate."""
    cb = blocks[case]; rows = cb['rows']
    maps = {c: pair_rows(cb, blocks[c]) for c in cands}
    flags = J.sp2_flags(rows, c_threshold, sp2_overrides)
    keep = []
    for i, r in enumerate(rows):
        if r['exp'] is None:
            continue
        ok = True
        for c in cands:
            j = maps[c].get(i)
            if j is None or blocks[c]['rows'][j]['sh'][level] is None:
                ok = False; break
        if ok:
            keep.append(i)
    nuc = np.array([rows[i]['nuc'] for i in keep]); dexp = np.array([rows[i]['exp'] for i in keep]); sp2 = np.array([flags[i] for i in keep], bool)
    sh = {c: np.array([blocks[c]['rows'][maps[c][i]]['sh'][level] for i in keep]) for c in cands}
    return keep, nuc, dexp, sp2, sh


def dp4plus_case(blocks, case, cands, level, par, truncate_nu=False, sp2_overrides=None, c_threshold=100.0):
    keep, nuc, dexp, sp2, sh = case_arrays(blocks, case, cands, level, sp2_overrides, c_threshold)
    if len(keep) < 3:
        return None, None, len(keep)
    terms = [D.log_terms(sh[c], dexp, nuc, sp2, par[level], truncate_nu=truncate_nu) for c in cands]
    return D.probabilities(terms, nuc), terms, len(keep)
