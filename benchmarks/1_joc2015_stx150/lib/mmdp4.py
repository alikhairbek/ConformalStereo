"""MM-DP4+ and DP4+ as DP4+App computes them (dp4_module.py of version 2.2.3): unscaled shift = TMS - shielding;
scaled shift from the straight line of unscaled against experimental shifts, 13C and 1H apart; each error goes through
a Student t with the parameters of its set (sp2 or sp3 for the unscaled errors, one set for the scaled errors); the
degrees of freedom are used as they are (no rounding, unlike the 2015 Excel sheet). The formula is the one of
lib/dp4plus.py; this module adds the parameter sets shipped with the program and the level names of the 2023 paper."""
import csv, os, sys
import numpy as np
FUNC = {'B3L': 'B3LYP', 'M06': 'M062x', 'MPW': 'mPW1PW91', 'WB9': 'wB97XD'}
BASIS = {'A': '6-31G(d,p)', 'B': '6-31+G(d,p)', 'C': '6-311+G(d,p)'}
SOLV = {'GAS': '', 'PCM': '.PCM', 'SMD': '.SMD'}
def app_level(code):
    """B3L_A_GAS -> B3LYP.6-31G(d,p) ; MPW_B_PCM -> mPW1PW91.6-31+G(d,p).PCM"""
    f, b, s = code.split('_'); return FUNC[f] + '.' + BASIS[b] + SOLV[s]
def load_parameters(path, family):
    return {r['level']: {k: float(v) for k, v in r.items() if k not in ('family', 'level')} for r in csv.DictReader(open(path)) if r['family'] == family}
def probabilities(D, shieldings, exp, nuc, sp2, par):
    """shieldings: list of arrays (one per candidate), rows as exp/nuc/sp2. Returns the dict of dp4plus.probabilities."""
    nuc = np.asarray(nuc); exp = np.asarray(exp, float); sp2 = np.asarray(sp2, bool)
    terms = [D.log_terms(np.asarray(s, float), exp, nuc, sp2, par, truncate_nu=False) for s in shieldings]
    return D.probabilities(terms, nuc)

import re
def proton_pairs(rows):
    """Pairs of protons on the same carbon, from the labels: two 1H rows with the same atom number (8a/8b, 6n/6x, 3 and
    3, 2' and 2''). rows: dicts with nuc and label. Returns a list of index pairs."""
    groups = {}
    for i, r in enumerate(rows):
        if r['nuc'] != 'H': continue
        m = re.match(r"^(\d+)", r['label'].strip())
        if m: groups.setdefault(m.group(1) + ("'" if "'" in r['label'] and not r['label'].strip()[len(m.group(1)):].strip("'abcdefxnsαβ") else ''), []).append(i)
    return [tuple(v) for v in groups.values() if len(v) == 2]
def exchange(shield, exp, pairs):
    """DP4+App, correlation_module.diasterotopic: inside a pair marked as exchangeable the larger experimental shift
    goes with the smaller shielding, for each candidate on its own."""
    s = list(shield)
    for i, j in pairs:
        lo, hi = min(s[i], s[j]), max(s[i], s[j])
        if exp[i] >= exp[j]: s[i], s[j] = lo, hi
        else: s[i], s[j] = hi, lo
    return s
