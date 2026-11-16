"""Pairs of molecules that share a constitution (two diastereomers, each with its own spectra): does the analysis of A's
spectra prefer the structure reported for B, and the other way round? A double preference would point to exchanged labels.
Usage: swap_test.py <res/core> <out_dir>      (needs RDKit)"""
import sys, os, re, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dp4lib as L
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
core, out = sys.argv[1:3]
S = L.Small(core, 'dp4ai'); ps = S.paths()
def keys(case):
    d = {}
    for p in ps:
        m = re.match(r'^%s\.zip!/[^/]+/(MM|Opt)/[^/]+_(\d+)\.sdf$' % re.escape(case), p)
        if m and int(m.group(2)) not in d:
            mol = Chem.MolFromMolBlock(S.read(p), removeHs=False); Chem.AssignStereochemistryFrom3D(mol)
            d[int(m.group(2))] = re.sub(r'/m[01]', '', Chem.MolToInchi(mol))
    return d
P = collections.defaultdict(dict)
for r in csv.DictReader(open(os.path.join(out, 'dp4ai_probabilities.csv'))):
    if r['data'] == 'all': P[(r['case'], r['protocol'])][int(r['isomer'])] = float(r['probability'])
print('%-6s %-6s | isomer of A that is the structure reported for B | P(own structure) / P(partner structure) in the analysis of A, per protocol' % ('A', 'B'))
for a, b in (('JB2', 'JB3'), ('JB13A', 'JB13B'), ('NL1A', 'NL1B'), ('NL2A', 'NL2B'), ('NP3A', 'NP3B'), ('TS3A', 'TS3B')):
    for x, y in ((a, b), (b, a)):
        kx, ky = keys(x), keys(y)
        partner = [n for n, k in kx.items() if k == ky[1]]
        cells = []
        for prot in ('MM', 'DFT_EB3LYP', 'DFT_EmPW', 'DFT_EM062X'):
            pr = P.get((x, prot))
            if pr and partner: cells.append('%s %5.1f/%5.1f' % (prot.replace('DFT_', ''), 100 * pr.get(1, 0), 100 * pr.get(partner[0], 0)))
        print('%-6s %-6s | isomer %-2s | %s' % (x, y, partner[0] if partner else '?', ' | '.join(cells)))
