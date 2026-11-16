"""2023 SI, section 8: MM-DP4+ probability (%) of the correct isomer for 34 validation molecules at the 36 levels of
theory, and the score table that follows it. Level labels are printed vertically; their order is taken from their
horizontal position. Writes jnp_section8.csv (molecule, level, probability, score).
Usage: jnp_section8.py <np3c00566_si_001.txt> <out csv>"""
import sys, re, csv, collections
lines = open(sys.argv[1], encoding='utf-8', errors='replace').read().split('\n')
a = next(i for i, l in enumerate(lines) if '8. MM-DP4+ probabilities (%) and scores' in l and i > 1000)
b = next(i for i, l in enumerate(lines) if l.strip('\x0c ').startswith('9. Experimental chemical shifts') and i > a)
LEV = re.compile(r'\b((?:B3L|M06|MPW|WB9)_[ABC]_(?:GAS|PCM|SMD))\b')
ROW = re.compile(r'^\s*(\d{2,3}[a-z])\s+(.*\S)\s*$')
blocks = []; head = []; cur = None
for i in range(a, b):
    l = lines[i]; m = ROW.match(l)
    if m and len(m.group(2).split()) >= 30:
        if cur is None or head:
            cur = dict(levels=[x for _, x in sorted(head)], rows=collections.OrderedDict(), pct='%' in m.group(2)); blocks.append(cur); head = []
        cur['rows'][m.group(1)] = [t.rstrip('%') for t in m.group(2).split()]
        continue
    for mm in LEV.finditer(l): head.append((mm.start(), mm.group(1)))
print('blocks:', [(len(b_['rows']), len(b_['levels']), 'percent' if b_['pct'] else 'score') for b_ in blocks])
for b_ in blocks: print('   levels in column order:', ' '.join(b_['levels']))
prob = [b_ for b_ in blocks if b_['pct']]; score = [b_ for b_ in blocks if not b_['pct']]
out = []
for b_ in prob:
    assert len(b_['levels']) == 36 and all(len(v) == 36 for v in b_['rows'].values())
    for mol, vals in b_['rows'].items():
        for lev, v in zip(b_['levels'], vals): out.append(dict(molecule=mol, level=lev, probability=float(v), score=''))
sc = {}
for b_ in score:
    if len(b_['levels']) != 36: b_['levels'] = prob[0]['levels']
    for mol, vals in b_['rows'].items():
        for lev, v in zip(b_['levels'], vals): sc[(mol, lev)] = v
for o in out: o['score'] = sc.get((o['molecule'], o['level']), '')
w = csv.DictWriter(open(sys.argv[2], 'w', newline=''), fieldnames=['molecule', 'level', 'probability', 'score']); w.writeheader(); w.writerows(out)
mols = list(collections.OrderedDict.fromkeys(o['molecule'] for o in out))
print('molecules: %d | values: %d | molecules: %s' % (len(mols), len(out), ' '.join(mols)))
tab = collections.Counter((o['score'], 'p>=95' if o['probability'] >= 95 else ('50<=p<95' if o['probability'] >= 50 else 'p<50')) for o in out)
print('score against probability:', dict(sorted(tab.items())))
