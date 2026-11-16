"""MM-DP4+ paper (J. Nat. Prod. 2023) SI, section 10: probability of the correct isomer for every validation case,
six formulations at three MMFF-geometry levels, plus DP4 and DP4+ (DFT geometries). Writes jnp_probabilities.csv.
Usage: jnp_probs.py <np3c00566_si_001.txt (pdftotext -layout)> <out_dir>"""
import sys, os, re, csv, collections
src, out = sys.argv[1:3]
lines = open(src, encoding='utf-8', errors='replace').read().split('\n')
start = max(i for i, l in enumerate(lines) if 'ScTMS+ScMSTD+Uns' in l and i < 28000) if False else next(i for i, l in enumerate(lines) if 'ScTMS+ScMSTD+Uns' in l and i > 20000)
LEV = [('B3LYP/6-31G**//MMFF', 'MMFF_B3L'), ('PCM/mPW1PW91/6-31+G**//MMFF', 'MMFF_mPW'), ('SMD/ωB97XD/6-31+G**//MMFF', 'MMFF_wB9')]
FORMS = ['ScTMS', 'ScTMS+Uns', 'ScMSTD', 'ScMSTD+Uns', 'ScTMS+ScMSTD', 'ScTMS+ScMSTD+Uns']
ROW = re.compile(r'^\s*(\d{2,3}[a-z])\s+((?:\d+(?:\.\d+)?%\s+){5}\d+(?:\.\d+)?%)(.*)$')
OTHER = re.compile(r'(\d{2,3}[a-z])\s+(\d+(?:\.\d+)?)%')
data = collections.OrderedDict(); other = collections.OrderedDict(); level = None; om = None
def tail(rest):
    global om
    if re.search(r'\bDP4\+', rest): om = 'DP4+'
    elif re.search(r'\bDP4\b', rest): om = 'DP4'
    for c, v in OTHER.findall(rest):
        if om: other.setdefault(c, {})[om] = float(v)
for l in lines[start - 2:]:
    if 'ScTMS' in l and 'Uns' in l: level = None; om = None; continue
    hit = next((code for name, code in LEV if name in l), None)
    if hit:
        level = hit; tail(l.split('//MMFF', 1)[1]); continue
    m = ROW.match(l)
    if m and level:
        vals = [float(x.rstrip('%')) for x in m.group(2).split()]
        data.setdefault(m.group(1), {})[level] = vals; tail(m.group(3))
    elif l.strip():
        tail(l)
rows = []
for c, d in data.items():
    r = dict(case=c)
    for lv in ('MMFF_B3L', 'MMFF_mPW', 'MMFF_wB9'):
        for f, v in zip(FORMS, d.get(lv, [None] * 6)): r['%s|%s' % (lv, f)] = v
    r['DP4'] = other.get(c, {}).get('DP4'); r['DP4+'] = other.get(c, {}).get('DP4+')
    rows.append(r)
w = csv.DictWriter(open(os.path.join(out, 'jnp_probabilities.csv'), 'w', newline=''), fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
nums = sorted(set(int(re.match(r'\d+', c).group()) for c in data))
print('cases: %d | compounds %d to %d | cases with all three levels: %d | with DP4: %d | with DP4+: %d' % (len(data), nums[0], nums[-1], sum(len(d) == 3 for d in data.values()),
      sum(r['DP4'] is not None for r in rows), sum(r['DP4+'] is not None for r in rows)))
print('first and last:', rows[0]['case'], rows[-1]['case'], '| example 121a:', {k: v for k, v in next(r for r in rows if r['case'] == '121a').items() if 'Uns' in k and 'ScMSTD+' not in k or k in ('DP4', 'DP4+')})
