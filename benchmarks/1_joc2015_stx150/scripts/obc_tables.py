"""Org. Biomol. Chem. 2025 (Luciano, Amichetti, Sarotti, Zanardi), SI part I, section 2: experimental shifts and
Boltzmann-averaged shieldings of compounds 1-39 at Levels A (SMD/wB97XD/6-31+G**//MMFF), B (PCM/mPW1PW91/6-31+G**//
B3LYP/6-31G*), C (PCM/mPW1PW91/6-31+G**//MMFF) and, for some compounds, D (PCM/mPW1PW91/6-31+G**//xTB).

Version 2. Reads every layout of the section:
  standard      index, nucleus, exp, levels of isomer a, levels of isomer b
  continuation  levels of isomers c and d on the same rows (no experimental column)
  pair          exp of one compound, its levels, exp of its diastereomer, its levels (each is the other's alternative)
  side by side  two such blocks printed next to each other, each with its own index and nucleus columns
The nucleus column of the SI has misprints (proton rows printed as C and the reverse). A row is taken as 1H when
its experimental value is at most 14 ppm and the shieldings of the compound's own isomer lie between 15 and 35 ppm
(median), and as 13C otherwise; changes are listed. Cells that cannot be the shielding of their nucleus are listed too.
Writes obc_compounds.json: {compound number: {table_line, reference, levels, rows: [{label, nucleus, nucleus_printed,
exp}], isomers: {a: {level: [shieldings]}, b: ..., (c, d)}, pair_with (for pairs: the other compound)}}.
In a pair, isomer 'a' is always the compound itself and 'b' its diastereomer.
Usage: obc_tables.py <d4ob01959h1_suppl.txt> <out dir>"""
import sys, os, re, json, collections
src, out = sys.argv[1:3]
lines = open(src, encoding='utf-8', errors='replace').read().split('\n')
a = next(i for i, l in enumerate(lines) if l.strip().startswith('2. Experimental chemical shifts') and i > 60)
b = next(i for i, l in enumerate(lines) if l.strip().startswith('3. MM-DP4+ probabilities') and i > a)
NUM = re.compile(r'^-?\d+(?:\.\d+)?$')
def groups(line):
    """-> [(label, nucleus, [numbers])] for every 'label nucleus numbers...' group of the line"""
    tok = line.split(); out = []; i = 0
    while i + 2 < len(tok) + 1 and i + 1 < len(tok):
        if tok[i + 1] in ('C', 'H', 'N') and i + 2 < len(tok) and NUM.match(tok[i + 2]) and len(tok[i]) <= 10:
            j = i + 2; vals = []
            while j < len(tok) and NUM.match(tok[j]) and not (j + 2 < len(tok) and tok[j + 1] in ('C', 'H', 'N') and NUM.match(tok[j + 2])): vals.append(float(tok[j])); j += 1
            out.append((tok[i], tok[i + 1], vals)); i = j
        else: return []
    return out if i == len(tok) else []
frags = []; cur = None; header = True; lev = []; n_exp_head = 0
for i in range(a, b):
    l = lines[i]; g = groups(l)
    if g and all(len(v) >= 3 for _, _, v in g):
        if cur is None or header:
            cur = dict(line=i + 1, rows=[], levels=[x for _, x in sorted(lev)], n_exp_head=n_exp_head, reference=''); frags.append(cur); header = False; lev = []; n_exp_head = 0
        cur['rows'].append(g); continue
    if 'Exp. data' in l or 'Level' in l:
        header = True; n_exp_head += l.count('Exp. data')
        for mm in re.finditer(r'Level ([ABCD])', l): lev.append((mm.start(), mm.group(1)))
    elif cur is not None and not cur['reference'] and re.search(r'\b(19|20)\d\d\b', l) and re.search(r'[A-Z][A-Za-z]+\.', l): cur['reference'] = l.strip()
def fix_nucleus(printed, exp, sh):
    """sh: the cells of the compound's own isomer in this row"""
    med = sorted(sh)[len(sh) // 2]
    return 'H' if (exp <= 14 and 15 <= med <= 35) else 'C'
def bad_cell(nuc, exp, x):
    """a cell that cannot be the shielding of this nucleus: 1H outside 15-35 ppm or more than 3 ppm from the
    experimental shift; 13C more than 40 ppm from it (references 31.9 and 190 ppm)"""
    return (not 15 <= x <= 35 or abs(31.9 - x - exp) > 3) if nuc == 'H' else abs(190.0 - x - exp) > 40
comp = collections.OrderedDict(); number = 0; fixes = []; kinds = collections.Counter()
def add(number, frag, rows, isomers, levels, pair=None):
    comp[number] = dict(table_line=frag['line'], reference=frag['reference'], levels=levels, rows=rows, isomers=isomers, pair_with=pair)
last_first = None
for f in frags:
    side = len(f['rows'][0]) == 2
    ncol = len(f['rows'][0][0][2])
    def looks_like_exp(k):                               # shift + shielding is roughly constant for a nucleus
        v = [r[0][2] for r in f['rows'] if not side]
        ok = [1 for x in v if (150 <= x[k] + x[k + 1] <= 230) or (28 <= x[k] + x[k + 1] <= 36)]
        return len(ok) >= 0.8 * len(v)
    if side: nlev = len(f['rows'][0][0][2]) - 1; kind = 'side'
    elif ncol in (7, 9): nlev = (ncol - 1) // 2; kind = 'standard'
    elif ncol == 6: nlev = 3; kind = 'continuation'
    elif ncol == 10: nlev = 4; kind = 'pair'
    elif ncol == 8: kind = 'pair' if (looks_like_exp(0) and looks_like_exp(4)) else 'continuation'; nlev = 3 if kind == 'pair' else 4
    else: raise SystemExit('layout not understood at line %d: %d columns' % (f['line'], ncol))
    levels = ['A', 'B', 'C', 'D'][:nlev]
    if side:                                             # two blocks next to each other: exp + levels each
        kinds['side by side'] += 1
        assert all(len(r) == 2 and len(r[0][2]) == 1 + nlev and len(r[1][2]) == 1 + nlev for r in f['rows']), f['line']
        for k in (0, 1):
            number += 1; rows = []; own = {l: [] for l in levels}; other = {l: [] for l in levels}
            for r in f['rows']:
                lab, nuc, v = r[k]; lab2, nuc2, v2 = r[1 - k]
                nn = fix_nucleus(nuc, v[0], v[1:])
                if nn != nuc: fixes.append((number, lab, nuc, nn, v[0]))
                rows.append(dict(label=lab, nucleus=nn, nucleus_printed=nuc, exp=v[0]))
                for q, l in enumerate(levels): own[l].append(v[1 + q]); other[l].append(v2[1 + q])
            add(number, f, rows, {'a': own, 'b': other}, levels, pair=number + (1 if k == 0 else -1))
        last_first = None; continue
    flat = [r[0] for r in f['rows']]
    assert len({len(v) for _, _, v in flat}) == 1, ('ragged table', f['line'])
    if kind == 'continuation':                           # isomers c and d
        kinds['continuation'] += 1
        assert last_first is not None and len(flat) == len(comp[last_first]['rows']), ('continuation does not fit', f['line'])
        for q, iso in enumerate('cd'):
            comp[last_first]['isomers'][iso] = {l: [v[q * nlev + p] for _, _, v in flat] for p, l in enumerate(levels)}
        continue
    if kind == 'standard':
        kinds['standard'] += 1; number += 1; rows = []
        for lab, nuc, v in flat:
            nn = fix_nucleus(nuc, v[0], v[1:1 + nlev])
            if nn != nuc: fixes.append((number, lab, nuc, nn, v[0]))
            rows.append(dict(label=lab, nucleus=nn, nucleus_printed=nuc, exp=v[0]))
        iso = {s: {l: [v[1 + q * nlev + p] for _, _, v in flat] for p, l in enumerate(levels)} for q, s in enumerate('ab')}
        add(number, f, rows, iso, levels); last_first = number; continue
    if kind == 'pair':                                   # exp, levels, exp, levels
        kinds['pair'] += 1
        for k in (0, 1):
            number += 1; rows = []; base = k * (1 + nlev); ob = (1 - k) * (1 + nlev)
            for lab, nuc, v in flat:
                nn = fix_nucleus(nuc, v[base], v[base + 1:base + 1 + nlev])
                if nn != nuc: fixes.append((number, lab, nuc, nn, v[base]))
                rows.append(dict(label=lab, nucleus=nn, nucleus_printed=nuc, exp=v[base]))
            own = {l: [v[base + 1 + p] for _, _, v in flat] for p, l in enumerate(levels)}; other = {l: [v[ob + 1 + p] for _, _, v in flat] for p, l in enumerate(levels)}
            add(number, f, rows, {'a': own, 'b': other}, levels, pair=number + (1 if k == 0 else -1))
        last_first = None; continue
cells = []
for n, c in comp.items():
    c['bad_cells'] = []
    for iso, d in c['isomers'].items():
        for l, col in d.items():
            for i, x in enumerate(col):
                if bad_cell(c['rows'][i]['nucleus'], c['rows'][i]['exp'], x): c['bad_cells'].append([iso, l, i]); cells.append((n, iso, l, c['rows'][i]['label'], c['rows'][i]['nucleus'], c['rows'][i]['exp'], x))
json.dump(comp, open(os.path.join(out, 'obc_compounds.json'), 'w'))
print('table fragments: %d %s -> compounds: %d' % (len(frags), dict(kinds), len(comp)))
print('isomers per compound:', dict(collections.Counter(len(c['isomers']) for c in comp.values())), '| with Level D:', [n for n, c in comp.items() if 'D' in c['levels']])
print('pairs (each compound is the alternative of the other):', [(n, c['pair_with']) for n, c in comp.items() if c['pair_with'] and n < c['pair_with']])
print('rows: %d (13C %d, 1H %d) | nucleus column corrected in %d rows of compounds %s' % (sum(len(c['rows']) for c in comp.values()), sum(r['nucleus'] == 'C' for c in comp.values() for r in c['rows']), sum(r['nucleus'] == 'H' for c in comp.values() for r in c['rows']), len(fixes), sorted({f[0] for f in fixes})))
for n, lab, p, q, e in fixes: print('   compound %2d row %-6s printed %s, taken as %s (experimental %s)' % (n, lab, p, q, e))
ncol = sum(len(c['isomers']) * len(c['levels']) for c in comp.values())
badcols = sorted({(n, iso, l) for n, iso, l, *_ in cells})
print('shielding columns (compound, isomer, level): %d | with cells that cannot be the shielding of their nucleus: %d (%d cells)' % (ncol, len(badcols), len(cells)))
bycomp = collections.OrderedDict()
for n, iso, l in badcols: bycomp.setdefault(n, []).append(iso + '/' + l)
for n, v in bycomp.items(): print('   compound %2d: %s | e.g. %s' % (n, ' '.join(v), '; '.join('%s%s exp %s cell %s (%s/%s)' % (nu, lab, e, x, iso, l) for nn, iso, l, lab, nu, e, x in cells if nn == n)[:110]))
