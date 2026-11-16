"""MM-DP4+ paper (J. Nat. Prod. 2023) SI, sections 7 and 9: experimental shifts and Boltzmann-averaged shieldings
(MMFF geometries) of the validation molecules. Reads every table and writes jnp_tables.json:
[{line, labels (compound codes if printed), levels (in column order), rows: [[atom label, nucleus, v1, v2, ...]]}].
Atom labels are kept as printed (1, 3a, 6x, 1', ...). A table that runs over a page break is put together again:
a fragment is a continuation when it has no carbon rows, or when the table before it has no proton rows yet.
Usage: jnp_tables.py <np3c00566_si_001.txt> <out_dir> [first_line last_line]"""
import sys, os, re, json, collections
src, out = sys.argv[1:3]
lines = open(src, encoding='utf-8', errors='replace').read().split('\n')
a = int(sys.argv[3]) if len(sys.argv) > 3 else next(i for i, l in enumerate(lines) if re.match(r'^\s*S156\s*$', l))
b = int(sys.argv[4]) if len(sys.argv) > 4 else next(i for i, l in enumerate(lines) if 'ScTMS+ScMSTD+Uns' in l and i > a)
NUM = r'-?\d+(?:\.\d+)?'
ROW = re.compile(r'^\s*(\S{1,10})\s+([CHN])\s+((?:' + NUM + r'\s+){4,}' + NUM + r')\s*$')
LEV = re.compile(r'\b((?:B3L|M06|MPW|WB9)_[ABC]_(?:GAS|PCM|SMD))\b')
CODE = re.compile(r'\b(\d{2,3}[a-z])\b(?:\s*\(([^)]*)\))?')
frags = []; cur = None; head_levels = []; head_codes = []; last_row_line = -99
for i in range(a, b):
    l = lines[i]
    if 'Downloaded from' in l: l = l[:l.index('Downloaded from')]
    m = ROW.match(l)
    if m:
        vals = [float(x) for x in m.group(3).split()]
        if cur is None or i - last_row_line > 12 or head_levels:
            cur = dict(line=i + 1, labels=list(head_codes), levels=[n for _, n in sorted(head_levels)], rows=[])
            frags.append(cur); head_levels = []; head_codes = []
        cur['rows'].append([m.group(1), m.group(2)] + vals); last_row_line = i
        continue
    for mm in LEV.finditer(l): head_levels.append((mm.start(), mm.group(1)))
    cs = [(mm.start(), mm.group(1) + ((' ' + mm.group(2)) if mm.group(2) else '')) for mm in CODE.finditer(l)]
    if len(cs) >= 1 and not LEV.search(l) and len(l.strip()) < 90 and not re.search(r'[A-Za-z]{4,}', re.sub(r'\([^)]*\)', '', l)):
        head_codes = [c for _, c in sorted(cs)]
tables = []; joined = 0
for f in frags:
    prev = tables[-1] if tables else None
    width = lambda t: collections.Counter(len(r) for r in t['rows']).most_common(1)[0][0]
    cont = prev is not None and width(prev) == width(f) and (not any(r[1] == 'C' for r in f['rows']) or not any(r[1] == 'H' for r in prev['rows']))
    if cont and (not f['levels'] or f['levels'] == prev['levels']):
        prev['rows'].extend(f['rows']); joined += 1
    else: tables.append(f)
json.dump(tables, open(os.path.join(out, 'jnp_tables.json'), 'w'))
w = collections.Counter(len(t['rows'][0]) - 2 for t in tables)
odd = sum(1 for t in tables if len({len(r) for r in t['rows']}) > 1)
print('fragments: %d | joined over a page break: %d | tables: %d | lines %d to %d' % (len(frags), joined, len(tables), a, b))
print('numeric columns per row -> tables: %s | tables whose rows differ in length: %d' % (dict(sorted(w.items())), odd))
print('rows: %d (C %d, H %d, N %d) | tables without proton rows: %d | without carbon rows: %d' % (sum(len(t['rows']) for t in tables),
      sum(r[1] == 'C' for t in tables for r in t['rows']), sum(r[1] == 'H' for t in tables for r in t['rows']), sum(r[1] == 'N' for t in tables for r in t['rows']),
      sum(1 for t in tables if not any(r[1] == 'H' for r in t['rows'])), sum(1 for t in tables if not any(r[1] == 'C' for r in t['rows']))))
