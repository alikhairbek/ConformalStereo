"""Shielding columns that are printed twice in the STX150 SI (validation set): (1) a candidate whose column at one
level equals its column at another level, (2) two candidates of one set with the same column at one level.
Usage: stx_copied_columns.py <stx out dir>"""
import sys, os, csv, collections
stx_out = sys.argv[1]
LEVELS = ['MMFF_PBE', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9', 'B3L_PBE', 'B3L_wB9']
blocks = collections.OrderedDict()
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_val_blocks.csv'))): blocks.setdefault(r['block'], []).append(r)
cases = list(csv.DictReader(open(os.path.join(stx_out, 'stx_val_cases.csv'))))
set_of = {}
for c in cases:
    for b in c['candidates'].split(): set_of.setdefault(b, c['candidates'])
vec = lambda rows, lev: [r[lev] for r in rows]
print('candidate blocks: %d' % len(blocks))
print('columns that are exact copies of another column of the same candidate (all nuclei, as printed):')
n = 0
for name, rows in blocks.items():
    for i, l1 in enumerate(LEVELS):
        for l2 in LEVELS[i + 1:]:
            a, b = vec(rows, l1), vec(rows, l2)
            same = sum(1 for x, y in zip(a, b) if x == y and x != '')
            if same >= 0.8 * len(a) and len(a) > 5:
                n += 1; print('  %-10s %s = %s for %d of %d nuclei | candidate set: %s | SI page %s' % (name, l1, l2, same, len(a), set_of.get(name, '?'), rows[0]['si_pdf_page']))
print('total', n)
print('same column printed for two candidates of one set:')
seen = set()
for c in cases:
    cand = c['candidates'].split()
    for i, a in enumerate(cand):
        for b in cand[i + 1:]:
            if (a, b) in seen: continue
            seen.add((a, b))
            for l in LEVELS:
                va, vb = vec(blocks[a], l), vec(blocks[b], l)
                if len(va) == len(vb) and len(va) > 5 and sum(1 for x, y in zip(va, vb) if x == y and x != '') >= 0.8 * len(va): print('  %s / %s identical at %s' % (a, b, l))
