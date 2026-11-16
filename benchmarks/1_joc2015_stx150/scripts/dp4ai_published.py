"""Published per-isomer probabilities of the DP4-AI paper (spreadsheet of its supplementary information):
read the three sheets, repair the molecule column, and write one long table.
Usage: dp4ai_published.py <d0sc00442a1_suppl.xlsx> <out_dir>      (needs openpyxl)"""
import sys, os, csv, collections
import openpyxl
src, out = sys.argv[1:3]
wb = openpyxl.load_workbook(src, data_only=True, read_only=True)
LEVEL = {'MM_3g': 'MM', 'Opt_crossvalidation': 'Opt', 'Opt_E_crossvalidation': 'Opt_E'}
rows = []
def num(v):
    if v is None: return None
    try:
        f = float(v)
        return None if f != f else f
    except (TypeError, ValueError):
        return None
for ws in wb.worksheets:
    if ws.title not in LEVEL: continue
    data = list(ws.iter_rows(values_only=True))
    h0 = next(i for i, r in enumerate(data) if r and any(str(c).strip().lower() == 'stereochemistry' for c in r if c is not None))
    head = [str(c).strip().lower() if c is not None else '' for c in data[h0]]
    ix = {}
    for i, h in enumerate(head):
        if h and h not in ix: ix[h] = i
    if 'molecule' not in ix: ix['molecule'] = ix['stereochemistry'] - 1          # unnamed first column
    mol = None; k = 0
    for r in data[h0 + 1:]:
        if r is None or all(c is None for c in r[:9]): continue
        if r[ix['molecule']] not in (None, ''): mol = str(r[ix['molecule']]).strip(); 
        st = r[ix['stereochemistry']]
        if st in (None, ''): continue
        k += 1
        rows.append(dict(level=LEVEL[ws.title], molecule=mol, stereo=str(st).strip(), row=k,
                         ai_all=num(r[ix['overall']]), ai_C=num(r[ix['carbon']]), ai_H=num(r[ix['proton']]),
                         marked=1 if str(r[ix['correct']] or '').strip().lower() == 'x' else 0,
                         pair_all=num(r[ix['pair_overall']]), pair_C=num(r[ix['pair_carbon']]), pair_H=num(r[ix['pair_proton']])))
w = csv.DictWriter(open(os.path.join(out, 'dp4ai_published_raw.csv'), 'w', newline=''), fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for lvl in ('MM', 'Opt', 'Opt_E'):
    rs = [r for r in rows if r['level'] == lvl]
    mols = collections.OrderedDict()
    for r in rs: mols.setdefault(r['molecule'], []).append(r)
    marks = collections.Counter(sum(x['marked'] for x in v) for v in mols.values())
    pos = collections.Counter(next((i for i, x in enumerate(v) if x['marked']), None) for v in mols.values())
    sums = [round(sum(x['ai_all'] or 0 for x in v)) for v in mols.values()]
    print('%-5s rows %3d | molecules %2d | marks per molecule %s | position of the mark within the molecule block %s | sums of the probabilities != 100 in %d molecules' % (
        lvl, len(rs), len(mols), dict(marks), dict(sorted(pos.items(), key=lambda kv: str(kv[0]))), sum(abs(s - 100) > 1 for s in sums)))
names = collections.OrderedDict((r['molecule'], 1) for r in rows if r['level'] == 'MM')
print('molecules:', ' '.join(names))
