"""Internal checks of the shielding tables printed in the STX150 SI (validation set, six levels):
(1) is a column a copy of another column? (2) how far apart are the columns of one compound?
Usage: stx_column_checks.py <out dir of the SI extraction> [case ...]"""
import sys, os, csv, math, collections
stx_out = sys.argv[1]; focus = sys.argv[2:]
LEVELS = ['MMFF_PBE', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9', 'B3L_PBE', 'B3L_wB9']
blocks = collections.OrderedDict()
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_val_blocks.csv'))):
    blocks.setdefault(r['block'], []).append(r)
def col(rows, lev, nuc):
    return [float(r[lev]) if r[lev] != '' else None for r in rows if r['nucleus'] == nuc]
def rms(a, b):
    d = [x - y for x, y in zip(a, b) if x is not None and y is not None]
    return math.sqrt(sum(v * v for v in d) / len(d)) if d else None
copies = []; stat = collections.defaultdict(list)
for name, rows in blocks.items():
    for i, l1 in enumerate(LEVELS):
        for l2 in LEVELS[i + 1:]:
            rc = rms(col(rows, l1, 'C'), col(rows, l2, 'C')); rh = rms(col(rows, l1, 'H'), col(rows, l2, 'H'))
            if rc is None: continue
            stat[(l1, l2)].append(rc)
            if rc < 0.02 and (rh is None or rh < 0.005): copies.append((name, l1, l2, rc))
print('compounds: %d | pairs of columns that are identical for a compound: %d %s' % (len(blocks), len(copies), copies[:12]))
print('typical distance between two columns of one compound (13C, RMS in ppm; median over compounds):')
for k, v in stat.items():
    v.sort(); print('   %-9s vs %-9s median %5.2f | smallest %5.2f | largest %5.2f' % (k[0], k[1], v[len(v) // 2], v[0], v[-1]))
for name in focus:
    rows = blocks.get(name)
    if not rows: print(name, 'not found'); continue
    print('%s: 13C RMS between columns: %s' % (name, ' '.join('%s/%s=%.2f' % (a[:4] + a[5:], b[:4] + b[5:], rms(col(rows, a, 'C'), col(rows, b, 'C'))) for i, a in enumerate(LEVELS) for b in LEVELS[i + 1:] if (a, b) in (('MMFF_PBE', 'xTB_PBE'), ('MMFF_PBE', 'B3L_PBE'), ('xTB_PBE', 'B3L_PBE'), ('xTB_PBE', 'xTB_mPW'), ('B3L_PBE', 'B3L_wB9')))))
