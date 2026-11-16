"""Case 89a at MMFF_PBE: the STX150 table says 0%, the shieldings printed in the STX150 SI give 99%.
Are the printed MMFF_PBE columns of 89a and 89b the isomers of the 2023 paper, in the same order?
Usage: stx_jnp_89.py <jnp_tables.json> <stx out dir>"""
import sys, os, json, csv, math, collections
T = json.load(open(sys.argv[1])); stx_out = sys.argv[2]
blocks = collections.OrderedDict()
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_val_blocks.csv'))): blocks.setdefault(r['block'], []).append(r)
def corr(x, y):
    n = len(x); mx, my = sum(x) / n, sum(y) / n; sxx = sum((a - mx) ** 2 for a in x); syy = sum((b - my) ** 2 for b in y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)
def line_fit(x, y):
    n = len(x); mx, my = sum(x) / n, sum(y) / n; sxx = sum((a - mx) ** 2 for a in x); sxy = sum((a - mx) * (b - my) for a, b in zip(x, y)); s = sxy / sxx; i = my - s * mx
    res = [b - (s * a + i) for a, b in zip(x, y)]
    return math.sqrt(sum(v * v for v in res) / n), max(abs(v) for v in res)
def key(t): return sorted(round(r[2], 1) for r in t['rows'] if r[1] == 'C')
def table_for(block, level):
    ks = sorted(round(float(r['exp']), 1) for r in blocks[block] if r['nucleus'] == 'C')
    return next(t for t in T if level in t['levels'] and key(t) == ks)
SL = ['MMFF_PBE', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9', 'B3L_PBE', 'B3L_wB9']
A = table_for('89a', 'MPW_B_PCM'); B = table_for('89b', 'MPW_B_PCM')
col = lambda t, lev: 2 + 1 + t['levels'].index(lev)
print('rows: STX 89a %d, 89b %d | 2023 SI 89a %d, 89b %d | atom labels of the two 2023 tables equal: %s' % (len(blocks['89a']), len(blocks['89b']), len(A['rows']), len(B['rows']), [r[:2] for r in A['rows']] == [r[:2] for r in B['rows']]))
print('[STX 89a - STX 89b] against [2023 89a - 2023 89b] (MMFF geometries), r(13C)/r(1H):')
for jl in ('MPW_B_PCM', 'B3L_B_PCM', 'WB9_B_PCM'):
    out = []
    for sl in SL:
        rr = []
        for nuc in 'CH':
            idx = [i for i, r in enumerate(A['rows']) if r[1] == nuc]
            rr.append(corr([float(blocks['89a'][i][sl]) - float(blocks['89b'][i][sl]) for i in idx], [A['rows'][i][col(A, jl)] - B['rows'][i][col(B, jl)] for i in idx]))
        out.append('%s %.2f/%.2f' % (sl, rr[0], rr[1]))
    print('  2023 level %s: %s' % (jl, ' | '.join(out)))
print('fit of each STX block (MMFF_PBE) against each 2023 isomer (MPW_B_PCM), RMS of the residuals, 13C / 1H (ppm):')
for sb in ('89a', '89b'):
    out = []
    for name, t in (('89a', A), ('89b', B)):
        rr = []
        for nuc in 'CH':
            idx = [i for i, r in enumerate(t['rows']) if r[1] == nuc]
            rr.append(line_fit([t['rows'][i][col(t, 'MPW_B_PCM')] for i in idx], [float(blocks[sb][i]['MMFF_PBE']) for i in idx])[0])
        out.append('2023 %s: %.2f / %.3f' % (name, rr[0], rr[1]))
    print('  STX block %s -> %s' % (sb, ' | '.join(out)))
print('experimental data of 89a: error of each candidate after linear scaling, RMS (largest), ppm:')
for src, lev in (('STX', 'MMFF_PBE'), ('STX', 'xTB_PBE'), ('STX', 'xTB_mPW'), ('STX', 'B3L_PBE'), ('2023', 'MPW_B_PCM'), ('2023', 'WB9_C_SMD'), ('2023', 'B3L_A_GAS')):
    out = []
    for cand in ('89a', '89b'):
        for nuc in 'CH':
            idx = [i for i, r in enumerate(A['rows']) if r[1] == nuc]
            if src == 'STX': sig = [float(blocks[cand][i][lev]) for i in idx]
            else:
                t = table_for(cand, lev); sig = [t['rows'][i][col(t, lev)] for i in idx]
            rms, mx = line_fit(sig, [A['rows'][i][2] for i in idx])
            out.append('%s %s %.2f (%.2f)' % (cand, nuc, rms, mx))
    print('  %-4s %-10s %s' % (src, lev, ' | '.join(out)))
