"""STX150 validation cases above 121 against the shielding tables of the 2025 SI (obc_tables.py, version 2).
Cases are matched to compounds through the experimental 13C shifts; rows are aligned through nucleus and
experimental value. For the correct isomer and the first alternative printed in the STX150 SI:
(1) is the STX xTB_mPW column the Level D column of the 2025 SI (same numbers)?
(2) does [STX correct - STX alternative] follow [2025 isomer a - isomer b] (r near +1) or its reverse (r near -1)?
    MMFF_PBE is compared with Level C, the xTB columns with Level D, the B3LYP-geometry columns with Level B.
Usage: stx_obc_identify.py <obc_compounds.json> <stx out dir> <out csv>"""
import sys, os, csv, json, math, collections, re
C = json.load(open(sys.argv[1])); stx_out = sys.argv[2]; out_csv = sys.argv[3]
blocks = collections.OrderedDict()
for r in csv.DictReader(open(os.path.join(stx_out, 'stx_val_blocks.csv'))): blocks.setdefault(r['block'], []).append(r)
cases = list(csv.DictReader(open(os.path.join(stx_out, 'stx_val_cases.csv'))))
num = lambda c: int(re.match(r'\d+', c).group())
def key(v): return sorted(round(x, 1) for x in v)
def overlap(k, ks):
    pool = collections.Counter(k); same = 0
    for a in ks:
        if pool.get(a, 0) > 0: pool[a] -= 1; same += 1
    return same
def groups(seq):
    out = []
    for i, k in enumerate(seq):
        if out and out[-1][0] == k: out[-1][1].append(i)
        else: out.append((k, [i]))
    return out
def align(crows, srows):
    gj = groups([(r['nucleus'], round(r['exp'], 2)) for r in crows]); gs = groups([(r['nucleus'], round(float(r['exp']), 2) if r['exp'] != '' else None) for r in srows])
    pairs = []; i = 0
    for k, idx in gj:
        j = i
        while j < len(gs) and gs[j][0] != k: j += 1
        if j == len(gs): continue
        sidx = gs[j][1]; i = j + 1
        if len(sidx) == len(idx): pairs += [(a, [b]) for a, b in zip(idx, sidx)]
        elif len(sidx) % len(idx) == 0:
            m = len(sidx) // len(idx); pairs += [(a, sidx[n * m:(n + 1) * m]) for n, a in enumerate(idx)]
    return pairs
def corr(x, y):
    n = len(x)
    if n < 4: return None
    mx, my = sum(x) / n, sum(y) / n; sxx = sum((a - mx) ** 2 for a in x); syy = sum((b - my) ** 2 for b in y)
    if sxx < 1e-12 or syy < 1e-12: return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)
PAIRS = [('MMFF_PBE', 'C'), ('xTB_PBE', 'D'), ('xTB_mPW', 'D'), ('xTB_wB9', 'D'), ('B3L_PBE', 'B'), ('B3L_wB9', 'B')]
ckey = {n: key([r['exp'] for r in c['rows'] if r['nucleus'] == 'C']) for n, c in C.items()}
out = []
for cs in cases:
    name = cs['case']
    if num(name) <= 121: continue
    rows = blocks[name]; ks = key([float(r['exp']) for r in rows if r['nucleus'] == 'C' and r['exp'] != ''])
    if not ks: out.append(dict(case=name, compound_2025='', note='no experimental data printed in the STX150 SI')); continue
    # several compounds can carry the same experimental column (the 2025 SI repeats that of 30 and 31 for 32 and 33):
    # among equally good matches take the one whose number fits the numbering of the two papers (offset 120 or 121)
    sc, _, n = max((overlap(k, ks) / max(len(k), len(ks)), int(num(name) - int(n) in (120, 121)), n) for n, k in ckey.items())
    if sc < 0.8: out.append(dict(case=name, compound_2025='', note='no table with these experimental values (best overlap %.2f)' % sc)); continue
    c = C[n]; pairs = align(c['rows'], rows); cand = cs['candidates'].split(); alt = [x for x in cand if x != name][0]
    bad = {(iso, l, i) for iso, l, i in c['bad_cells']}
    rec = dict(case=name, compound_2025=int(n), overlap_13C=round(sc, 2), isomers_2025=len(c['isomers']), K_printed_stx=len(cand), rows_paired='%d of %d' % (len(pairs), len(c['rows'])), pair_with=c['pair_with'] or '')
    def s(block, lev, bs):
        vals = [blocks[block][b][lev] for b in bs]
        return None if any(v == '' for v in vals) else sum(float(v) for v in vals) / len(vals)
    if len(blocks[alt]) != len(rows): rec['note'] = 'candidate blocks differ in length'; out.append(rec); continue
    if 'D' in c['levels']:
        for who, blk, iso in (('correct', name, 'a'), ('alternative', alt, 'b')):
            d = [abs(s(blk, 'xTB_mPW', bs) - c['isomers'][iso]['D'][a]) for a, bs in pairs if s(blk, 'xTB_mPW', bs) is not None]
            rec['xTB_mPW_is_level_D_' + who] = '%d of %d' % (sum(x <= 0.006 for x in d), len(d))
    for sl, ol in PAIRS:
        if ol not in c['levels']: continue
        rr = []
        for nuc in 'CH':
            x = []; y = []
            for a, bs in pairs:
                if c['rows'][a]['nucleus'] != nuc or ('a', ol, a) in bad or ('b', ol, a) in bad: continue
                u, v = s(name, sl, bs), s(alt, sl, bs)
                if u is None or v is None: continue
                x.append(u - v); y.append(c['isomers']['a'][ol][a] - c['isomers']['b'][ol][a])
            rr.append(corr(x, y))
        rec['r13C_%s_vs_%s' % (sl, ol)] = '' if rr[0] is None else round(rr[0], 2); rec['r1H_%s_vs_%s' % (sl, ol)] = '' if rr[1] is None else round(rr[1], 2)
    out.append(rec)
keys = []
for r in out:
    for k in r:
        if k not in keys: keys.append(k)
w = csv.DictWriter(open(out_csv, 'w', newline=''), fieldnames=keys); w.writeheader(); w.writerows(out)
m = [r for r in out if r['compound_2025'] != '']
print('STX150 cases above 121: %d | matched to a compound of the 2025 SI: %d | compounds not used by any case: %s' % (len(out), len(m), sorted(set(int(n) for n in C) - {r['compound_2025'] for r in m})))
offs = collections.Counter(num(r['case']) - r['compound_2025'] for r in m); print('STX number minus 2025 number:', dict(offs), '| cases with the smaller offset start at', min((r['case'] for r in m if num(r['case']) - r['compound_2025'] == min(offs)), key=num))
tot = [0, 0, 0, 0]
for r in m:
    for i, w_ in enumerate(('correct', 'alternative')):
        v = r.get('xTB_mPW_is_level_D_' + w_, '')
        if v: a, b = map(int, v.split(' of ')); tot[2 * i] += a; tot[2 * i + 1] += b
print('xTB_mPW of STX150 equal to Level D of the 2025 SI: correct isomer %d of %d values, alternative %d of %d' % tuple(tot))
g = lambda r, k: ('%5s' % r[k]) if r.get(k, '') != '' else '    -'
print('case  2025 | Level D identity (correct; alternative) | r(13C) r(1H): MMFF_PBE~C | xTB_PBE~D | xTB_mPW~D | xTB_wB9~D | B3L_PBE~B | B3L_wB9~B')
for r in out:
    if r['compound_2025'] == '': print('%-5s %s' % (r['case'], r['note'])); continue
    print('%-5s %4s | %-9s %-9s | %s' % (r['case'], r['compound_2025'], r.get('xTB_mPW_is_level_D_correct', '-'), r.get('xTB_mPW_is_level_D_alternative', '-'),
          ' | '.join('%s %s' % (g(r, 'r13C_%s_vs_%s' % p), g(r, 'r1H_%s_vs_%s' % p)) for p in PAIRS)))
