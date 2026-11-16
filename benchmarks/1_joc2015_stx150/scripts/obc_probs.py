"""Org. Biomol. Chem. 2025 SI part I, sections 3-5: probability of the correct isomer for compounds 1-39 at
Level A (MM-DP4+, SMD/wB97XD/6-31+G**//MMFF), Level B (DP4+, PCM/mPW1PW91/6-31+G**//B3LYP/6-31G*) and
Level C (MM-DP4+, PCM/mPW1PW91/6-31+G**//MMFF), against the STX150 validation table for the cases above 121
(columns MMFF_wB9, B3L_mPW, MMFF_mPW). The numbering of the two papers is found from the values themselves.
The probabilities fix the split point only to a range, because the compounds around it are all at 100%; the default,
144, comes from the experimental shifts (stx_obc_identify.py: the 2025 table between those of 143a and 144a has no STX case).
Usage: obc_probs.py <d4ob01959h1_suppl.txt> <stx110_pcorrect.csv> <out csv> [split]"""
import sys, re, csv, collections
lines = open(sys.argv[1], encoding='utf-8', errors='replace').read().split('\n')
ROW = re.compile(r'^\s*(\d{1,2})\s+((?:\d+\.\d+%\s+){8}\d+\.\d+%)\s*$')
sec = None; P = collections.defaultdict(dict)
for i, l in enumerate(lines):
    if i < 60: continue
    s = l.strip()
    if s.startswith('3. MM-DP4+ probabilities'): sec = 'A'
    elif s.startswith('4. DP4+ probabilities'): sec = 'B'
    elif s.startswith('5. MM-DP4+ probabilities'): sec = 'C'
    elif s.startswith('6.'): sec = None
    m = ROW.match(l)
    if m and sec: P[int(m.group(1))][sec] = [float(x.rstrip('%')) for x in m.group(2).split()]   # ..., H full, C full, All full
print('compounds with probabilities at Levels A, B, C:', {k: sum(1 for c in P.values() if k in c) for k in 'ABC'})
stx = collections.OrderedDict((r['case'], r) for r in csv.DictReader(open(sys.argv[2])))
num = lambda c: int(re.match(r'\d+', c).group())
hi = [c for c in stx if num(c) > 121]
def agree(c, k): return k in P and abs(float(stx[c]['MMFF_wB9']) - P[k]['A'][8]) <= 1 and abs(float(stx[c]['B3L_mPW']) - P[k]['B'][8]) <= 1
print('STX cases above 121: %d (%s to %s). Numbering: STX = 2025 compound + 121 below a split point, + 120 from there on' % (len(hi), hi[0], hi[-1]))
scores = {sp: sum(agree(c, num(c) - (121 if num(c) < sp else 120)) for c in hi) for sp in range(122, 162)}
best = max(scores.values()); good = [sp for sp, v in scores.items() if v == best]
print('  cases that agree at Levels A and B: %d of %d with one offset of 121; %d of %d with the split anywhere from %d to %d' % (scores[161], len(hi), best, len(hi), good[0], good[-1]))
split = int(sys.argv[4]) if len(sys.argv) > 4 else 144
assert split in good, 'split point outside the range the probabilities allow'
def n_diff(sp): return sum(abs(float(stx[c]['MMFF_mPW']) - P[num(c) - (121 if num(c) < sp else 120)]['C'][8]) > 1 for c in hi)
print('  split used: %d | MMFF_mPW cases that differ from Level C, for every admissible split: %s' % (split, sorted({n_diff(sp) for sp in good})))
out = []; cnt = collections.Counter()
for c in hi:
    k = num(c) - (121 if num(c) < split else 120)
    for lev, col in (('A', 'MMFF_wB9'), ('B', 'B3L_mPW'), ('C', 'MMFF_mPW')):
        a, b = float(stx[c][col]), P[k][lev][8]; ok = abs(a - b) <= 1.0; cnt[(col, ok)] += 1
        out.append(dict(stx_case=c, compound_2025=k, stx_column=col, level_2025=lev, stx_table=a, si_2025_all_full=b, same=int(ok)))
w = csv.DictWriter(open(sys.argv[3], 'w', newline=''), fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
for col in ('MMFF_wB9', 'B3L_mPW', 'MMFF_mPW'): print('  STX column %-8s equals the 2025 value within 1 point in %d of %d cases' % (col, cnt[(col, True)], cnt[(col, True)] + cnt[(col, False)]))
both = [stx[c] for c in hi]
print('  in the STX table the MMFF_mPW value equals the MMFF_wB9 value in %d of %d rows above 121 (%d of them not 100/100)' % (sum(r['MMFF_mPW'] == r['MMFF_wB9'] for r in both), len(both), sum(1 for r in both if r['MMFF_mPW'] == r['MMFF_wB9'] and float(r['MMFF_mPW']) < 99.5)))
low = [r for r in stx.values() if num(r['case']) <= 121 and r['MMFF_mPW'] != '' and r['MMFF_wB9'] != '' and not (float(r['MMFF_mPW']) >= 99.5 and float(r['MMFF_wB9']) >= 99.5)]
print('  for comparison, rows up to 121 where the two are not both 100: equal in %d of %d' % (sum(r['MMFF_mPW'] == r['MMFF_wB9'] for r in low), len(low)))
d = [o for o in out if o['stx_column'] == 'MMFF_mPW']
a = [o['stx_table'] for o in d]; b = [o['si_2025_all_full'] for o in d]
print('  MMFF_mPW (STX table) against Level C (2025 SI): different in %d of %d cases; correct isomer above 50%%: %d against %d' % (sum(1 - o['same'] for o in d), len(d), sum(x > 50 for x in a), sum(x > 50 for x in b)))
print('  on different sides of 50%% (case table/Level C/Level A):', ' '.join('%s %.0f/%.1f/%.1f' % (o['stx_case'], o['stx_table'], o['si_2025_all_full'], P[o['compound_2025']]['A'][8]) for o in d if (o['stx_table'] > 50) != (o['si_2025_all_full'] > 50)))
