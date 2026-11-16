"""How much does the class of one carbon (sp2 or sp3 parameters) change a recomputed STX150 probability?
The STX150 SI does not print the classes; the recomputation takes every carbon at or above 100 ppm as sp2.
Usage: stx_sp2_sensitivity.py <cs1 scripts dir> <cs1 out dir> <case> [label ...]   (labels of carbons to flip)"""
import sys
sys.path.insert(0, sys.argv[1])
import stx_data as S
blocks, cases, par = S.load(sys.argv[2])
case = sys.argv[3]; flips = sys.argv[4:]
cands = cases[case]['cands']
rows = blocks[case]['rows']
high = [(r['label'], r['exp']) for r in rows if r['nuc'] == 'C' and r['exp'] is not None and r['exp'] >= 85]
print('case %s, candidates %s | carbons at 85 ppm or above: %s' % (case, ' '.join(cands), high))
def run(ov):
    out = []
    for lev in S.LEVELS:
        p, terms, n = S.dp4plus_case(blocks, case, cands, lev, par, sp2_overrides=ov)
        out.append('%s=%.0f' % (lev, 100 * p[('DP4+', 'All')][cands.index(case)]) if p is not None else lev + '=-')
    return ' '.join(out)
print('   threshold rule          :', run(None))
for lab in flips or [l for l, _ in high]:
    cur = [r['exp'] for r in rows if r['nuc'] == 'C' and r['label'] == lab][0] >= 100
    print('   carbon %-3s as %s        :' % (lab, 'sp3' if cur else 'sp2'), run({('C', lab): not cur}))
