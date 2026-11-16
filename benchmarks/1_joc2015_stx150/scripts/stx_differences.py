"""Recomputed (rows as printed, inputs by rule) against published, every STX150 record with both values, by the cause
of the difference: a defect of the printed data, an input that the SI does not print, or rounding (from stx_records.csv).
Usage: stx_differences.py <out dir>"""
import sys, os, csv, collections
rec = [r for r in csv.DictReader(open(os.path.join(sys.argv[1], 'stx_records.csv'))) if r['published'] != '' and r['recomputed_as_printed'] != '']
SOURCE_MORE = lambda r: int(r['K_source']) > int(r['K_printed'])
def cause(r):
    t = r['reason']
    if 'exchanged' in t: return 'exchanged isomers'
    if 'copied column' in t: return 'copied column'
    if 'share a column' in t: return 'shared column'
    if 'credible fit' in t: return 'candidate without credible fit'
    if r['inputs'].startswith('recovered'): return 'unprinted input recovered'
    if SOURCE_MORE(r) and float(r['recomputed_as_printed']) - float(r['published']) > -10: return 'fewer candidates'
    if 'rounding range' in t and 'outside' not in t: return 'rounding of printed shieldings'
    return 'no cause found'
CA = ['exchanged isomers', 'copied column', 'shared column', 'candidate without credible fit', 'fewer candidates', 'unprinted input recovered', 'rounding of printed shieldings', 'no cause found']
def size(r):
    d = abs(float(r['recomputed_as_printed']) - float(r['published'])); return 'within 1' if d <= 1 else ('1 to 10' if d <= 10 else ('10 to 50' if d <= 50 else 'above 50'))
tab = collections.Counter((size(r), cause(r)) for r in rec)
print('records with both values: %d' % len(rec))
print('%-10s %s  all' % ('', ' | '.join('%s' % c for c in CA)))
for s in ('within 1', '1 to 10', '10 to 50', 'above 50'): print('%-10s %s  %d' % (s, ' | '.join('%*d' % (len(c), tab[(s, c)]) for c in CA), sum(tab[(s, c)] for c in CA)))
DEF = CA[:5]
print('cases touched by a defect of the printed data (first five columns): %d of %d' % (len({r['case'] for r in rec if cause(r) in DEF}), len({r['case'] for r in rec})))
for title, test in (('published 0-1%, recomputed 95% or more', lambda r: float(r['published']) <= 1 and float(r['recomputed_as_printed']) >= 95), ('published 99-100%, recomputed 5% or less', lambda r: float(r['published']) >= 99 and float(r['recomputed_as_printed']) <= 5)):
    sel = [r for r in rec if test(r)]
    print('%s: %d records, %d with a cause found' % (title, len(sel), sum(cause(r) != 'no cause found' for r in sel)))
    for r in sel: print('   %-5s %-9s published %3.0f recomputed %5.1f | %s' % (r['case'], r['level'], float(r['published']), float(r['recomputed_as_printed']), cause(r)))
