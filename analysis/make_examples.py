"""Example files of the companion program, taken from the pooled table.
Usage: python make_examples.py <pooled_candidates.csv> <folder examples of the program>"""
import sys, csv, io, os, math
pooled, out = sys.argv[1:3]
WF = ('JOC2015', 'DP4+', 'B3LYP', 'B3LYP/6-31G*/gas')
rows = [r for r in csv.DictReader(io.open(pooled, encoding='utf-8')) if (r['benchmark'], r['method'], r['geometry'], r['level']) == WF]
with io.open(os.path.join(out, 'calibration_example.csv'), 'w', encoding='utf-8', newline='') as f:
    f.write(u'molecule,family,shifts,candidate,log_likelihood,actual\n')
    for r in rows:
        f.write(u'%s,%s,%s,%s,%s,%s\n' % (r['case'], r['family'], r['n_nuclei'], r['candidate'], r['log_ratio'], r['is_true']))
case = [r for r in rows if r['case'] == '87a']
ell = [float(r['log_ratio']) for r in case]
p = [math.exp(x) for x in ell]
with io.open(os.path.join(out, 'typed_example.txt'), 'w', encoding='utf-8', newline='') as f:
    f.write(u'# Case 87a of the DP4+ benchmark (J. Org. Chem. 2015), DP4+ at B3LYP/6-31G* in the gas phase, recomputed.\n')
    f.write(u'# shifts: %s\n' % case[0]['n_nuclei'])
    f.write(u'# workflow: %s\n' % ' | '.join(WF))
    for r, x in zip(case, p):
        f.write(u'%s  %.10g\n' % (r['candidate'], 100 * x / sum(p)))
with io.open(os.path.join(out, 'calibration_template.csv'), 'w', encoding='utf-8', newline='') as f:
    f.write(u'molecule,family,shifts,candidate,probability,actual\n'
            u'compound 1,series A,24,isomer 1,97.3312,1\ncompound 1,series A,24,isomer 2,2.6688,0\n'
            u'compound 2,series A,26,isomer 1,0.0412,0\ncompound 2,series A,26,isomer 2,99.9588,1\n'
            u'compound 3,series B,31,isomer 1,64.1,0\ncompound 3,series B,31,isomer 2,35.2,1\ncompound 3,series B,31,isomer 3,0.7,0\n')
print('examples written:', len(rows), 'rows of', ' | '.join(WF), '; case 87a with', len(case), 'candidates')
