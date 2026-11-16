"""Tables 1 to 5 of the paper and Tables S1 to S13 of the Supporting Information, built from the tables of the analysis.

Usage: python paper_tables.py <folder with final_*.csv> <folder with the logs> <folder with exhibit_*.csv> <out tables.json> [folder for CSV copies]

The JSON file holds every table as it is typeset (cells with a light markup: *italic*, _{sub}, ^{sup}). With a fifth
argument each table is also written as plain CSV, together with an index that names the source files of every table."""
import sys, csv, json, re, os, io
final, logs, exhibit, out = sys.argv[1:5]
csv_dir = sys.argv[5] if len(sys.argv) > 5 else None
SOURCES = {
    'T1': 'final_benchmarks.csv', 'T2': 'final_scores.csv; exhibit_set_kinds.csv', 'T3': 'final_designs.csv; logs/01_final_method.txt (workflows below the level)',
    'T4': 'final_stress.csv', 'T5': 'final_kinds.csv; final_scores.csv (intervals of JOC2015 and Cam2017)', 'S1': 'exhibit_sources.csv', 'S2': 'final_units.csv',
    'S3': 'logs/03_checks.txt (labels drawn at random); final_random_halves.csv', 'S4': 'final_reliability.csv', 'S5': 'final_temperature_size.csv; exhibit_temperature_bins.csv',
    'S6': 'final_scores.csv', 'S7': 'final_workflows.csv', 'S8': 'final_conditional.csv', 'S9': 'final_misses.csv', 'S10': 'final_kinds_paired.csv', 'S11': 'final_energy.csv',
    'S12': 'exhibit_couplings.csv', 'S13': 'logs/04_explore_candidates.txt'}
def read(d, name):
    return list(csv.DictReader(open(os.path.join(d, name), encoding='utf-8')))
def pct(x, d=1):
    return ('%.' + str(d) + 'f%%') % (100.0 * float(x))
def iv(lo, hi, d=0, scale=100.0):
    return ('(%.' + str(d) + 'f to %.' + str(d) + 'f)') % (scale * float(lo), scale * float(hi))
def num(x, d=2):
    return ('%.' + str(d) + 'f') % float(x)
def cap(x):
    return x if x.startswith('xTB') else x[0].upper() + x[1:]
BENCH = ['JOC2015', 'STX150', 'DP4-AI', 'Cam2017', 'Cam2019']
GEO = {'B3LYP': 'B3LYP', 'MMFF': 'MMFF', 'xTB': 'xTB', 'PM6': 'PM6', 'PM7': 'PM7', 'HF': 'HF', 'M06-2X': 'M06-2X'}
T = {}

# ---------------------------------------------------------------- Table 1
fb = {(r['benchmark'], r['kind_of_calculation']): r for r in read(final, 'final_benchmarks.csv')}
order = [('JOC2015', 'DP4+, B3LYP geometries', 'DP4+; B3LYP'), ('STX150', 'DP4+, B3LYP geometries', 'DP4+; B3LYP'), ('STX150', 'xTB-DP4+, xTB geometries', 'xTB-DP4+; xTB'),
         ('STX150', 'MM-DP4+, MMFF geometries', 'MM-DP4+; MMFF'), ('DP4-AI', 'DP4, HF or DFT geometries', 'DP4; B3LYP'), ('DP4-AI', 'DP4, MMFF geometries', 'DP4; MMFF'),
         ('Cam2017', 'DP4, MMFF geometries', 'DP4; MMFF'), ('Cam2019', 'DP4, HF or DFT geometries', 'DP4; HF, B3LYP, M06-2X'), ('Cam2019', 'DP4, MMFF geometries', 'DP4; MMFF'),
         ('Cam2019', 'DP4, semi-empirical geometries', 'DP4; PM6, PM7')]
rows, tot, prev = [], dict(wf=0, rec=0, top=0, c=0, cr=0), None
for b, k, mg in order:
    r = fb[(b, k)]
    rec, top, c, cr = int(r['records']), int(r['top_right_records']), int(r['calls_99']), int(r['calls_99_right'])
    rows.append([b if b != prev else '', mg, r['workflows'], '{:,}'.format(rec), r['families'], num(r['mean_candidates'], 1), '%.1f%%' % (100.0 * top / rec),
                 '%.1f%% (%d)' % (100.0 * cr / c, c), num(r['temperature'], 1)])
    prev = b
    tot['wf'] += int(r['workflows']); tot['rec'] += rec; tot['top'] += top; tot['c'] += c; tot['cr'] += cr
assert (tot['wf'], tot['rec'], tot['top'], tot['c'], tot['cr']) == (76, 3100, 2215, 1791, 1566)
rows.append(['All', '', '76', '3,100', '152', '', '%.1f%%' % (100.0 * tot['top'] / tot['rec']), '%.1f%% (%s)' % (100.0 * tot['cr'] / tot['c'], '{:,}'.format(tot['c'])), ''])
T['T1'] = dict(label='Table 1', title='The five benchmarks.',
    caption=' Number of workflows, records and families, mean number of candidates per record, share of the records whose top candidate is the actual structure, share of the top candidates stated at 99% or more that are the actual structure (with the number of such calls), and the temperature *T*.',
    head=['Benchmark', 'Method; geometries', 'Workflows', 'Records', 'Families', 'Candidates', 'Top candidate right', 'Stated ≥ 99%: right (calls)', '*T*'],
    rows=rows, widths=[1050, 1790, 1000, 760, 780, 1040, 900, 1186, 520], align='llccccccc', lastrule=True,
    note='Sources: JOC2015 {{c:grimblat2015}}, STX150 {{c:luciano2026}}, DP4-AI {{c:howarth2020}}, Cam2017 {{c:ermanis2017}} and Cam2019 {{c:ermanis2019}}. Candidates are averaged per family within a workflow and then over these units; the shares refer to records. *T* is the factor by which all log-likelihoods have to be divided for the stated probabilities to fit the outcomes (Section 2.4); a value of 1 means calibrated probabilities. The last row counts every family once; families are shared between JOC2015 and STX150 and between Cam2017 and Cam2019 (Section 2.2).')

# ---------------------------------------------------------------- Table 2
sc = read(final, 'final_scores.csv')
SK = {(r['benchmark'], '0.' + r['target'][0]): r for r in read(exhibit, 'exhibit_set_kinds.csv')}
def S(target, method, b):
    rr = [r for r in sc if r['target'] == target and r['method'] == method and r['benchmark'] == b]
    assert len(rr) == 1, (target, method, b)
    return rr[0]
rows = []
for b in BENCH:
    for t in ('0.9', '0.8'):
        c = S(t, 'likelihood ratio per shift', b); s = S(t, 'stated probabilities add up to the target', b); f = S(t, 'the same after fitting a temperature', b)
        size = num(c['mean_set_size']) + (' (%s to %s)' % (num(c['size_lo']), num(c['size_hi'])) if float(c['size_hi']) - float(c['size_lo']) > 0.004 else '')
        rows.append([b if t == '0.9' else '', '%d%%' % round(100 * float(t)), num(c['mean_candidates'], 1), size,
                     '%.0f%% (%.0f%%)' % (100 * float(SK[(b, t)]['one_structure']), 100 * float(c['single_structure_sets_right'])),
                     '%.0f%%' % (100 * float(SK[(b, t)]['every_candidate'])),
                     '%s; %s' % (pct(s['coverage']), num(s['mean_set_size'])), '%s; %s' % (pct(f['coverage']), num(f['mean_set_size']))])
T['T2'] = dict(label='Table 2', title='Calibrated sets and their price.',
    caption=' Sets of the likelihood ratio per shift with one family left out: mean number of structures with its 95% interval, share of sets with one structure (and how often it is the actual one), share of sets with every candidate. For comparison, coverage and mean size of the smallest set whose stated probabilities add up to the level, before and after a temperature fitted on the calibration families.',
    head=['Benchmark', 'Level', 'Candidates', 'Structures in the set', 'One structure (right)', 'Every candidate', 'Stated set: covers; size', 'After fitted *T*: covers; size'],
    rows=rows, widths=[980, 600, 960, 1700, 1190, 940, 1320, 1336], align='lccccccc',
    note='All values are averages per family (Section 2.4). The share of right single structures rests on few sets where single structures are rare (Cam2019 at 90%).')

# ---------------------------------------------------------------- Table 3
de = read(final, 'final_designs.csv')
low = {}
for line in open(os.path.join(logs, '01_final_method.txt'), encoding='utf-8'):
    m = re.match(r'\s+(\d+)%\s+(.+?)\s+lowest\s+([\d.]+)% .*below target by more than 5:\s+(\d+) of (\d+)', line)
    if m:
        low[(m.group(1), m.group(2).strip())] = (m.group(3), m.group(4), m.group(5))
rows = []
names = [('own workflow', 'Own workflow'), ('whole benchmark', 'Whole benchmark'), ('same kind of calculation', 'Same kind of calculation, all benchmarks'), ('everything', 'All records'), ('other benchmarks', 'Other benchmarks only')]
for key, label in names:
    for t in ('0.9', '0.8'):
        cells = []
        for b in BENCH:
            rr = [r for r in de if r['target'] == t and r['cutoff_from'] == key and r['benchmark'] == b]
            cells.append(pct(rr[0]['coverage']))
        lo = low[(str(round(100 * float(t))), key)]
        assert lo[2] == '76'
        rows.append([label if t == '0.9' else '', '%d%%' % round(100 * float(t))] + cells + [lo[1], lo[0] + '%'])
T['T3'] = dict(label='Table 3', title='Where the cut-off may come from.',
    caption=' Coverage of the sets when the cut-off is calibrated on the own workflow, on all workflows of the benchmark, on all workflows of the same kind of calculation in every benchmark, on all records, or on the other benchmarks only. The family of the record is left out in every case. The last columns give the number of the 76 workflows whose coverage is more than five points below the level, and the lowest coverage of a workflow.',
    head=['Cut-off from', 'Level', 'JOC2015', 'STX150', 'DP4-AI', 'Cam2017', 'Cam2019', 'Workflows below', 'Lowest'],
    rows=rows, widths=[2300, 620, 860, 860, 860, 860, 860, 986, 820], align='lcccccccc',
    note='The first origin reaches the level by construction (Section 2.4). The intervals of these values are in the file final_designs.csv of the repository.')

# ---------------------------------------------------------------- Table 4
st = read(final, 'final_stress.csv')
def stress(cal, tst, t, score):
    rr = [r for r in st if r['calibrated_on'] == cal and r['tested_on'] == tst and r['target'] == t and r['score'] == score]
    assert len(rr) == 1, (cal, tst, t, score, len(rr))
    return rr[0]
pairs = []
for r in st:
    k = (r['calibrated_on'], r['tested_on'])
    if k not in pairs:
        pairs.append(k)
want = ['size', 'candidates', 'collection (STX150)', 'collection (Cam2017)']
pairs.sort(key=lambda k: want.index(next(r['test'] for r in st if (r['calibrated_on'], r['tested_on']) == k)))
rows = []
for cal, tst in pairs:
    for t in ('0.9', '0.8'):
        a = stress(cal, tst, t, 'likelihood ratio per shift'); b = stress(cal, tst, t, 'likelihood ratio')
        rows.append([(cap(cal) + ' → ' + tst) if t == '0.9' else '', '%d%%' % round(100 * float(t)), a['test_families'],
                     '%s %s' % (pct(a['coverage']), iv(a['coverage_lo'], a['coverage_hi'])), '%s %s' % (pct(b['coverage']), iv(b['coverage_lo'], b['coverage_hi']))])
T['T4'] = dict(label='Table 4', title='Stress tests.',
    caption=' Coverage when the cut-off of a workflow is calibrated on one part of its families and applied to a part that differs in a known way, for the likelihood ratio per shift and for the undivided ratio (95% intervals from resampling families).',
    head=['Calibrated on → tested on', 'Level', 'Families tested', 'Ratio per shift', 'Undivided ratio'],
    rows=rows, widths=[4000, 620, 1000, 1703, 1703], align='lcccc',
    note='Only workflows with enough calibration families for the level take part, so the number of families tested differs between the levels.')

# ---------------------------------------------------------------- Table 5
kd = read(final, 'final_kinds.csv')
rows, prev = [], None
for r in kd:
    def sz(t):
        if not r['mean_set_size_' + t]:
            return '–'
        if r['benchmark'] in ('JOC2015', 'Cam2017'):                       # one kind of calculation: the interval of Table 2
            c = S({'95': '0.95', '90': '0.9', '80': '0.8'}[t], 'likelihood ratio per shift', r['benchmark'])
            assert abs(float(c['mean_set_size']) - float(r['mean_set_size_' + t])) < 5e-4
            return num(c['mean_set_size']) + ' (%s to %s)' % (num(c['size_lo']), num(c['size_hi']))
        return num(r['mean_set_size_' + t]) + (' (%s to %s)' % (num(r['size_lo_' + t]), num(r['size_hi_' + t])) if float(r['size_hi_' + t]) - float(r['size_lo_' + t]) > 0.004 else '')
    rows.append([r['benchmark'] if r['benchmark'] != prev else '', '%s; %s' % (r['method'], r['geometry']), r['families'], pct(r['top_candidate_right'], 0),
                 sz('95'), sz('90'), sz('80'), num(r['median_cutoff_per_shift_90']), num(r['median_cutoff_per_shift_80'])])
    prev = r['benchmark']
T['T5'] = dict(label='Table 5', title='Kinds of calculation after calibration.',
    caption=' Share of right top candidates, mean number of structures in the sets at three levels (95% interval), and the median over the workflows of the cut-off per shift.',
    head=['Benchmark', 'Method; geometries', 'Families', 'Top right', 'Structures at 95%', 'Structures at 90%', 'Structures at 80%', '*ĉ* at 90%', '*ĉ* at 80%'],
    rows=rows, widths=[900, 1350, 800, 700, 1450, 1450, 1450, 463, 463], align='llccccccc',
    note='Averages per family (Section 2.4). A dash marks a level for which the benchmark has too few families. A cut-off of 0.00 means that the set is the top candidate alone.')
assert all(sum(t['widths']) == 9026 for t in T.values()), {k: sum(t['widths']) for k, t in T.items()}

# ================================================================ Supporting Information
# S1 sources
src = read(exhibit, 'exhibit_sources.csv')
T['S1'] = dict(label='Table S1', title='Sources and reproduction of the benchmarks.', caption=' Data, rule of the recomputation, agreement with the published numbers, and records used.',
    head=['Benchmark', 'Data', 'Recomputed as', 'Agreement with the published numbers', 'Records used'],
    rows=[[r['benchmark'] + ': ' + r['paper'], cap(r['data']), cap(r['recomputed_as']), cap(r['agreement_with_the_published_numbers']), cap(r['records_used'])] for r in src],
    widths=[1500, 2000, 1800, 2100, 1626], align='lllll', size=16)
# S2 units
un = read(final, 'final_units.csv')
rows = []
for cut, left, label in (('families weighted', 'one family', 'Families weighted, one family left out (the method)'), ('records unweighted', 'one family', 'Records unweighted, one family left out'), ('records unweighted', 'one record', 'Records unweighted, one record left out')):
    for t in ('0.9', '0.8'):
        cells = []
        for b in BENCH:
            r = next(x for x in un if x['target'] == t and x['cutoff'] == cut and x['left_out'] == left and x['benchmark'] == b)
            cells.append('%s; %s' % (pct(r['coverage_by_family']), num(r['mean_set_size'])))
        rows.append([label if t == '0.9' else '', '%d%%' % round(100 * float(t))] + cells)
T['S2'] = dict(label='Table S2', title='Families or records as the units.', caption=' Coverage averaged per family and mean number of structures in the sets (likelihood ratio per shift, cut-off from the own workflow) for three ways of calibrating.',
    head=['Calibration', 'Level'] + BENCH, rows=rows, widths=[2666, 610, 1150, 1150, 1150, 1150, 1150], align='lcccccc')
# S3 checks
rl = {}
for line in open(os.path.join(logs, '03_checks.txt'), encoding='utf-8'):
    m = re.match(r'\s+(\d+)%: (JOC2015 .*)', line)
    if m:
        for part in m.group(2).split('|'):
            q = re.match(r'\s*(\S+)\s+([\d.]+)% \(lowest ([\d.]+)%\), kept (\d+)%', part)
            rl[(m.group(1), q.group(1))] = (q.group(2), q.group(3), q.group(4))
rh = read(final, 'final_random_halves.csv')
rows = []
for t in ('90', '80'):
    rows.append(['Labels drawn at random' if t == '90' else '', t + '%'] +
                ['%s%% (%s%%); %s%%' % rl[(t, b)] for b in BENCH])
for t in ('0.9', '0.8'):
    cells = []
    for b in BENCH:
        r = next(x for x in rh if x['target'] == t and x['score'] == 'likelihood ratio per shift' and x['benchmark'] == b)
        cells.append('%s %s; %s' % (pct(r['coverage']), iv(r['coverage_5'], r['coverage_95']), num(r['mean_set_size'])))
    rows.append(['Families split at random' if t == '0.9' else '', '%d%%' % round(100 * float(t))] + cells)
T['S3'] = dict(label='Table S3', title='Checks of the procedure.', caption=' Coverage with meaningless labels, and coverage when the families are split at random (likelihood ratio per shift).',
    head=['Check', 'Level'] + BENCH, rows=rows, widths=[1676, 550, 1360, 1360, 1360, 1360, 1360], align='lcccccc', size=16, fixed=True,
    note='Labels drawn at random: the actual structure of every record is replaced by a candidate drawn at random (20 draws) and one family is left out; coverage of the drawn label, with the lowest draw in brackets, and share of the candidates kept. Families split at random: 500 splits into halves, with the cut-off from one half applied to the other; coverage, with the range from 5% to 95% of the splits in brackets, and mean number of structures in the set.')
# S4 reliability
re_ = read(final, 'final_reliability.csv')
rows, prev = [], None
for r in re_:
    k = r['kind_of_calculation']
    rows.append([(k if k != 'all' else 'All records') if k != prev else '', r['stated_probability_of_top_candidate'], r['records'], r['families'], pct(r['mean_stated'], 2 if float(r['mean_stated']) > 0.99 else 1), r['right'],
                 '%s %s' % (pct(r['share_right']), iv(r['share_right_lo'], r['share_right_hi'], 1))])
    prev = k
T['S4'] = dict(label='Table S4', title='Reliability of the stated probabilities.', caption=' Records by the stated probability of their top candidate: number of records and families, mean stated probability, number and share of top candidates that are the actual structure (95% interval from resampling families).',
    head=['Kind of calculation', 'Stated probability', 'Records', 'Families', 'Mean stated', 'Right', 'Share right'], rows=rows, widths=[2500, 1700, 800, 800, 1000, 700, 1526], align='llccccc')
# S5 temperature and size
ts = read(final, 'final_temperature_size.csv')
tb = {}
for r in read(exhibit, 'exhibit_temperature_bins.csv'):
    tb.setdefault(r['benchmark'], []).append(r)
def thirds(name):
    return '; '.join(num(r['temperature'], 1) for r in tb[name]) if name in tb else '–'
def med(v):
    return ('%.1f' % float(v)).replace('.0', '')
T['S5'] = dict(label='Table S5', title='Temperature and size of the molecule.', caption=' Temperature of the families of a benchmark in three groups of equal size by number of shifts, from the fewest to the most shifts; exponent *b* of *T* = *T*_{30}(*N*/30)^i{b} with one *T*_{30} per workflow (95% interval from resampling families; the grid runs from −1 to 3); share of the resamples with a positive exponent; and gain in log-likelihood over *b* = 0 for the fitted exponent and for the exponents 1/2 and 1.',
    head=['Workflows', 'Number', 'Families', '*T* by thirds', 'Exponent', 'Resamples above zero', 'Gain, fitted *b*', 'Gain, *b* = 1/2', 'Gain, *b* = 1'],
    rows=[[cap(r['workflows']), r['number'], r['families'], thirds(r['workflows']), '%s (%s to %s)' % (num(r['exponent']), num(r['lo']), num(r['hi'])), pct(r['share_of_resamples_above_zero']),
           num(r['gain_in_log_likelihood_over_b0']), num(r['gain_of_b_half']), num(r['gain_of_b_one'])] for r in ts],
    widths=[2026, 700, 750, 1300, 1650, 900, 600, 550, 550], align='lcccccccc',
    note='Median number of shifts of the three groups: ' + '; '.join('%s, %s, %s and %s' % ((b,) + tuple(med(r['median_shifts']) for r in tb[b])) for b in BENCH) + '.')
# S6 scores and baselines
meth = [('likelihood ratio per shift', 'Likelihood ratio per shift (the method)'), ('likelihood ratio / sqrt(N)', 'Likelihood ratio over the square root of *N*'), ('likelihood ratio', 'Likelihood ratio, undivided'),
        ('stated probability', 'Stated probability'), ('cumulative probability', 'Cumulative probability'), ('cumulative probability, fitted temperature', 'Cumulative probability after a fitted temperature'),
        ('rank', 'Rank'), ('stated probabilities add up to the target', 'No guarantee: stated probabilities add up to the level'), ('the same after fitting a temperature', 'No guarantee: the same after a fitted temperature'),
        ('top candidate alone', 'No guarantee: top candidate alone')]
rows = []
for t in ('0.9', '0.8'):
    for key, label in meth:
        rows.append(['%d%%' % round(100 * float(t)) if key == meth[0][0] else '', label] + ['%s; %s' % (pct(S(t, key, b)['coverage']), num(S(t, key, b)['mean_set_size'])) for b in BENCH])
T['S6'] = dict(label='Table S6', title='Scores and baselines.', caption=' Coverage and mean number of structures in the sets for the conformal scores (one family left out, cut-off from the own workflow) and for three constructions without a guarantee.',
    head=['Level', 'Score'] + BENCH, rows=rows, widths=[560, 3016, 1090, 1090, 1090, 1090, 1090], align='clccccc')
# S7 workflows
wf = read(final, 'final_workflows.csv')
rows, prev = [], None
for r in wf:
    k = (r['benchmark'], r['method'], r['geometry'])
    c = lambda t: (num(r['cutoff_per_shift_' + t], 3) if r['cutoff_per_shift_' + t] not in ('', 'inf') else '–')
    z = lambda t: (num(r['mean_set_size_' + t]) if r['mean_set_size_' + t] else '–')
    rows.append([r['benchmark'] if r['benchmark'] != (prev[0] if prev else None) else '', ('%s; %s' % (r['method'], r['geometry'])) if k != prev else '', r['level'], r['families'], r['median_shifts'], pct(r['top_candidate_right'], 0), num(r['temperature'], 1),
                 c('95'), c('90'), c('80'), z('90'), z('80')])
    prev = k
T['S7'] = dict(label='Table S7', title='Cut-offs of the 76 workflows.', caption=' Families, median number of shifts, share of right top candidates (average per family), temperature, cut-off per shift at three levels from all families of the workflow, and mean number of structures in the sets with one family left out. Levels are labelled as in the sources.',
    head=['Benchmark', 'Method; geometries', 'Level', 'Families', 'Shifts', 'Top right', '*T*', '*ĉ* 95%', '*ĉ* 90%', '*ĉ* 80%', 'Set 90%', 'Set 80%'],
    rows=rows, widths=[850, 1250, 1900, 620, 520, 620, 520, 600, 600, 600, 473, 473], align='lllccccccccc', size=16)
# S8 coverage by class
cd = [r for r in read(final, 'final_conditional.csv') if r['score'] == 'likelihood ratio per shift']
rows, prev = [], None
for r in [x for x in cd if x['target'] == '0.9']:
    o = next(x for x in cd if x['target'] == '0.8' and x['split_by'] == r['split_by'] and x['class'] == r['class'])
    sp = {'record without data problems': 'Record without a data problem'}.get(r['split_by'], cap(r['split_by']))
    rows.append([sp if sp != prev else '', cap(r['class']), '{:,}'.format(int(r['records'])), r['families'], num(r['mean_candidates'], 1),
                 '%s %s' % (pct(r['coverage']), iv(r['coverage_lo'], r['coverage_hi'])), num(r['mean_set_size']), '%s %s' % (pct(o['coverage']), iv(o['coverage_lo'], o['coverage_hi'])), num(o['mean_set_size'])])
    prev = sp
T['S8'] = dict(label='Table S8', title='Coverage by class.', caption=' Coverage (95% interval) and mean number of structures in the sets with one family left out, by class of record. The class of a set with one structure, several or every candidate is taken at the respective level, so that its records differ between the levels; the numbers of records shown are those of the 90% level.',
    head=['Split by', 'Class', 'Records', 'Families', 'Candidates', 'Coverage at 90%', 'Set', 'Coverage at 80%', 'Set'], rows=rows, widths=[1750, 1900, 700, 700, 850, 1300, 263, 1300, 263], align='llccccccc')
T['S8']['widths'] = [1700, 1850, 680, 680, 800, 1330, 328, 1330, 328]
# S9 misses
ms = read(final, 'final_misses.csv')
def fam(s):
    return re.sub(r'\b(dp|ai|cam):', '', s)
T['S9'] = dict(label='Table S9', title='Who is missed.', caption=' Records whose set does not contain the actual structure (one family left out), families never missed, families missed in half of their records or more, share of the misses carried by the three families with most misses, and these families with their missed and total records. Families of JOC2015 and STX150 carry the compound numbers of the sources.',
    head=['Level', 'Benchmark', 'Families', 'Records missed', 'Never missed', 'Missed in half or more', 'Share in three families', 'The three families'],
    rows=[['%d%%' % round(100 * float(r['target'])) if r['benchmark'] == 'JOC2015' else '', r['benchmark'], r['families'], '%s of %s' % (r['units_missed'], '{:,}'.format(int(r['units']))), r['families_never_missed'],
           r['families_missed_in_half_their_records_or_more'], pct(r['share_of_misses_in_the_three_families_with_most'], 0), fam(r['families_with_most_misses'])] for r in ms],
    widths=[560, 900, 700, 1100, 700, 900, 900, 3266], align='clcccccl')
# S10 kinds paired
kp = read(final, 'final_kinds_paired.csv')
rows, prev = [], None
for r in kp:
    rows.append([r['benchmark'] if r['benchmark'] != prev else '', r['kind_a'].replace(' ', '; ', 1), r['kind_b'].replace(' ', '; ', 1), '%d%%' % round(100 * float(r['target'])), r['families'], num(r['size_a']), num(r['size_b']),
                 '%+.2f (%+.2f to %+.2f)' % (float(r['difference']), float(r['lo']), float(r['hi']))])
    prev = r['benchmark']
T['S10'] = dict(label='Table S10', title='Kinds of calculation compared.', caption=' Mean number of structures in the sets of two kinds of calculation on the families that both have, and their difference A − B (95% interval from resampling these families).',
    head=['Benchmark', 'A', 'B', 'Level', 'Families', 'Set A', 'Set B', 'Difference'], rows=rows, widths=[1000, 1700, 1700, 620, 800, 700, 700, 1806], align='lllccccc')
# S11 energy
en = read(final, 'final_energy.csv')
T['S11'] = dict(label='Table S11', title='Uncertain conformer energies.', caption=' Sets from the plain likelihood and from the likelihood averaged over a noise of 4 kJ mol^{−1} on the conformer energies, on the 1,108 DP4 records with conformer data (one family left out).',
    head=['Level', 'Benchmark', 'Records', 'Families', 'Coverage: plain; averaged', 'Set: plain; averaged', 'Difference', 'One structure: plain; averaged', 'Families with smaller; larger sets'],
    rows=[['%d%%' % round(100 * float(r['target'])) if r['benchmark'] == 'DP4-AI' else '', cap(r['benchmark']), '{:,}'.format(int(r['records'])), r['families'], '%s; %s' % (pct(r['coverage_plain']), pct(r['coverage_averaged'])),
           '%s; %s' % (num(r['size_plain']), num(r['size_averaged'])), '%+.2f (%+.2f to %+.2f)' % (float(r['difference']), float(r['lo']), float(r['hi'])),
           '%s; %s' % (pct(r['single_plain'], 0), pct(r['single_averaged'], 0)), '%s; %s' % (r['families_smaller'], r['families_larger'])] for r in en],
    widths=[560, 900, 700, 700, 1500, 1200, 1500, 1066, 900], align='clccccccc')
# S12 couplings
cp = read(exhibit, 'exhibit_couplings.csv')
rows = []
for r in [x for x in cp if x['target'] == '0.9']:
    o = next(x for x in cp if x['target'] == '0.8' and x['data'] == r['data'])
    d = lambda x: ('' if x['data'] == 'shifts' else '%+.2f (%+.2f to %+.2f)' % (float(x['difference_of_mean_set_size']), float(x['diff_lo']), float(x['diff_hi'])))
    rows.append([cap(r['data']), r['actual_first'], '%s (%s to %s)' % (num(r['mean_set_size']), num(r['size_lo']), num(r['size_hi'])), d(r), '%s (%s to %s)' % (num(o['mean_set_size']), num(o['size_lo']), num(o['size_hi'])), d(o)])
T['S12'] = dict(label='Table S12', title='Proton–proton coupling constants.', caption=' The 33 compounds of the Cam2017 deposit with calculated couplings (22 families, B3LYP/6-31G(d,p) on MMFF geometries): records with the actual structure first, mean number of structures in the sets (95% interval) and difference to the same data without couplings. The number of couplings is added to the number of shifts in the score.',
    head=['Data', 'Actual first', 'Set at 90%', 'Difference', 'Set at 80%', 'Difference'], rows=rows, widths=[2800, 700, 1400, 1400, 1326, 1400], align='lccccc',
    note='Coupling models: the linear correction and the Gaussian of the PyDP4 source for the total coupling (standard deviation 0.79 Hz); the same for the Fermi-contact part of the coupling, with the parameters of the source; and the scaling and the error density recovered from the terms printed in the deposited outputs, which has heavier tails. The differences of the rows with averaged energies refer to the shifts with averaged energies, except for "Shifts, energies averaged", which is compared with the plain shifts.')
# S13 number of candidates in the score
rows, t = [], None
lab = {'ratio': 'Ratio', 'ratio / log2(K)': 'Ratio / log_{2}*K*', 'ratio / sqrt(N)': 'Ratio / √*N*', 'ratio / sqrt(N) / log2(K)': 'Ratio / √*N* / log_{2}*K*', 'ratio / N': 'Ratio / *N* (the method)', 'ratio / N / log2(K)': 'Ratio / *N* / log_{2}*K*'}
for line in open(os.path.join(logs, '04_explore_candidates.txt'), encoding='utf-8'):
    m = re.match(r'target (\d+)%', line)
    if m:
        t = m.group(1); first = True; continue
    parts = [p.strip() for p in line.rstrip('\n').split('|')]
    if len(parts) == 9 and t:
        cells = []
        for p in parts[1:8]:
            q = re.match(r'.*:\s+([\d.]+)%\s+([\d.]+)$', p)
            cells.append('%s%%; %s' % (q.group(1), q.group(2)))
        rows.append([t + '%' if first else '', lab[parts[0]]] + cells)
        first = False
assert len(rows) == 12
T['S13'] = dict(label='Table S13', title='The number of candidates in the score.', caption=' Coverage and mean number of structures in the sets (one family left out) when the likelihood ratio is divided by powers of the number of shifts *N* and by the logarithm of the number of candidates *K*, by class of record.',
    head=['Level', 'Score', '*K* = 2', '*K* = 3 or 4', '*K* = 5 to 8', '*K* ≥ 9', '*N* < 25', '*N* = 25 to 39', '*N* ≥ 40'], rows=rows, widths=[560, 1900, 938, 938, 938, 938, 938, 938, 938], align='clccccccc')
T['T5']['size'] = 17
T['S5']['size'] = 16
T['S2']['note'] = 'Coverage is the share of covered records, averaged per family within a workflow and then over these units.'
T['S7']['note'] = 'The cut-offs are those for use with a new molecule; the sets with one family left out use cut-offs from the remaining families. A dash marks a level for which the workflow has too few families. Level labels: JOC2015, functional, basis set and solvent model of the NMR calculation; STX150, geometries and NMR functional as abbreviated in our tables [full names to be copied from the Supporting Information of the source]; DP4-AI, Cam2017 and Cam2019, the labels of the deposits.'
assert all(sum(t['widths']) == 9026 for t in T.values()), {k: sum(t['widths']) for k, t in T.items() if sum(t['widths']) != 9026}
def plain(x, literal=False):
    """a cell or a caption without the markup of the typesetting"""
    if not literal:
        x = x.replace('**', '').replace('*', '')
    x = re.sub(r'\{\{[cx]:([^}]*)\}\}', r'\1', x)
    return re.sub(r'([_^])i?\{([^}]*)\}', r'\1\2', x).strip()
assert set(SOURCES) == set(T)
if csv_dir:
    if not os.path.isdir(csv_dir):
        os.makedirs(csv_dir)
    index = [['table', 'file', 'title', 'caption', 'note', 'built from']]
    for k, t in T.items():
        name = t['label'].replace(' ', '_') + '.csv'
        lit = k in ('S1', 'S7')                            # their asterisks belong to the names of basis sets
        with io.open(os.path.join(csv_dir, name), 'w', encoding='utf-8', newline='') as f:
            w = csv.writer(f, lineterminator='\n')
            w.writerow([plain(h) for h in t['head']])
            w.writerows([[plain(c, lit) for c in r] for r in t['rows']])
        index.append([t['label'], name, plain(t['title']), plain(t['caption']), plain(t.get('note', '')), SOURCES[k]])
    with io.open(os.path.join(csv_dir, 'index.csv'), 'w', encoding='utf-8', newline='') as f:
        csv.writer(f, lineterminator='\n').writerows(index)
for k in ('S1', 'S7'):                                 # literal asterisks of basis sets; elsewhere they are markup
    T[k]['rows'] = [[c.replace('*', '\\*') for c in r] for r in T[k]['rows']]
json.dump(T, open(out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
for k, t in T.items():
    print('%-4s %-48s %2d rows x %2d columns' % (k, t['title'][:48], len(t['rows']), len(t['head'])))
