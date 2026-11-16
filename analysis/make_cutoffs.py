"""Cut-offs of the 76 workflows for the companion program, computed with the functions of the program itself.
Usage: python make_cutoffs.py <pooled_candidates.csv> <folder of the program (app)> <cutoffs.csv> [final_workflows.csv to compare with]"""
import sys, csv, io, collections
pooled, app, out = sys.argv[1:4]
sys.path.insert(0, app)
import conformalstereo as cs

groups = collections.OrderedDict()
with io.open(pooled, encoding='utf-8') as f:
    for r in csv.DictReader(f):
        groups.setdefault((r['benchmark'], r['method'], r['geometry'], r['level']), collections.OrderedDict()).setdefault(r['case'], []).append(r)
rows = []
for (bench, method, geo, level), cases in groups.items():
    items, ns, ks, right = [], [], [], 0
    for case, rr in cases.items():
        n = int(rr[0]['n_nuclei'])
        ell = [float(r['log_ratio']) for r in rr]
        k = [int(r['is_true']) for r in rr].index(1)
        items.append(((max(ell) - ell[k]) / n, rr[0]['family']))
        ns.append(n); ks.append(len(rr)); right += ell[k] == max(ell) and ell.index(max(ell)) == k
    row = collections.OrderedDict(id=' | '.join((bench, method, geo, level)), benchmark=bench, method=method, geometry=geo, level=level,
                                  records=len(cases), families=len(set(f for _, f in items)), shifts_min=min(ns), shifts_max=max(ns),
                                  candidates_min=min(ks), candidates_max=max(ks), top_right='%.4f' % (right / float(len(cases))))
    for lv, key in ((0.8, 'cutoff_80'), (0.9, 'cutoff_90'), (0.95, 'cutoff_95')):
        c = cs.cutoff(items, lv)
        row[key] = '' if c == cs.INF else '%.8f' % c
    rows.append(row)
with io.open(out, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator='\n')
    w.writeheader(); w.writerows(rows)
print('%d workflows written to %s' % (len(rows), out))
if len(sys.argv) > 4:
    ref = {(r['benchmark'], r['method'], r['geometry'], r['level']): r for r in csv.DictReader(io.open(sys.argv[4], encoding='utf-8'))}
    worst = 0.0; n = 0
    for r in rows:
        q = ref[(r['benchmark'], r['method'], r['geometry'], r['level'])]
        assert int(q['records']) == r['records'] and int(q['families']) == r['families'], r['id']
        for a, b in (('cutoff_80', 'cutoff_per_shift_80'), ('cutoff_90', 'cutoff_per_shift_90'), ('cutoff_95', 'cutoff_per_shift_95')):
            assert (r[a] == '') == (q[b] == ''), (r['id'], a)
            if r[a] != '':
                worst = max(worst, abs(float(r[a]) - float(q[b]))); n += 1
    print('compared with the tables of the analysis: %d cut-offs, largest difference %.1e' % (n, worst))
    assert worst < 1e-8
