"""Quality control: do the energy and the shieldings of each conformer belong to the same geometry?
Compares the principal moments of inertia stored for every Gaussian job (NMR job vs. first and last geometry of the
optimisation vs. single point) and the dates of the jobs. Writes per-conformer flags for the DP4 run.
Usage: qc_geometry.py <res/core> <work_dir>"""
import sys, os, re, glob, gzip, json, collections, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dp4ai_data as A
core, work = sys.argv[1:3]
J = collections.defaultdict(dict)
for f in sorted(glob.glob(os.path.join(core, 'gaussian_dp4ai_*.jsonl.gz'))):
    for line in gzip.open(f, 'rt'):
        r = json.loads(line); m = A.P1.match(r['path'])
        if not m: continue
        d = m.groupdict()
        if d['sub'] not in (None, 'nmr', 'opt', 'e'): continue
        kind = 'nmr' if 'iso' in r else ('e' if r.get('scf_m') == 'RM062X' else ('opt' if r.get('scf_m') == 'RB3LYP' else None))
        if kind: J[(d['zip'], d['level'], int(d['iso']), int(d['conf']))][kind] = r
def same(a, b):                                   # moments are stored with 3 decimals; 0.01% is far below any real change of geometry
    return a is not None and b is not None and max(abs(x - y) / max(abs(x), 1e-9) for x, y in zip(a, b)) < 1e-4
def when(r):
    m = re.search(r'at \w{3} (\w{3})\s+(\d+) (\d\d:\d\d:\d\d) (\d{4})', r.get('tail', '') or '')
    try: return datetime.datetime.strptime('%s %s %s %s' % m.groups(), '%b %d %H:%M:%S %Y') if m else None
    except ValueError: return None
stat = collections.Counter(); flags = {}; per_case = collections.defaultdict(collections.Counter)
for (case, level, iso, conf), d in sorted(J.items()):
    n = d.get('nmr')
    if n is None or level == 'MM': continue
    o, e = d.get('opt'), d.get('e')
    geom_ok = o is not None and same(n.get('pm'), o.get('pm_last') or o.get('pm'))
    fl = dict(geom_ok=bool(geom_ok))
    stat['%-5s NMR job on the final optimised geometry' % level if geom_ok else '%-5s NMR job NOT on the final optimised geometry' % level] += 1
    if not geom_ok and o is not None:
        per_case[(case, level)]['NMR not on optimised geometry'] += 1
        if same(n.get('pm'), o.get('pm')): stat['%-5s    ... it is on the starting geometry of the optimisation' % level] += 1
        mm = J.get((case, 'MM', iso, conf), {}).get('nmr')
        if mm is not None and max(abs(a - b) for a, b in zip(mm['iso'], n['iso'])) < 5e-4: stat['%-5s    ... and its shieldings equal those in the MM folder' % level] += 1
        tn, to = when(n), when(o)
        if tn and to and tn < to: stat['%-5s    ... and it finished before the optimisation did' % level] += 1
    if level == 'Opt_E':
        e_ok = e is not None and same(n.get('pm'), e.get('pm'))
        fl['e_ok'] = bool(e_ok)
        stat['Opt_E single point on the geometry of the NMR job' if e_ok else 'Opt_E single point NOT on the geometry of the NMR job'] += 1
        if not e_ok: per_case[(case, level)]['single point on another geometry'] += 1
        if e is not None and o is not None and not same(e.get('pm'), o.get('pm_last') or o.get('pm')): stat['Opt_E single point not on the final optimised geometry'] += 1
    per_case[(case, level)]['conformers'] += 1
    flags['%s|%s|%d|%d' % (case, level, iso, conf)] = fl
for k, v in sorted(stat.items()): print('%-66s %6d' % (k, v))
print('\ncases affected (conformers with a problem / conformers):')
for (case, level), c in sorted(per_case.items()):
    if len(c) > 1: print('   %-6s %-5s %s' % (case, level, ', '.join('%s %d/%d' % (k, v, c['conformers']) for k, v in c.items() if k != 'conformers')))
json.dump(flags, open(os.path.join(work, 'qc_flags.json'), 'w'))
