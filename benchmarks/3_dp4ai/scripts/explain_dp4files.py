"""For every original .dp4 file: which shieldings (level) and which conformer energies reproduce its 'calc' column?
The TMS constant is solved from the file itself (one constant per nucleus), so this test does not depend on the reference used.
Usage: explain_dp4files.py <res/core> <work_dir>"""
import sys, os, json, re, glob, gzip, math, collections, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dp4lib as L
core, work = sys.argv[1:3]
R, T, HARTREE = 8.3145, 298.15, 2625.499629554010
P1 = re.compile(r'^(?P<zip>[^!]+)\.zip!/(?P<top>[^/]+)/(?P<level>MM|Opt|Opt_E)/(?:(?P<sub>[^/]+)/)?(?P<stem>[^/]+?)_(?P<iso>\d+)ginp(?P<conf>\d+)\.(?:out|log)$')
E = collections.defaultdict(dict)       # (case, level, kind) -> {(iso, conf): energy}
SIG = collections.defaultdict(dict)     # (case, level) -> {(iso, conf): (elements, sigma)}
for f in sorted(glob.glob(os.path.join(core, 'gaussian_dp4ai_*.jsonl.gz'))):
    for line in gzip.open(f, 'rt'):
        r = json.loads(line); m = P1.match(r['path'])
        if not m: continue
        d = m.groupdict(); k = (int(d['iso']), int(d['conf']))
        if d['sub'] not in (None, 'nmr', 'opt', 'e'): continue
        if 'iso' in r:
            if any(v is None for v in r['iso']): continue
            SIG[(d['zip'], d['level'])][k] = (r['el'].split(), r['iso']); kind = 'nmr'
        elif r.get('scf_m') == 'RM062X': kind = 'e'
        elif r.get('scf_m') == 'RB3LYP': kind = 'opt'
        else: continue
        if r.get('scf') is not None: E[(d['zip'], d['level'], kind)][k] = r['scf']
def averaged(case, level, iso, how):
    confs = sorted(c for (i, c) in SIG[(case, level)] if i == iso)
    if not confs: return None
    if how == 'uniform': en = {c: 0.0 for c in confs}
    else:
        src = E.get((case, level, how)) or {}
        en = {c: src[(iso, c)] for c in confs if (iso, c) in src}
    if not en: return None
    use = sorted(en); e0 = min(en.values())
    w = [math.exp(-(en[c] - e0) * HARTREE * 1000 / (R * T)) for c in use]
    if how == 'lowest': w = [1.0 if en[c] == e0 else 0.0 for c in use]
    q = sum(w); els = SIG[(case, level)][(iso, use[0])][0]
    return els, [sum(wi / q * SIG[(case, level)][(iso, c)][1][a] for wi, c in zip(w, use)) for a in range(len(els))]
S = L.Small(core, 'dp4ai')
rows = []
for p in sorted(x for x in S.paths() if x.endswith('.dp4')):
    m = re.match(r'^([^!]+)\.zip!/[^/]+/(MM|Opt|Opt_E)/([^/]+)$', p)
    if not m: continue
    case, folder, name = m.groups(); d = L.parse_dp4(S.read(p)); best = None; allres = {}
    for level in ('MM', 'Opt', 'Opt_E'):
        for how in ('nmr', 'opt', 'e', 'uniform'):
            tms = {'C': [], 'H': []}; bad = False
            for k, tabs in d['isomers'].items():
                av = averaged(case, level, k, how)
                if av is None: bad = True; break
                els, sig = av
                for nuc in 'CH':
                    for label, calc, scaled, exp, err in tabs.get(nuc, []):
                        a = int(re.sub(r'\D', '', label)) - 1
                        if a >= len(sig) or els[a] != nuc: bad = True; break
                        tms[nuc].append((calc + sig[a]) / (1 + calc * 1e-6))      # reference constant implied by this atom
                    if bad: break
                if bad: break
            if bad or not (tms['C'] or tms['H']): continue
            res = max(max(abs(x - statistics.median(v)) for x in v) for v in tms.values() if v)
            cand = (round(res, 3), level, how, round(statistics.median(tms['C']), 3) if tms['C'] else None, round(statistics.median(tms['H']), 4) if tms['H'] else None)
            allres[(level, how)] = round(res, 3)
            if best is None or cand[0] < best[0] - 1e-9: best = cand
    natoms = sum(len(t.get('C', [])) + len(t.get('H', [])) for t in d['isomers'].values())
    rows.append(dict(case=case, folder=folder, file=name, n_iso=len(d['isomers']), n_values=natoms, layout='header' if d['header'] else 'plain',
                     res=best[0] if best else None, level=best[1] if best else None, energies=best[2] if best else None,
                     tms_c=best[3] if best else None, tms_h=best[4] if best else None, probs=d['probs'],
                     residuals={'%s/%s' % k: v for k, v in allres.items()}))
json.dump(rows, open(os.path.join(work, 'dp4files_explained.json'), 'w'), indent=0)
ok = [r for r in rows if r['res'] is not None and r['res'] <= 0.0101]
print('.dp4 files: %d | calc column reproduced within 0.01 ppm by some (shieldings, energies) choice: %d | values in those: %d' % (len(rows), len(ok), sum(r['n_values'] for r in ok)))
print('how they were made (shieldings level, energies from) -> files:')
for k, n in collections.Counter((r['folder'], r['level'], r['energies']) for r in ok).most_common():
    print('   folder %-5s shieldings %-5s energies %-7s %3d' % (k + (n,)))
print('reference constants implied (C, H) -> files:', collections.Counter((r['tms_c'], r['tms_h']) for r in ok).most_common(12))
print('not reproduced: %d' % (len(rows) - len(ok)))
for r in rows:
    if r not in ok:
        print('   %-6s %-5s %-16s iso %2d | best: %s' % (r['case'], r['folder'], r['file'], r['n_iso'], ('%s/%s residual %.2f' % (r['level'], r['energies'], r['res'])) if r['level'] else 'labels do not fit any level'))

# which energies: count only the files where the choice of energies makes a visible difference
dft = [r for r in ok if r['level'] in ('Opt', 'Opt_E')]
tab = collections.Counter()
for r in dft:
    fits = sorted(h for h in ('nmr', 'opt', 'e') if r['residuals'].get('%s/%s' % (r['level'], h), 9) <= 0.0101)
    tab['+'.join(fits)] += 1
print('reproduced outputs that use DFT-geometry shieldings: %d | conformer energies that reproduce them: %s' % (len(dft), dict(tab)))
mmcopies = [r for r in ok if r['folder'] != 'MM' and r['level'] == 'MM']
print('outputs stored in an Opt or Opt_E folder that hold MM-level results: %d' % len(mmcopies))
