"""DP4-AI deposit: one record per candidate isomer and level with Boltzmann-averaged shieldings.

Energies follow the original PyDP4 code (Gaussian.ReadEnergies): last SCF energy of the NMR job at level MM,
of the optimisation at level Opt, of the M06-2X single point at level Opt_E. T = 298.15 K.
Usage: build_dp4ai.py <res/core> <work_dir>"""
import sys, os, re, gzip, json, glob, math, collections
core, work = sys.argv[1], sys.argv[2]
R, T, HARTREE = 8.3145, 298.15, 2625.499629554010
P1 = re.compile(r'^(?P<zip>[^!]+)\.zip!/(?P<top>[^/]+)/(?P<level>MM|Opt|Opt_E)/(?:(?P<sub>[^/]+)/)?(?P<stem>[^/]+?)_(?P<iso>\d+)ginp(?P<conf>\d+)\.(?:out|log)$')
nmr = collections.defaultdict(dict); eopt = {}; esp = {}; skipped = collections.Counter(); meta = {}
for f in sorted(glob.glob(os.path.join(core, 'gaussian_dp4ai_*.jsonl.gz'))):
    for line in gzip.open(f, 'rt'):
        r = json.loads(line)
        m = P1.match(r['path'])
        if not m:
            skipped[re.sub(r'\d+', 'N', r['path'].split('!/', 1)[1])[:60]] += 1
            continue
        d = m.groupdict(); key = (d['zip'], d['level'], d['stem'], int(d['iso']), int(d['conf']))
        if 'iso' in r:
            nmr[(d['zip'], d['level'], d['sub'] or 'flat')][key[2:]] = r
        elif r.get('scf_m') == 'RM062X':
            esp[key] = r
        elif r.get('scf_m') == 'RB3LYP':
            eopt[key] = r
        else:
            skipped['unclassified job: %s %s' % (d['level'], r.get('scf_m'))] += 1
print('NMR sets (zip, level, folder):', len(nmr), '| optimisation jobs:', len(eopt), '| M06-2X jobs:', len(esp), '| skipped:', dict(skipped))
out = open(os.path.join(work, 'dp4ai_isomers.jsonl'), 'w')
jobs = open(os.path.join(work, 'dp4ai_conformers.csv'), 'w')
jobs.write('case,level,nmr_set,stem,isomer,conformer,energy_hartree,energy_from,rel_kJ,population,nmr_ended_normally,energy_job_ended_normally,opt_converged,cpu_s_nmr,cpu_s_opt,cpu_s_e\n')
summary = collections.OrderedDict(); problems = collections.Counter()
for (zp, level, sub), confs in sorted(nmr.items()):
    by_iso = collections.defaultdict(list)
    for (stem, iso, conf), r in confs.items():
        by_iso[(stem, iso)].append((conf, r))
    for (stem, iso), lst in sorted(by_iso.items()):
        lst.sort()
        rows = []
        for conf, r in lst:
            k = (zp, level, stem, iso, conf)
            if level == 'MM':
                e, src, ejob = r.get('scf'), 'nmr', r
            elif level == 'Opt':
                ejob = eopt.get(k); e, src = (ejob or {}).get('scf'), 'opt'
            else:
                ejob = esp.get(k); e, src = (ejob or {}).get('scf'), 'e'
            if e is None:
                problems['%s: conformer without an energy' % level] += 1
                continue
            if any(v is None for v in r['iso']):
                problems['shielding overflow'] += 1
                continue
            rows.append((conf, r, e, src, ejob))
        if not rows:
            problems['isomer without usable conformers'] += 1
            continue
        els = [r['el'] for _, r, _, _, _ in rows]
        if len(set(els)) != 1:
            problems['conformers disagree on the atom list'] += 1
            continue
        emin = min(e for _, _, e, _, _ in rows)
        rel = [(e - emin) * HARTREE for _, _, e, _, _ in rows]
        w = [math.exp(-x * 1000 / (R * T)) for x in rel]
        q = sum(w); pop = [x / q for x in w]
        n_at = len(rows[0][1]['iso'])
        sigma = [sum(p * r['iso'][a] for p, (_, r, _, _, _) in zip(pop, rows)) for a in range(n_at)]
        r0 = rows[0][1]
        basis = 'def2svp' if 'def2svp' in r0.get('route', '').lower() else ('6-311g(d)' if '6-311g(d)' in r0.get('route', '').lower() else r0.get('route', ''))
        rec = dict(case=zp, level=level, nmr_set=sub, stem=stem, isomer=iso, n_conf=len(lst), n_used=len(rows), elements=r0['el'].split(),
                   sigma=[round(x, 4) for x in sigma], basis=basis, solvent=(r0.get('solv') or '').lower(), nbf=r0.get('nbf'),
                   max_pop=round(max(pop), 4), n_eff=round(1 / sum(p * p for p in pop), 2), energy_from=rows[0][3],
                   conf_sigma=[[c] + r['iso'] for c, r, _, _, _ in rows], conf_pop=[round(p, 5) for p in pop])
        out.write(json.dumps(rec, separators=(',', ':')) + '\n')
        for (conf, r, e, src, ejob), x, p in zip(rows, rel, pop):
            ko = eopt.get((zp, level, stem, iso, conf)) or {}; ke = esp.get((zp, level, stem, iso, conf)) or {}
            jobs.write('%s,%s,%s,%s,%d,%d,%.9f,%s,%.4f,%.6f,%d,%d,%s,%s,%s,%s\n' % (zp, level, sub, stem, iso, conf, e, src, x, p, 1 if r.get('norm') else 0,
                       1 if ejob.get('norm') else 0, ko.get('opt_done', ''), r.get('cpu_s', ''), ko.get('cpu_s', ''), ke.get('cpu_s', '')))
        s = summary.setdefault((zp, level, sub), dict(iso=0, conf=0, used=0, nat=n_at, basis=basis, solv=rec['solvent']))
        s['iso'] += 1; s['conf'] += len(lst); s['used'] += len(rows)
out.close(); jobs.close()
print('problems:', dict(problems))
cases = sorted(set(k[0] for k in summary), key=lambda s: (re.sub(r'\d.*', '', s), int(re.search(r'\d+', s).group()), s))
print('%d cases | isomers x conformers per level (NMR set "nmr" unless noted)' % len(cases))
print('%-6s %4s %-9s %-10s | %-12s %-12s %-12s | other NMR sets' % ('case', 'atoms', 'basis', 'solvent', 'MM', 'Opt', 'Opt_E'))
tot = collections.Counter()
for c in cases:
    cells = []
    for lvl in ('MM', 'Opt', 'Opt_E'):
        s = summary.get((c, lvl, 'nmr')) or summary.get((c, lvl, 'flat'))
        cells.append('%2d x %5d%s' % (s['iso'], s['used'], '' if s['used'] == s['conf'] else '*') if s else '-')
        if s: tot[lvl + ' isomers'] += s['iso']; tot[lvl + ' conformers'] += s['used']
    other = ['%s/%s:%dx%d' % (k[1], k[2], v['iso'], v['used']) for k, v in summary.items() if k[0] == c and k[2] not in ('nmr', 'flat')]
    any_s = next(v for k, v in summary.items() if k[0] == c)
    print('%-6s %4d %-9s %-10s | %-12s %-12s %-12s | %s' % (c, any_s['nat'], any_s['basis'], any_s['solv'][:10], cells[0], cells[1], cells[2], ' '.join(other)))
print('totals:', dict(tot))
