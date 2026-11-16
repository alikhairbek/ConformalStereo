"""DP4 for every molecule of the DP4-AI deposit at every protocol, recomputed from the Gaussian outputs.

Protocols (geometry | conformer energies; NMR is mPW1PW91 throughout):
  MM            MMFF geometries              | energy of the NMR job                (folder MM)
  DFT_EB3LYP    B3LYP/6-31G(d,p) geometries  | B3LYP energy of the optimisation     (folder Opt, the code's rule for a workflow with 'o')
  DFT_EmPW      same geometries              | energy of the NMR job                (folder Opt; what most stored outputs used)
  DFT_EM062X    same geometries              | M06-2X/def2-TZVP single point        (folder Opt_E, the code's rule for a workflow with 'e')
Experimental shifts: the description file <case>NMR of the deposit (the version found in most level folders);
pairing and statistics as in the 2020 code (dp4core). The reported structure of each case comes from truth.py
(isomer 1 in all 46 cases).
Flags per case and protocol: qc_clean = every conformer that carries more than 1% of the population has its shieldings (and its
energy) on the intended geometry (qc_geometry.py); exact_counts = the description lists exactly one value per carbon and per proton,
so the sorted pairing needs no further rule.
Usage: dp4ai_run.py <res/core> <work_dir> <out_dir> [NMR|NMRold]"""
import sys, os, re, json, csv, collections, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dp4lib as L, dp4core as D, dp4ai_data as A
core, work, out = sys.argv[1:4]; which = sys.argv[4] if len(sys.argv) > 4 else 'NMR'
PROTOCOLS = [('MM', 'MM', 'nmr'), ('DFT_EB3LYP', 'Opt', 'opt'), ('DFT_EmPW', 'Opt', 'nmr'), ('DFT_EM062X', 'Opt_E', 'e')]
TMS = (188.452125, 32.1243166667)          # any constant gives the same DP4: the shifts are rescaled against experiment
dep = A.Deposit(core); S = L.Small(core, 'dp4ai'); ps = S.paths()
truth = json.load(open(os.path.join(work, 'dp4ai_truth.json')))
QC = json.load(open(os.path.join(work, 'qc_flags.json')))
def description(case):
    """Most frequent usable version of <case><which> across the three level folders; returns (dict, note)."""
    seen = collections.OrderedDict()
    for folder in ('MM', 'Opt', 'Opt_E'):
        for p in ps:
            if re.match(r'^%s\.zip!/[^/]+/%s/%s%s$' % (re.escape(case), folder, re.escape(case), which), p):
                try: d = L.parse_desc(S.read(p))
                except Exception: continue
                if d['bad'] or len(d['H']) < 2: continue
                key = (tuple(d['C']), tuple(d['H']))
                seen.setdefault(key, []).append((folder, d))
    if not seen: return None, 'no usable description'
    key, lst = max(seen.items(), key=lambda kv: (len(kv[1]), kv[1][0][0] == 'MM'))
    note = '' if len(seen) == 1 else 'versions differ between level folders (%s)' % '; '.join('/'.join(f for f, _ in v) for v in seen.values())
    return lst[0][1], note
cases = dep.cases(); long_rows = []; case_rows = []; notes = collections.OrderedDict()
for case in cases:
    desc, note = description(case)
    if desc is None: notes[case] = note; continue
    if note: notes[case] = note
    fam = truth[case]['_summary']['family']
    n_cand = max(len(dep.isomers(case, f)) for f in ('MM', 'Opt', 'Opt_E'))
    for prot, folder, how in PROTOCOLS:
        isos = dep.isomers(case, folder)
        TRUE = truth[case]['_summary']['true_isomer']
        if TRUE not in isos: notes['%s %s' % (case, prot)] = 'no NMR calculations for the reported structure'; continue
        built = []; conf = []; ok = True; bad_pop = 0.0
        for k in isos:
            av = dep.averaged(case, folder, k, how)
            if av is None: ok = False; break
            els, sig, pop = av
            if folder != 'MM':
                need = ('geom_ok', 'e_ok') if how == 'e' else ('geom_ok',)
                bad_pop = max(bad_pop, sum(p for c, p in pop.items() if not all(QC.get('%s|%s|%d|%d' % (case, folder, k, c), {}).get(n, False) for n in need)))
            one = {}
            for nuc in 'CH':
                idx = [a for a, e in enumerate(els) if e == nuc]
                one[nuc] = (['%s%d' % (nuc, a + 1) for a in idx], [L.to_shift(sig[a], nuc, TMS) for a in idx])
            built.append(one); conf.append((len(pop), max(pop.values()), 1 / sum(p * p for p in pop.values())))
        if not ok: notes['%s %s' % (case, prot)] = 'an isomer has no conformer with both shieldings and energy'; continue
        if len(isos) < n_cand: notes['%s %s' % (case, prot)] = 'only %d of %d candidates have NMR calculations' % (len(isos), n_cand)
        r = D.dp4(built, desc['C'], desc['H'], rule='2020')
        n_c_atoms, n_h_atoms = len(built[0]['C'][0]), len(built[0]['H'][0])
        exact = (len(desc['C']) in (n_c_atoms,) or len(desc['C']) < 2) and len(desc['H']) == n_h_atoms
        for data in ('all', 'C', 'H'):
            if data not in r: continue
            pr = r[data]; order = sorted(range(len(pr)), key=lambda i: -pr[i]); t = isos.index(TRUE)
            case_rows.append(dict(case=case, family=fam, protocol=prot, data=data, n_candidates=len(isos), n_candidates_full=n_cand,
                                  n_C=r.get('C_n', 0) if data != 'H' else 0, n_H=r.get('H_n', 0) if data != 'C' else 0,
                                  p_true=round(pr[t], 6), rank_true=order.index(t) + 1, top_isomer=isos[order[0]], p_top=round(pr[order[0]], 6),
                                  p_best_wrong=round(max(p for i, p in enumerate(pr) if i != t), 6), top1_correct=int(order[0] == t),
                                  conformers=sum(c[0] for c in conf), conformers_true=conf[t][0], n_eff_true=round(conf[t][2], 2),
                                  unused_C_peaks=r.get('C_unused_peaks', 0), unused_H_peaks=r.get('H_unused_peaks', 0),
                                  exp_C=len(desc['C']), exp_H=len(desc['H']), atoms_C=n_c_atoms, atoms_H=n_h_atoms,
                                  exact_counts=int(exact), qc_clean=int(bad_pop <= 0.01), qc_bad_population=round(bad_pop, 4)))
            for k, p in zip(isos, pr):
                long_rows.append(dict(case=case, family=fam, protocol=prot, data=data, isomer=k, is_true=int(k == TRUE), probability='%.6g' % p,
                                      exact_counts=int(exact), qc_clean=int(bad_pop <= 0.01)))
tag = '' if which == 'NMR' else '_' + which
for name, rows in (('dp4ai_cases%s.csv' % tag, case_rows), ('dp4ai_probabilities%s.csv' % tag, long_rows)):
    w = csv.DictWriter(open(os.path.join(out, name), 'w', newline=''), fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print('description file: <case>%s | cases with results: %d of %d' % (which, len(set(r['case'] for r in case_rows)), len(cases)))
for k, v in notes.items(): print('   note %-18s %s' % (k, v))
flagged = sorted(set('%s %s (%.0f%% of the population)' % (r['case'], r['protocol'], 100 * r['qc_bad_population']) for r in case_rows if not r['qc_clean']))
print('not clean at a DFT protocol (%d): %s' % (len(flagged), '; '.join(flagged)))
inexact = sorted(set('%s (C %d/%d, H %d/%d)' % (r['case'], r['exp_C'], r['atoms_C'], r['exp_H'], r['atoms_H']) for r in case_rows if not r['exact_counts']))
print('peaks listed != atoms (%d cases): %s' % (len(inexact), '; '.join(inexact)))
def table(title, keep):
    print('\n%s' % title)
    print('%-11s %-4s | cases | top-1 correct | median P(true) | top >= 99%%: right/all | true <= 1%%' % ('protocol', 'data'))
    for prot, _, _ in PROTOCOLS:
        for data in ('all', 'C', 'H'):
            rs = [r for r in case_rows if r['protocol'] == prot and r['data'] == data and keep(r)]
            if not rs: continue
            pt = sorted(r['p_true'] for r in rs); conf = [r for r in rs if r['p_top'] >= 0.99]
            print('%-11s %-4s | %5d | %3d (%4.1f%%)   | %14.3f | %10d/%-9d | %10d' % (prot, data, len(rs), sum(r['top1_correct'] for r in rs), 100 * sum(r['top1_correct'] for r in rs) / len(rs),
                  pt[len(pt) // 2], sum(r['top1_correct'] for r in conf), len(conf), sum(r['p_true'] <= 0.01 for r in rs)))
table('all cases, as deposited', lambda r: True)
table('clean cases only (geometry check passed, one experimental value per atom)', lambda r: r['qc_clean'] and r['exact_counts'])
