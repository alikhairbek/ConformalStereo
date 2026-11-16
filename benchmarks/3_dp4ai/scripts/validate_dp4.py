"""Full check of my DP4 against the original output files of the deposit.
For every .dp4 file whose shifts I reproduce: (1) is the printed pairing of atoms and peaks the sorted pairing?
(2) do my probabilities (my shifts, my pairing, experimental values at full precision where the file echoes them) equal the printed ones?
Usage: validate_dp4.py <res/core> <work_dir>      (run explain_dp4files.py and desc_inventory.py first)"""
import sys, os, re, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dp4lib as L, dp4core as D, dp4ai_data as A
core, work = sys.argv[1:3]
dep = A.Deposit(core); S = L.Small(core, 'dp4ai')
expl = {(r['case'], r['folder'], r['file']): r for r in json.load(open(os.path.join(work, 'dp4files_explained.json')))}
DESCS = json.load(open(os.path.join(work, 'dp4ai_descriptions.json')))
stat = collections.Counter(); worst = collections.defaultdict(float); rows = []; notes = []
for p in sorted(x for x in S.paths() if x.endswith('.dp4')):
    m = re.match(r'^([^!]+)\.zip!/[^/]+/(MM|Opt|Opt_E)/([^/]+)$', p); case, folder, name = m.groups()
    ex = expl[(case, folder, name)]
    if ex['res'] is None or ex['res'] > 0.0101: stat['shifts not reproduced'] += 1; continue
    d = L.parse_dp4(S.read(p))
    isos = []
    for k in d['isomers']:
        els, sig, _ = dep.averaged(case, ex['level'], k, ex['energies'])
        tms = (ex['tms_c'] or 188.0, ex['tms_h'] or 32.0)
        one = {}
        for nuc in 'CH':
            idx = [a for a, e in enumerate(els) if e == nuc]
            one[nuc] = (['%s%d' % (nuc, a + 1) for a in idx], [L.to_shift(sig[a], nuc, tms) for a in idx])
        isos.append(one)
    desc = L.parse_desc(d['desc_text']) if d['desc_text'].strip() else None
    first = next(iter(d['isomers'].values()))
    if desc is None or desc['bad']:
        exp = {nuc: [r[3] for r in first.get(nuc, [])] for nuc in 'CH'}; stat['experimental values taken from the tables (2 decimals)'] += 1
    else:
        exp = {'C': desc['C'], 'H': desc['H']}; stat['experimental values from the echoed description'] += 1
        if desc['omits'] or desc['equivalents']: stat['description has OMIT or equivalence lines'] += 1
    rule = '2020' if d['header'] else '2019'                      # the two output layouts come from two generations of the code
    if desc is None or desc['bad']:                               # header layout: look for the stored description with the same values
        want = {nuc: sorted(round(v, 2) for v in exp[nuc]) for nuc in 'CH'}
        for cand in DESCS.get(case, []):
            if all(len(cand[nuc]) == len(want[nuc]) and all(abs(a - b) < 0.0051 for a, b in zip(sorted(cand[nuc]), want[nuc])) for nuc in 'CH' if want[nuc]):
                exp = {nuc: (cand[nuc] if want[nuc] else []) for nuc in 'CH'}; desc = dict(bad=0, omits=[], equivalents=[]); stat['... of which matched to a stored description (full precision)'] += 1
                break
    out = D.dp4(isos, exp['C'], exp['H'], rule=rule)
    # (1) pairing
    same_pairs = True; n_rows = 0
    for nuc in 'CH':
        if nuc + '_detail' not in out: continue
        for (k, tabs), det in zip(d['isomers'].items(), out[nuc + '_detail']):
            printed = {r[0]: r[3] for r in tabs.get(nuc, [])}; mine = {l: e for l, c, s, e in det}
            n_rows += len(printed)
            if set(printed) != set(mine) or any(abs(printed[l] - mine[l]) > 0.0051 for l in printed): same_pairs = False
    stat['pairing identical' if same_pairs else 'pairing differs'] += 1
    # (2) probabilities
    dev = {}
    for key in ('C', 'H', 'all'):
        if key in out and key in d['probs'] and len(out[key]) == len(d['probs'][key]):
            dev[key] = max(abs(100 * a - b) for a, b in zip(out[key], d['probs'][key]))
    full = desc is not None and not desc['bad']
    for key, v in dev.items():
        if same_pairs and full: worst[key] = max(worst[key], v)
    ok = same_pairs and all(v <= 0.11 for v in dev.values()) and dev
    stat['probabilities reproduced within 0.1 point' if ok else 'probabilities differ'] += 1
    if not ok: notes.append('%s/%s/%s pairs=%s full=%s dev=%s' % (case, folder, name, same_pairs, full, {k: round(v, 1) for k, v in dev.items()}))
    rows.append(dict(case=case, folder=folder, file=name, level=ex['level'], energies=ex['energies'], same_pairs=same_pairs, full_precision=full, dev=dev, n_values=n_rows, ok=bool(ok)))
json.dump(rows, open(os.path.join(work, 'dp4files_validated.json'), 'w'), indent=0)
print('files:'); [print('   %-70s %d' % kv) for kv in stat.items()]
print('largest difference to a printed probability where pairing agrees and experimental values have full precision:', {k: round(v, 3) for k, v in worst.items()})
print('atom-peak pairs compared in files that agree: %d' % sum(r['n_values'] for r in rows if r['ok']))
print('files that differ (%d):' % len(notes)); print('\n'.join('   ' + n for n in notes[:40]))

# What the authors' own outputs say, taken as printed (no recomputation): how often is isomer 1 ranked first?
seen = set(); t = collections.Counter()
for p in sorted(x for x in S.paths() if x.endswith('.dp4')):
    m = re.match(r'^([^!]+)\.zip!/[^/]+/(MM|Opt|Opt_E)/([^/]+)$', p); case, folder, name = m.groups()
    ex = expl[(case, folder, name)]
    if ex['res'] is None or ex['res'] > 0.0101: continue
    d = L.parse_dp4(S.read(p))
    desc = L.parse_desc(d['desc_text']) if d['desc_text'].strip() else None
    if desc is not None and not desc['bad'] and len(desc['H']) >= 3 and sorted(desc['H']) == [float(i) for i in range(1, len(desc['H']) + 1)]:
        t['outputs made with placeholder proton data (left out)'] += 1; continue
    pr = d['probs'].get('all') or d['probs'].get('H') or d['probs'].get('C')
    if not pr: continue
    key = (case, bool(d['header']), tuple(pr))
    if key in seen: continue
    seen.add(key)
    top = max(range(len(pr)), key=lambda i: pr[i])
    t['distinct outputs'] += 1; t['isomer 1 first'] += top == 0; t['calls at 99% or more'] += pr[top] >= 99; t['... of which isomer 1'] += (pr[top] >= 99 and top == 0)
print('\nstored outputs as printed, reproduced ones with real experimental data:', dict(t))
