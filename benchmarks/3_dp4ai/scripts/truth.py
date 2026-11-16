"""Which candidate isomer is the structure the authors gave as input (file <case>_.sdf, the reported structure)?
Compares stereo InChIs computed from the 3D structures; enantiomers count as the same candidate (NMR cannot tell them apart).
Usage: truth.py <res/core> <work_dir>      (needs RDKit)"""
import sys, os, re, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dp4lib as L
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
core, work = sys.argv[1:3]
S = L.Small(core, 'dp4ai'); ps = S.paths()
def inchi_of(molblock):
    m = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=True)
    if m is None: return None
    Chem.AssignStereochemistryFrom3D(m)
    return Chem.MolToInchi(m)
def key(inchi):                                   # drop the enantiomer flag: /m0 and /m1 are mirror images
    return re.sub(r'/m[01]', '', inchi) if inchi else None
cases = sorted(set(p.split('.zip!/')[0] for p in ps if '.zip!/' in p), key=lambda s: (re.sub(r'\d.*', '', s), int(re.search(r'\d+', s).group()), s))
out = {}; flags = collections.Counter()
for case in cases:
    res = {}
    for folder in ('MM', 'Opt', 'Opt_E'):
        inp = [p for p in ps if re.match(r'^%s\.zip!/[^/]+/%s/[^/]+_\.sdf$' % (re.escape(case), folder), p)]
        iso = {int(re.search(r'_(\d+)\.sdf$', p).group(1)): p for p in ps if re.match(r'^%s\.zip!/[^/]+/%s/[^/]+_\d+\.sdf$' % (re.escape(case), folder), p)}
        if not inp or not iso: continue
        t = S.read(inp[0]); k0 = key(inchi_of(t)); title = t.split('\n')[0].strip()
        ks = {n: key(inchi_of(S.read(p))) for n, p in iso.items()}
        dep = {}
        for n, p in iso.items():
            q = p[:-4] + '.inchi'
            if q in S.index: dep[n] = S.read(q).strip()
        hits = [n for n, k in ks.items() if k is not None and k == k0]
        res[folder] = dict(title=title, n_iso=len(iso), hits=hits, distinct=len(set(ks.values())), unreadable=sum(k is None for k in ks.values()),
                           dep_distinct=len(set(re.sub(r'/m[01]', '', v) for v in dep.values())), inchi_true=k0)
    one = sorted(p for p in ps if re.match(r'^%s\.zip!/[^/]+/(MM|Opt|Opt_E)/[^/]+_1\.sdf$' % re.escape(case), p))
    i1 = inchi_of(S.read(one[0])) if one else None
    m1 = Chem.MolFromMolBlock(S.read(one[0]), removeHs=False) if one else None
    if m1 is not None: Chem.AssignStereochemistryFrom3D(m1)
    n_sym = len(set(r for a, r in zip(m1.GetAtoms(), Chem.CanonicalRankAtoms(m1, breakTies=False, includeChirality=True)) if a.GetSymbol() == 'C')) if m1 is not None else None
    const = i1.split('/t')[0].split('/b')[0] if i1 else None
    res['_summary'] = dict(inchi_isomer1=i1, constitution=const, carbon_classes_isomer1=n_sym, n_carbons=sum(a.GetSymbol() == 'C' for a in m1.GetAtoms()) if m1 is not None else None)
    hits_all = set(tuple(v['hits']) for f, v in res.items() if f != '_summary')
    if len(hits_all) == 1 and len(next(iter(hits_all))) == 1:
        res['_summary'].update(true_isomer=next(iter(hits_all))[0], true_from='input structure file')
    else:                                           # no input file: can the number of 13C signals decide?
        cands = {int(re.search(r'_(\d+)\.sdf$', p).group(1)): p for p in ps if re.match(r'^%s\.zip!/[^/]+/MM/[^/]+_\d+\.sdf$' % re.escape(case), p)}
        classes = {}
        for n, p in cands.items():
            m = Chem.MolFromMolBlock(S.read(p), removeHs=False); Chem.AssignStereochemistryFrom3D(m)
            classes[n] = len(set(r for a, r in zip(m.GetAtoms(), Chem.CanonicalRankAtoms(m, breakTies=False, includeChirality=True)) if a.GetSymbol() == 'C'))
        dpath = [p for p in ps if re.match(r'^%s\.zip!/[^/]+/MM/%sNMR$' % (re.escape(case), re.escape(case)), p)]
        n_sig = len(set(L.parse_desc(S.read(dpath[0]))['C'])) if dpath else None
        fit = [n for n, c in classes.items() if c >= (n_sig or 10 ** 9)]
        note = 'no input structure file; candidates have %s distinct carbons by symmetry, the spectrum lists %s distinct 13C shifts' % (classes, n_sig)
        if len(fit) == 1:
            res['_summary'].update(true_isomer=fit[0], true_from='symmetry: ' + note)
        else:
            res['_summary'].update(true_isomer=None, true_from='undetermined: ' + note)
        print('       %s: %s -> reported structure = isomer %s' % (case, note, res['_summary']['true_isomer']))
    out[case] = res
    hs = set(tuple(v['hits']) for f, v in res.items() if f != '_summary')
    flag = 'ok' if len(hs) == 1 and len(next(iter(hs))) == 1 else 'CHECK'
    flags[flag] += 1
    first = next((v for f, v in res.items() if f != '_summary'), {})
    print('%-6s %-5s true isomer per folder: %s | candidates %s distinct %s (deposit inchi distinct %s) | input title: %s' % (
        case, flag, ' '.join('%s=%s' % (f, ','.join(map(str, v['hits'])) or '-') for f, v in res.items() if f != '_summary'), first.get('n_iso'), first.get('distinct'), first.get('dep_distinct'), first.get('title', '')[:28]))
print('input structure file identifies the reported structure:', dict(flags))
print('reported structure is isomer 1 in %d of %d cases' % (sum(r['_summary'].get('true_isomer') == 1 for r in out.values()), len(out)))
fam = collections.defaultdict(list)
for c, r in out.items(): fam[r['_summary']['constitution']].append(c)
shared = [v for v in fam.values() if len(v) > 1]
print('constitutions: %d for %d cases | cases sharing one candidate set: %s' % (len(fam), len(out), shared))
for c, r in out.items(): r['_summary']['family'] = min(fam[r['_summary']['constitution']])
json.dump(out, open(os.path.join(work, 'dp4ai_truth.json'), 'w'), indent=0)
