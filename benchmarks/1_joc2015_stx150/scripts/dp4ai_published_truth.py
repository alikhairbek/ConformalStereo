"""Cross-check of the reported structure: the isomer marked as correct in the published spreadsheet against
isomer 1 of the deposit (identified from the input structure files, see truth.py of the benchmark package).
The spreadsheet names isomers by the signs of the InChI stereo layer; the deposit keeps one .inchi file per isomer.
Usage: dp4ai_published_truth.py <out_dir with dp4ai_published_raw.csv> <res/core of run 1> <scripts dir of the benchmark package>"""
import sys, os, re, csv, collections
out, core, pkg = sys.argv[1:4]
sys.path.insert(0, pkg)
import dp4lib as L
rows = list(csv.DictReader(open(os.path.join(out, 'dp4ai_published_raw.csv'))))
by = collections.defaultdict(lambda: collections.OrderedDict())
for r in rows: by[r['level']].setdefault(r['molecule'], []).append(r)
truth = {m: v[0]['stereo'] for m, v in by['MM'].items()}                 # MM sheet: the marked isomer is the first of each block
chk = collections.Counter()
for m, v in by['Opt_E'].items():
    mk = [x['stereo'] for x in v if x['marked'] == '1']
    chk['Opt_E mark = MM mark' if mk == [truth[m]] else 'Opt_E mark differs'] += 1
for m, v in by['Opt'].items():
    chk['Opt sheet: first row is the MM-marked isomer' if v[0]['stereo'] == truth[m] else 'Opt sheet: first row differs'] += 1
    mk = [i for i, x in enumerate(v) if x['marked'] == '1']
    chk['Opt sheet: mark sits one row below' if mk == [1] else 'Opt sheet: mark elsewhere (%s)' % mk] += 1
print(dict(chk))
S = L.Small(core, 'dp4ai'); ps = S.paths()
agree = collections.Counter(); notes = []; mapping = {}
for m in truth:
    iso = {}
    for p in ps:
        mm = re.match(r'^%s\.zip!/[^/]+/(MM|Opt|Opt_E)/[^/]+_(\d+)\.inchi$' % re.escape(m), p)
        if mm and int(mm.group(2)) not in iso:
            t = S.read(p).strip(); layer = re.search(r'/t([^/]+)', t)
            iso[int(mm.group(2))] = ''.join(re.findall(r'[+\-?]', layer.group(1))) if layer else ''
    if not iso: agree['not in the deposit'] += 1; notes.append('%s: not in the deposit' % m); continue
    hit = [n for n, s in iso.items() if s == truth[m]]
    mapping[m] = {s: n for n, s in iso.items()}
    if hit == [1]: agree['published correct isomer = isomer 1 of the deposit'] += 1
    else: agree['different'] += 1; notes.append('%s: published %s, deposit isomers %s' % (m, truth[m], iso))
    if set(iso.values()) != set(x['stereo'] for x in by['MM'][m]): notes.append('%s: stereo labels differ between sheet %s and deposit %s' % (m, sorted(x['stereo'] for x in by['MM'][m])[:4], sorted(iso.values())[:4]))
print(dict(agree)); print('\n'.join('   ' + n[:220] for n in notes[:12]))
w = csv.writer(open(os.path.join(out, 'dp4ai_published.csv'), 'w', newline=''))
w.writerow(['level', 'molecule', 'stereo', 'deposit_isomer', 'is_true', 'ai_all', 'ai_C', 'ai_H', 'pair_all', 'pair_C', 'pair_H'])
for r in rows:
    w.writerow([r['level'], r['molecule'], r['stereo'], mapping.get(r['molecule'], {}).get(r['stereo'], ''), int(r['stereo'] == truth[r['molecule']]),
                r['ai_all'], r['ai_C'], r['ai_H'], r['pair_all'], r['pair_C'], r['pair_H']])
