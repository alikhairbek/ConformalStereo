"""Inventory of the experimental NMR descriptions kept in the DP4-AI deposit (text files and the copies echoed in old outputs).
Usage: desc_inventory.py <res/core> <work_dir>"""
import sys, os, re, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dp4lib as L
core, work = sys.argv[1:3]
S = L.Small(core, 'dp4ai')
def looks_like_desc(t):
    first = t.lstrip().split('\n')[0]
    return bool(re.match(r'^\s*-?\d+(\.\d+)?\s*\(', first))
def placeholder(vals):
    return len(vals) >= 3 and all(abs(v - round(v)) < 1e-9 for v in vals) and sorted(vals) == [float(i) for i in range(1, len(vals) + 1)]
by_case = collections.defaultdict(lambda: collections.defaultdict(list)); bad = collections.Counter()
for p in sorted(S.paths()):
    m = re.match(r'^([^!]+)\.zip!/[^/]+/(MM|Opt|Opt_E)/([^/]+)$', p)
    if not m: continue
    case, folder, name = m.groups()
    if name.endswith('.dp4'):
        t = L.parse_dp4(S.read(p))['desc_text']; src = 'echo'
        if not t.strip(): continue
    elif re.search(r'NMR', name) and '.' not in name:
        t = S.read(p); src = 'file'
        if not looks_like_desc(t): bad['not a description: ' + re.sub(r'\d+', 'N', name)] += 1; continue
    else: continue
    try: d = L.parse_desc(t)
    except Exception as e: bad['unreadable'] += 1; continue
    if d['bad']: bad['malformed number in %s' % src] += 1; continue
    key = (tuple(sorted(d['C'])), tuple(sorted(d['H'])), tuple(sorted(map(str, d['omits']))), tuple(sorted(str(e) for e in d['equivalents'])))
    by_case[case][key].append((folder, name, src))
out = {}
print('skipped:', dict(bad))
print('%-6s %s' % ('case', 'distinct usable descriptions: nC/nH [flags] x places (example names)'))
for case in sorted(by_case, key=lambda s: (re.sub(r'\d.*', '', s), int(re.search(r'\d+', s).group()), s)):
    cells = []; out[case] = []
    for key, places in sorted(by_case[case].items(), key=lambda kv: -len(kv[1])):
        c, h, om, eq = key
        flags = ('Hplaceholder ' if placeholder(h) else '') + ('Cplaceholder ' if len(c) < 2 else '') + ('omit ' if om else '') + ('equiv ' if eq else '') + ('labels ' if False else '')
        names = sorted(set(n for _, n, s in places if s == 'file'))
        cells.append('%d/%d %sx%d (%s)' % (len(c), len(h), '[' + flags.strip() + '] ' if flags else '', len(places), ','.join(names)[:38]))
        out[case].append(dict(C=list(c), H=list(h), omits=list(om), equivalents=list(eq), places=places, h_placeholder=placeholder(h), c_placeholder=len(c) < 2))
    print('%-6s %s' % (case, ' | '.join(cells)))
json.dump(out, open(os.path.join(work, 'dp4ai_descriptions.json'), 'w'))
