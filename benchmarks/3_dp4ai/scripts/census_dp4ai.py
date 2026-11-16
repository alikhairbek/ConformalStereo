"""Census of the Gaussian job records of the DP4-AI deposit: folder patterns, job types, methods, anomalies.
Usage: census_dp4ai.py <res/core>"""
import sys, os, re, gzip, json, glob, collections
core = sys.argv[1]
pat = re.compile(r'^(?P<zip>[^!]+)\.zip!/(?P<top>[^/]+)/(?P<level>[^/]+)/(?P<sub>[^/]+)/(?P<mol>.+?)_(?P<iso>\d+)ginp(?P<conf>\d+)(?P<suf>[^.]*)\.(?P<ext>out|log)$')
shape = collections.Counter(); odd = collections.Counter(); combos = collections.Counter(); notnorm = collections.Counter(); multi = 0
routes = collections.Counter(); n = 0; sufs = collections.Counter(); solv = collections.Counter(); vers = collections.Counter()
for f in sorted(glob.glob(os.path.join(core, 'gaussian_dp4ai_*.jsonl.gz'))):
    for line in gzip.open(f, 'rt'):
        r = json.loads(line); n += 1
        m = pat.match(r['path'])
        if not m:
            odd[re.sub(r'\d+', 'N', r['path'].split('!/', 1)[1])[:80]] += 1
            continue
        d = m.groupdict()
        if r.get('njobs', 1) > 1: multi += 1
        cls = 'nmr' if 'iso' in r else ('opt' if r.get('n_geom', 1) > 1 or 'opt' in r.get('route', '').lower() else 'sp')
        combos[(d['level'], d['sub'], cls, r.get('scf_m', '-'))] += 1
        routes[(d['level'], d['sub'], re.sub(r'\s+', ' ', r.get('route', '-').lower()))] += 1
        sufs[d['suf']] += 1
        solv[r.get('solv', 'none')] += 1
        vers[r.get('ver', '-')] += 1
        if not r.get('norm'): notnorm[(d['level'], d['sub'], 'err' if r.get('err') else 'cut', r.get('tail', '')[:40])] += 1
        if d['zip'] != d['top'] or d['top'] != d['mol']: shape['zip/top/mol differ: %s %s %s' % (d['zip'], d['top'], d['mol'])] += 1
print('job-step records:', n, '| not matching the expected path pattern:', sum(odd.values()), '| records from multi-step files:', multi)
for k, v in odd.most_common(12): print('   odd path', v, k)
for k, v in shape.most_common(8): print('  ', v, k)
print('level / folder / job / SCF method:')
for k, v in sorted(combos.items()): print('   %-7s %-4s %-4s %-14s %6d' % (k + (v,)))
print('file-name suffixes after the conformer number:', dict(sufs))
print('solvents:', dict(solv)); print('versions:', dict(vers))
print('routes (most common per folder):')
seen = collections.Counter()
for (lvl, sub, rt), v in routes.most_common():
    if seen[(lvl, sub)] < 3:
        seen[(lvl, sub)] += 1
        print('   %-6s %-4s %6d  %s' % (lvl, sub, v, rt[:150]))
print('not ended normally:', sum(notnorm.values()))
for k, v in notnorm.most_common(10): print('  ', v, k)
