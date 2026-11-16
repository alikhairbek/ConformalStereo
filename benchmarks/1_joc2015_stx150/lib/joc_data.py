"""Shared loaders for the DP4+ (JOC 2015) tables."""
import csv, re, collections
import numpy as np

def load_shieldings(path):
    comp = collections.OrderedDict()
    rd = csv.DictReader(open(path))
    levels = [c for c in rd.fieldnames if c.count('/') == 2]
    for r in rd:
        d = comp.setdefault(r['compound'], [])
        d.append(dict(nuc=r['nucleus'], label=r['label'], exp=(float(r['exp']) if r['exp'] != '' else None),
                      sh={l: float(r[l]) for l in levels}))
    return comp, levels

def load_params(path):
    return {r['level']: {k: float(v) for k, v in r.items() if k != 'level'} for r in csv.DictReader(open(path))}

def load_published(path):
    pub = collections.OrderedDict()   # case -> level -> (method, data) -> {candidate: pct}
    order = collections.OrderedDict() # case -> candidate order
    for r in csv.DictReader(open(path)):
        pub.setdefault(r['case'], collections.OrderedDict()).setdefault(r['level'], collections.OrderedDict()) \
           .setdefault((r['method'], r['data']), collections.OrderedDict())[r['candidate']] = float(r['prob_pct'])
        o = order.setdefault(r['case'], [])
        if r['candidate'] not in o: o.append(r['candidate'])
    return pub, order

def base_label(label):
    m = re.match(r'^(\d+)', label.strip())
    return m.group(1) if m else label.strip()

def sp2_flags(rows, c_threshold=100.0, overrides=None):
    """Heuristic hybridisation flags from the experimental shifts: carbons above the threshold are sp2 (or sp);
    a proton takes the flag of the carbon carrying the same atom number."""
    overrides = overrides or {}
    cflag = {}
    for r in rows:
        if r['nuc'] == 'C' and r['exp'] is not None:
            cflag[r['label']] = overrides.get(('C', r['label']), r['exp'] >= c_threshold)
    cbase = {}
    for lab, f in cflag.items():
        cbase.setdefault(base_label(lab), f)
    flags = []
    for r in rows:
        if r['nuc'] == 'C':
            flags.append(cflag.get(r['label'], False))
        else:
            if ('H', r['label']) in overrides: flags.append(overrides[('H', r['label'])])
            elif r['label'] in cflag: flags.append(cflag[r['label']])
            elif base_label(r['label']) in cbase: flags.append(cbase[base_label(r['label'])])
            else: flags.append(r['exp'] is not None and r['exp'] >= 5.0)
    return flags

def pair(case_rows, cand_rows):
    """Rows of the case that carry an experimental shift, matched to the candidate's rows by (nucleus, label, occurrence)."""
    idx = collections.defaultdict(list)
    for j, r in enumerate(cand_rows):
        idx[(r['nuc'], r['label'])].append(j)
    used = collections.Counter(); out = []
    for i, r in enumerate(case_rows):
        if r['exp'] is None:
            continue
        k = (r['nuc'], r['label']); lst = idx.get(k, [])
        j = lst[used[k]] if used[k] < len(lst) else None
        used[k] += 1
        out.append((i, j))
    return out
