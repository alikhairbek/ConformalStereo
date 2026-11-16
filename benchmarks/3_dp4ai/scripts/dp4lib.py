"""Shared helpers: TMS constants, shifts, reading the original .dp4 outputs and NMR descriptions of the deposits."""
import re, json, math, zipfile, gzip, collections
TMS_DEFAULT = (191.69255, 31.7518583)           # PyDP4 default (B3LYP/6-31G**), used when no entry matches
def tms_table(path):
    t = {}
    for line in open(path):
        b = line.split()
        if len(b) == 5:
            t[(b[0].lower(), b[1].lower(), b[2].lower())] = (float(b[3]), float(b[4]))
    return t
def to_shift(sigma, el, tms):
    c, h = tms
    return (c - sigma) / (1 - c / 1e6) if el == 'C' else (h - sigma) / (1 - h / 1e6)
SEC = re.compile(r'^Assigned (C|H) shifts(?: for isomer (\d+))?:', re.M)
ROW = re.compile(r'^([A-Z][a-z]?\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$', re.M)
def parse_dp4(text):
    """Original PyDP4 / DP4-AI output -> dict(isomers={n: {'C': [(label, calc, scaled, exp, err)], 'H': [...]}}, probs={'H': [...], 'C': [...], 'all': [...]})."""
    out = dict(isomers=collections.OrderedDict(), probs={}, header={})
    for k, v in re.findall(r'^(Solvent|Force Field|DFT NMR Functional|DFT NMR Basis|Number of isomers|Stats model|DFT Opt Functional|DFT Opt Basis|DFT Single Point Functional|DFT Single Point Basis)\s*=\s*(.+)$', text, re.M):
        out['header'][k] = v.strip()
    marks = [(m.start(), m.group(1), m.group(2)) for m in SEC.finditer(text)]
    iso_marks = [(m.start(), int(m.group(1))) for m in re.finditer(r'^Assigned shifts for isomer (\d+):', text, re.M)]
    for k, (pos, nuc, iso) in enumerate(marks):
        end = marks[k + 1][0] if k + 1 < len(marks) else len(text)
        stop = text.find('Results of DP4', pos)
        if 0 <= stop < end: end = stop
        if iso is None:                             # older layout: "Assigned shifts for isomer N:" precedes both tables
            prev = [n for p, n in iso_marks if p < pos]
            iso = prev[-1] if prev else 1
        rows = [(a, float(b), float(c), float(d), float(e)) for a, b, c, d, e in ROW.findall(text[pos:end])]
        out['isomers'].setdefault(int(iso), {})[nuc] = rows
    for titles, key in ((('Results of DP4 using Proton', 'Results of DP4 using proton data only'), 'H'),
                        (('Results of DP4 using Carbon', 'Results of DP4 using carbon data only'), 'C'),
                        (('Results of DP4', 'Results of DP4 using all available data'), 'all')):
        for title in titles:
            m = re.search(r'^%s:\s*\n((?:Isomer \d+:\s*[\d.]+%%\s*\n?)+)' % re.escape(title), text, re.M)
            if m:
                out['probs'][key] = [float(x) for x in re.findall(r'Isomer \d+:\s*([\d.]+)%', m.group(1))]
    out['raw'] = {}
    for nuc in 'CH':
        v = re.findall(r'^DP4 based on %s:\s*([-+0-9.eE]+|nan)' % nuc, text, re.M)
        if v: out['raw'][nuc] = [float(x) for x in v]
    head = text[:marks[0][0]] if marks else text
    cut = head.find('Assigned shifts for isomer')
    out['desc_text'] = head[:cut] if cut >= 0 else ('' if out['header'] else head)
    return out
def parse_desc(text):
    """NMR description (first line carbons, third line protons; labels in brackets, 'any' = unassigned)."""
    lines = text.replace('\r', '').split('\n')
    def one(line):
        vals = []
        for item in re.sub(r'\(.*?\)', '', line.strip()).split(','):
            item = item.strip()
            if item:
                try: vals.append(float(item))
                except ValueError: vals.append(None)
        return vals
    c = one(lines[0]) if lines else []
    h = one(lines[2]) if len(lines) > 2 else []
    rest = [l for l in lines[4:] if l.strip()]
    return dict(C=c, H=h, omits=[x for l in rest if 'OMIT' in l for x in l[5:].split(',')], equivalents=[l.split(',') for l in rest if 'OMIT' not in l], bad=sum(v is None for v in c + h))
class Small:
    """Small text files copied from a deposit (deduplicated): read by original path."""
    def __init__(self, core, key):
        self.index = {}
        for l in gzip.open('%s/small_index_%s.tsv.gz' % (core, key), 'rt', encoding='utf-8').read().split('\n')[1:]:
            if l:
                path, size, md5, stored, same = l.split('\t')
                self.index[path] = same or path
        import glob
        self.zips = [zipfile.ZipFile(p) for p in sorted(glob.glob('%s/small_%s_*.zip' % (core, key)))]
        self.where = {n: z for z in self.zips for n in z.namelist()}
    def paths(self): return list(self.index)
    def read(self, path):
        real = self.index[path].replace('!/', '/')
        return self.where[real].read(real).decode('latin-1')
