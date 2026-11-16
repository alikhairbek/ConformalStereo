"""Loads the DP4-AI deposit records written by the Kaggle notebook: shieldings and energies per conformer."""
import os, re, glob, gzip, json, math, collections
R, T, HARTREE = 8.3145, 298.15, 2625.499629554010          # constants of the original code
P1 = re.compile(r'^(?P<zip>[^!]+)\.zip!/(?P<top>[^/]+)/(?P<level>MM|Opt|Opt_E)/(?:(?P<sub>[^/]+)/)?(?P<stem>[^/]+?)_(?P<iso>\d+)ginp(?P<conf>\d+)\.(?:out|log)$')
class Deposit:
    def __init__(self, core, key='dp4ai'):
        self.E = collections.defaultdict(dict)      # (case, level, kind) -> {(iso, conf): hartree}
        self.SIG = collections.defaultdict(dict)    # (case, level) -> {(iso, conf): (elements, sigma)}
        self.meta = {}
        for f in sorted(glob.glob(os.path.join(core, 'gaussian_%s_*.jsonl.gz' % key))):
            for line in gzip.open(f, 'rt'):
                r = json.loads(line); m = P1.match(r['path'])
                if not m: continue
                d = m.groupdict(); k = (int(d['iso']), int(d['conf']))
                if d['sub'] not in (None, 'nmr', 'opt', 'e'): continue          # extra folders (NL1, NL2, man) are not part of a level
                if 'iso' in r:
                    if any(v is None for v in r['iso']): continue
                    self.SIG[(d['zip'], d['level'])][k] = (r['el'].split(), r['iso']); kind = 'nmr'
                    self.meta.setdefault((d['zip'], d['level']), dict(route=r.get('route', ''), solv=(r.get('solv') or '').lower(), nbf=r.get('nbf')))
                elif r.get('scf_m') == 'RM062X': kind = 'e'
                elif r.get('scf_m') == 'RB3LYP': kind = 'opt'
                else: continue
                if r.get('scf') is not None: self.E[(d['zip'], d['level'], kind)][k] = r['scf']
    def cases(self):
        return sorted(set(c for c, _ in self.SIG), key=lambda s: (re.sub(r'\d.*', '', s), int(re.search(r'\d+', s).group()), s))
    def isomers(self, case, level):
        return sorted(set(i for i, _ in self.SIG.get((case, level), {})))
    def averaged(self, case, level, iso, how, temperature=T):
        """Boltzmann-averaged shieldings of one isomer. how = 'nmr' | 'opt' | 'e' (source of the conformer energies) or 'uniform'.
        Returns (elements, sigma, populations{conf: p}) or None."""
        sig = self.SIG.get((case, level), {})
        confs = sorted(c for (i, c) in sig if i == iso)
        if not confs: return None
        if how == 'uniform': en = {c: 0.0 for c in confs}
        else:
            src = self.E.get((case, level, how)) or {}
            en = {c: src[(iso, c)] for c in confs if (iso, c) in src}
        if not en: return None
        use = sorted(en); e0 = min(en.values())
        w = [math.exp(-(en[c] - e0) * HARTREE * 1000 / (R * temperature)) for c in use]
        q = sum(w); els = sig[(iso, use[0])][0]
        if any(sig[(iso, c)][0] != els for c in use): return None
        return els, [sum(wi / q * sig[(iso, c)][1][a] for wi, c in zip(w, use)) for a in range(len(els))], {c: wi / q for wi, c in zip(w, use)}
