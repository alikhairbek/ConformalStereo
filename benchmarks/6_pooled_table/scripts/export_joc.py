"""DP4+ log-likelihoods of every candidate of the 2015 benchmark (48 cases, 24 levels), recomputed from the printed
shieldings with the routines and inputs of the first archive (lib: dp4plus.py, joc_data.py; inputs: joc_*.csv).
Usage: export_joc.py <lib dir> <inputs dir> <out.csv>      (needs numpy and scipy)
"""
import csv
import json
import os
import re
import sys

sys.path.insert(0, sys.argv[1])
import numpy as np
import dp4plus as D
import joc_data as J


def main(inp, out):
    comp, levels = J.load_shieldings(os.path.join(inp, 'joc_shieldings.csv'))
    par = J.load_params(os.path.join(inp, 'joc_dp4plus_parameters.csv'))
    pub, order = J.load_published(os.path.join(inp, 'joc_published_probs.csv'))
    ovf = os.path.join(inp, 'joc_sp2_overrides.json')
    OV = {c: {tuple(k.split('|')): v for k, v in d.items()} for c, d in json.load(open(ovf)).items()} if os.path.exists(ovf) else {}
    w = csv.writer(open(out, 'w', newline=''))
    w.writerow(['case', 'family', 'level', 'candidates', 'candidate', 'is_correct', 'n_C', 'n_H', 'logL_scaled', 'logL_unscaled'])
    n = 0
    worst = 0.0
    close = 0
    far = set()
    for case, cands in order.items():
        rows = comp[case]
        flags = J.sp2_flags(rows, 100.0, OV.get(case))
        pairs = {c: J.pair(rows, comp[c]) for c in cands}
        keep = [i for i, r in enumerate(rows) if r['exp'] is not None and all(dict(pairs[c]).get(i) is not None for c in cands)]
        nuc = np.array([rows[i]['nuc'] for i in keep])
        dexp = np.array([rows[i]['exp'] for i in keep])
        sp2 = np.array([flags[i] for i in keep])
        for lvl in levels:
            terms = []
            for c in cands:
                pm = dict(pairs[c])
                sh = np.array([comp[c][pm[i]]['sh'][lvl] for i in keep])
                terms.append(D.log_terms(sh, dexp, nuc, sp2, par[lvl], truncate_nu=True))
            p = D.probabilities(terms, nuc)[('DP4+', 'All')]
            q = pub[case][lvl].get(('DP4+', 'All'))
            if q:
                d = max(abs(100 * v - q[c]) for c, v in zip(cands, p))
                worst = max(worst, d)
                close += d <= 0.5
                far.add(case) if d > 10 else None
            for c, t in zip(cands, terms):
                w.writerow([case, re.match(r'\d+', case).group(), lvl, len(cands), c, int(c == case), int((nuc == 'C').sum()), int((nuc == 'H').sum()),
                            '%.6f' % t['s'].sum(), '%.6f' % t['u'].sum()])
            n += 1
    print('records %d (cases %d, levels %d); every candidate within 0.5 point of the published DP4+ probability in %d records; '
          'largest difference %.1f points; cases with a difference above 10 points at some level: %s'
          % (n, len(order), len(levels), close, worst, ', '.join(sorted(far)) or 'none'))


if __name__ == '__main__':
    main(*sys.argv[2:4])
