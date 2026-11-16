"""Inputs of a DP4+App calculation that the SI files do not print, recovered from the probabilities they do print.

DP4+App reads, next to the experimental shifts, two columns set by the user: 'sp2' (which nuclei take the sp2
parameters) and 'exchange' (pairs of nuclei whose assignment may be swapped; for each candidate the larger experimental
shift then goes with the smaller shielding). Neither is in the SI. Starting from the rule used elsewhere (carbons at or
above 100 ppm are sp2, protons follow their carbon, no exchange), single changes are tried one at a time, a class flip
of one nucleus or one more exchangeable pair, and the change that brings the recomputed probabilities closest to the
printed ones is kept, until no change helps. Candidate pairs: two rows of the same nucleus that are protons on the same
carbon (by label) or lie close in experimental shift (1H within 0.6 ppm, 13C within 3 ppm); a nucleus enters one pair
at most; pairs with equal experimental values are left out because exchanging them changes nothing. Class changes are
limited to what chemistry allows: carbons between 85 and 125 ppm with their protons, and protons between 4.3 and
7.5 ppm that have no carbon of their own in the table.
The result is accepted only if every printed probability is then reproduced within the tolerance."""
import re, itertools
import numpy as np
import mmdp4 as M

def candidate_pairs(rows, h_window=0.6, c_window=3.0):
    base = lambda lab: (re.match(r'^(\d+)', lab.strip()) or [None, None])[1] if re.match(r'^(\d+)', lab.strip()) else None
    out = []
    for i, j in itertools.combinations(range(len(rows)), 2):
        a, b = rows[i], rows[j]
        if a['nuc'] != b['nuc'] or a['exp'] == b['exp']: continue
        d = abs(a['exp'] - b['exp'])
        geminal = a['nuc'] == 'H' and base(a['label']) is not None and base(a['label']) == base(b['label'])
        if geminal or d <= (h_window if a['nuc'] == 'H' else c_window): out.append((i, j))
    return out

def evaluate(D, rows, shieldings, sp2, pairs, targets, pars):
    """shieldings: {level: [array per candidate]}; targets: [(level, key, candidate index, value)].
    Returns (sum of absolute deviations, largest, {level: probabilities dict})."""
    exp = [r['exp'] for r in rows]; nuc = [r['nuc'] for r in rows]; cache = {}; tot = 0.0; worst = 0.0
    for lev, key, k, val in targets:
        if lev not in cache:
            sh = [M.exchange(s, exp, pairs) for s in shieldings[lev]] if pairs else shieldings[lev]
            cache[lev] = M.probabilities(D, sh, exp, nuc, sp2, pars[lev])
        pr = cache[lev]
        if key not in pr: continue
        d = abs(100 * float(pr[key][k]) - val); tot += d; worst = max(worst, d)
    return tot, worst, cache

def flip_moves(rows):
    """Class changes that chemistry allows: a carbon between 85 and 125 ppm (acetal or alkene?), together with the
    protons that carry its label; a proton between 4.3 and 7.5 ppm that has no carbon of its own in the table."""
    base = lambda lab: (re.match(r'^(\d+)', lab.strip()).group(1) if re.match(r'^(\d+)', lab.strip()) else lab.strip())
    clabels = {base(r['label']) for r in rows if r['nuc'] == 'C'}
    moves = []
    for i, r in enumerate(rows):
        if r['nuc'] == 'C' and 85 <= r['exp'] <= 125:
            moves.append([i] + [j for j, q in enumerate(rows) if q['nuc'] == 'H' and base(q['label']) == base(r['label'])])
        if r['nuc'] == 'H' and 4.3 <= r['exp'] <= 7.5 and base(r['label']) not in clabels: moves.append([i])
    return moves

def infer(D, rows, shieldings, sp2_start, targets, pars, max_steps=8, min_gain=0.5, tol=1.0):
    """Greedy search on all targets. Returns dict(sp2, pairs, flips, start=(total, worst), end=(total, worst), steps)."""
    sp2 = list(sp2_start); pairs = []
    cur, worst, _ = evaluate(D, rows, shieldings, sp2, pairs, targets, pars); start = (cur, worst)
    cand = candidate_pairs(rows); fm = flip_moves(rows); flips = []; steps = 0
    while steps < max_steps and worst > tol:
        best = None; used = {x for p in pairs for x in p}
        for mv in fm:
            s2 = list(sp2)
            for i in mv: s2[i] = not sp2[mv[0]]
            v = evaluate(D, rows, shieldings, s2, pairs, targets, pars)[0]
            if best is None or v < best[0]: best = (v, 'flip', mv)
        for p in cand:
            if p in pairs or p[0] in used or p[1] in used: continue
            v = evaluate(D, rows, shieldings, sp2, pairs + [p], targets, pars)[0]
            if best is None or v < best[0]: best = (v, 'pair', p)
        for p in pairs:
            v = evaluate(D, rows, shieldings, sp2, [q for q in pairs if q != p], targets, pars)[0]
            if best is None or v < best[0]: best = (v, 'unpair', p)
        if best is None or cur - best[0] < min_gain: break
        if best[1] == 'flip':
            new = not sp2[best[2][0]]
            for i in best[2]: sp2[i] = new
            flips.append(best[2][0])
        elif best[1] == 'pair': pairs.append(best[2])
        else: pairs.remove(best[2])
        cur, worst, _ = evaluate(D, rows, shieldings, sp2, pairs, targets, pars); steps += 1
    net = sorted({i for i in flips if flips.count(i) % 2 == 1})
    return dict(sp2=sp2, pairs=pairs, flips=net, start=start, end=(cur, worst), steps=steps)
