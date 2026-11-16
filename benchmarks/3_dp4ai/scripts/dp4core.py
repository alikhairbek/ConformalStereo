"""DP4 as implemented in the original PyDP4 / DP4-AI code (single-Gaussian error model, internal scaling)."""
import math
SIGMA_C, SIGMA_H = 2.269372270818724, 0.18731058105269952
def linregress(x, y):
    n = len(x); mx = sum(x) / n; my = sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x); sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    slope = sxy / sxx
    return slope, my - slope * mx
def scale(calc, exp):
    """Internal scaling: regress calculated on experimental shifts, return (calc - intercept) / slope."""
    slope, intercept = linregress(exp, calc)
    return [(c - intercept) / slope for c in calc]
def log_prob(errors, sigma):
    """Sum over atoms of log(2 * Phi(-|e| / sigma)); erfc keeps the far tail finite."""
    return sum(math.log(max(math.erfc(abs(e) / sigma / math.sqrt(2)), 1e-300)) for e in errors)
def normalise(logs):
    m = max(logs); w = [math.exp(l - m) for l in logs]; q = sum(w)
    return [x / q for x in w]
def pair_sorted(calc, exp, ascending=False):
    """Unassigned ('any') data: the i-th largest experimental shift goes to the atom with the i-th largest calculated shift
    (2020 code; with fewer peaks than atoms the atoms with the lowest shifts stay without a peak). ascending=True is the
    pre-2020 rule: pairing starts from the smallest shifts, so the atoms with the highest shifts stay without a peak.
    Returns a list parallel to calc with the experimental value or None."""
    sign = 1 if ascending else -1
    order = sorted(range(len(calc)), key=lambda i: sign * calc[i])
    out = [None] * len(calc)
    for i, e in zip(order, sorted(exp, reverse=not ascending)):
        out[i] = e
    return out
def dp4(isomers, exp_c, exp_h, rule='2020', scaled=True):
    """isomers[k] = {'C': (labels, shifts), 'H': (labels, shifts)} -> probabilities for C, H and both.
    rule '2020': descending pairing, atoms left without a peak in any isomer are dropped from all isomers (DP4.ProcessIsomers).
    rule '2019': ascending pairing, each isomer on its own (older DP4.AssignExpNMR). A nucleus with fewer than 2 peaks is skipped."""
    res = {}; logs = {}
    for nuc, exp, sigma in (('C', exp_c, SIGMA_C), ('H', exp_h, SIGMA_H)):
        if not exp or len(exp) < 2 or not all(iso[nuc][0] for iso in isomers):
            continue
        paired = [pair_sorted(iso[nuc][1], exp, ascending=(rule == '2019')) for iso in isomers]
        dropped = set(l for iso, pr in zip(isomers, paired) for l, e in zip(iso[nuc][0], pr) if e is None) if rule == '2020' else set()
        ll = []; detail = []
        for iso, pr in zip(isomers, paired):
            keep = [(l, c, e) for l, c, e in zip(iso[nuc][0], iso[nuc][1], pr) if l not in dropped and e is not None]
            calc = [c for _, c, _ in keep]; ex = [e for _, _, e in keep]
            sc = scale(calc, ex) if scaled and len(set(ex)) > 1 else calc
            err = [s - e for s, e in zip(sc, ex)]
            ll.append(log_prob(err, sigma)); detail.append([(l, c, s, e) for (l, c, e), s in zip(keep, sc)])
        logs[nuc] = ll; res[nuc] = normalise(ll); res[nuc + '_detail'] = detail; res[nuc + '_n'] = len(detail[0])
        res[nuc + '_unused_peaks'] = max(0, len(exp) - min(len(iso[nuc][1]) for iso in isomers))
    if logs:
        tot = [sum(logs[n][k] for n in logs) for k in range(len(isomers))]
        res['all'] = normalise(tot); res['logs'] = logs
    return res
