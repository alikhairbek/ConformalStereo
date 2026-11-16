"""Pairing of an unassigned list of experimental shifts with the calculated shifts of one candidate.

With as many values as atoms the two sorted lists are paired in order; that is also the pairing with the highest DP4
likelihood. With fewer values than atoms (signals not observed) or more (extra signals), the 2020 code of PyDP4 pairs
from the top and leaves the atoms, or the signals, with the lowest shifts out. Here the order is kept and the items left
out are those whose omission gives the highest DP4 likelihood: a dynamic programme over the two sorted lists, repeated
with the internal scaling until the pairing no longer changes.
"""
import math

LOG2 = math.log(2.0)


def term(error, sigma):
    return math.log(max(math.erfc(abs(error) / sigma / math.sqrt(2)), 1e-300))


def linear(calc, exp):
    n = len(calc)
    mx, my = sum(exp) / n, sum(calc) / n
    sxx = sum((a - mx) ** 2 for a in exp)
    slope = sum((a - mx) * (b - my) for a, b in zip(exp, calc)) / sxx
    return slope, my - slope * mx


def align(calc, exp, sigma):
    """calc, exp: ascending lists on the same scale -> pairs [(i, j)] of the order-preserving pairing of
    min(len) items with the highest summed likelihood term"""
    n, m = len(calc), len(exp)
    if n == m:
        return [(i, i) for i in range(n)]
    swap = n < m
    a, b = (exp, calc) if swap else (calc, exp)       # a is the longer list; every item of b is paired
    na, nb = len(a), len(b)
    neg = float('-inf')
    f = [[neg] * (nb + 1) for _ in range(na + 1)]
    for i in range(na + 1):
        f[i][0] = 0.0
    for i in range(1, na + 1):
        for j in range(1, min(i, nb) + 1):
            skip = f[i - 1][j]
            take = f[i - 1][j - 1] + term(a[i - 1] - b[j - 1], sigma)
            f[i][j] = take if take >= skip else skip
    pairs = []
    i, j = na, nb
    while j > 0:
        if f[i][j] == f[i - 1][j - 1] + term(a[i - 1] - b[j - 1], sigma) and not (f[i - 1][j] > f[i][j] - 1e-15 and f[i - 1][j] == f[i][j] and i - 1 >= j):
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        else:
            i -= 1
    pairs.reverse()
    return [(q, p) for p, q in pairs] if swap else pairs


def optimal(calc, exp, sigma, rounds=10):
    """calc: calculated shifts of the atoms (any order, referenced so that they are close to the experimental scale),
    exp: experimental values (any order).
    -> (indices of the atoms kept, indices of the experimental values kept), each in ascending order of shift.
    The alternation of pairing and scaling can stop at a poor pairing, so it is started three times (calculated shifts
    as they are; pairing from the top, as the 2020 code; pairing from the bottom) and the best end point is kept."""
    ci = sorted(range(len(calc)), key=lambda i: calc[i])
    ei = sorted(range(len(exp)), key=lambda i: exp[i])
    c = [calc[i] for i in ci]
    e = [exp[i] for i in ei]
    k = min(len(c), len(e))
    if len(c) == len(e):
        return ci, ei
    starts = [align(c, e, sigma),
              [(len(c) - k + t, len(e) - k + t) for t in range(k)],
              [(t, t) for t in range(k)]]
    best = None
    for pairs in starts:
        value = float('-inf')
        for _ in range(rounds):
            slope, icpt = linear([c[i] for i, j in pairs], [e[j] for i, j in pairs])
            if slope <= 0:
                break
            scaled = [(x - icpt) / slope for x in c]
            value = sum(term(scaled[i] - e[j], sigma) for i, j in pairs)
            new = align(scaled, e, sigma)
            if new == pairs:
                break
            pairs = new
        else:
            slope, icpt = linear([c[i] for i, j in pairs], [e[j] for i, j in pairs])
            value = sum(term((c[i] - icpt) / slope - e[j], sigma) for i, j in pairs)
        if best is None or value > best[0]:
            best = (value, pairs)
    pairs = best[1]
    return [ci[i] for i, j in pairs], [ei[j] for i, j in pairs]
