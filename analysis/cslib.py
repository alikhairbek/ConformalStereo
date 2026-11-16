"""Shared routines of the final-method analysis: records, scores, family-weighted conformal thresholds, sets, summaries.

A record is one experimental data set assigned at one workflow (benchmark, method, kind of geometry, NMR level): its
candidates with their log-likelihoods relative to the best one, and which candidate is the actual structure.
A family is a set of records that must not be split between calibration and test (records that share candidates, or
the same molecule in two benchmarks).
"""
import csv
import collections
import math

import numpy as np

INF = float('inf')


def group_of(method, geo):
    if method != 'DP4':
        return '%s, %s geometries' % (method, geo)
    if geo == 'MMFF':
        return 'DP4, MMFF geometries'
    if geo in ('PM6', 'PM7'):
        return 'DP4, semi-empirical geometries'
    return 'DP4, HF or DFT geometries'


def kclass(K):
    return '2' if K == 2 else '3 or 4' if K <= 4 else '5 to 8' if K <= 8 else '9 or more'


KCLASSES = ['2', '3 or 4', '5 to 8', '9 or more']


def softmax(ell):
    w = np.exp(ell - ell.max())
    return w / w.sum()


def lse(ell):
    m = ell.max()
    return m + math.log(np.exp(ell - m).sum())


def tail_mass_score(ell, beta=1.0):
    """minus the log of the probability mass of the candidates that are no more likely than the candidate, at inverse
    temperature beta: 0 for the best candidate, -log p for the least likely one (adaptive prediction sets on a log scale)"""
    x = beta * ell
    x = x - lse(x)
    order = np.argsort(x)                                  # ascending log-probability
    xs = x[order]
    m = xs.max()
    cum = m + np.log(np.cumsum(np.exp(xs - m)))            # log of the mass of the candidates up to and including each
    out = np.empty_like(x)
    out[order] = -cum
    out[out < 1e-12] = 0.0                                 # the best candidate carries all the mass: exactly zero
    return out


SCORES = ['likelihood ratio per shift', 'likelihood ratio / sqrt(N)', 'likelihood ratio', 'stated probability', 'cumulative probability', 'rank']


def scores_of(ell, n):
    """every score of the declared family for the candidates of one record (larger = less plausible)"""
    gap = ell.max() - ell
    rank = np.empty(len(ell))
    rank[np.argsort(-ell, kind='stable')] = np.arange(1, len(ell) + 1)
    return {'likelihood ratio per shift': gap / n,           # log-likelihood below the best candidate, divided by the number of shifts
            'likelihood ratio / sqrt(N)': gap / math.sqrt(n),
            'likelihood ratio': gap,
            'stated probability': lse(ell) - ell,            # minus the log of the probability the method states
            'cumulative probability': tail_mass_score(ell),
            'rank': rank}


def load(table, column='log_ratio'):
    """-> list of records (dicts); records without a value in `column` are left out"""
    raw = collections.OrderedDict()
    for r in csv.DictReader(open(table)):
        if r[column] == '':
            continue
        raw.setdefault((r['benchmark'], r['method'], r['geometry'], r['level'], r['case']), []).append(r)
    recs = []
    for (bench, method, geo, level, case), rr in raw.items():
        ell = np.array([float(r[column]) for r in rr])
        truth = [int(r['is_true']) for r in rr]
        assert sum(truth) == 1, (bench, level, case)
        k = truth.index(1)
        top = int(np.argmax(ell))
        n = int(rr[0]['n_nuclei'])
        p = softmax(ell)
        sc = scores_of(ell, n)
        rec = dict(bench=bench, method=method, geo=geo, level=level, case=case, wf=(bench, method, geo, level), group=group_of(method, geo),
                   family=rr[0]['family'], clean=int(rr[0]['clean']), K=len(rr), N=n, ell=ell, k=k, top=top, first=(k == top),
                   p=p, ptop=float(p[top]), pa=float(p[k]), names=[r['candidate'] for r in rr],
                   t={s: float(v[k]) for s, v in sc.items()}, s={s: np.sort(v) for s, v in sc.items()}, tops={s: float(v[top]) for s, v in sc.items()})
        recs.append(rec)
    return recs


def threshold(items, target):
    """family-weighted conformal threshold. items: [(score of the actual structure, family)].
    A score weighs one over the number of records of its family; every family weighs one over the number of families
    plus one. Infinite if the families are too few for the target."""
    if not items:
        return INF
    cnt = collections.Counter(f for _, f in items)
    G = len(cnt)
    acc = 0.0
    for s, f in sorted(items, key=lambda x: x[0]):
        acc += 1.0 / cnt[f] / (G + 1)
        if acc >= target - 1e-12:
            return s
    return INF


def record_threshold(scores, target):
    """threshold that treats the records as the exchangeable units (split conformal prediction)"""
    n = len(scores)
    k = math.ceil((n + 1) * target - 1e-9)
    return sorted(scores)[k - 1] if 1 <= k <= n else INF


def apply_set(rec, score, thr):
    """-> (covered, size) of {candidates with score <= thr}, the top candidate always included"""
    size = int(np.searchsorted(rec['s'][score], thr, side='right'))
    cov = rec['t'][score] <= thr
    if rec['tops'][score] > thr:
        size += 1
        cov = cov or rec['first']
    return bool(cov), size


def stated_set(p, k, target):
    """smallest set of candidates whose probabilities add up to the target -> (covered, size)"""
    order = np.argsort(-p, kind='stable')
    cum = np.cumsum(p[order])
    m = int(np.searchsorted(cum, target - 1e-12, side='left')) + 1
    m = min(m, len(p))
    return bool(k in order[:m]), m


class Fit:
    """temperature of one pool of records: the factor by which the log-likelihoods have to be divided so that the
    probabilities of the actual structures are as large as possible (family-weighted maximum likelihood)"""

    def __init__(self, recs):
        self.recs = recs
        R, Kmax = len(recs), max(r['K'] for r in recs)
        self.L = np.full((R, Kmax), -np.inf)
        for i, r in enumerate(recs):
            self.L[i, :r['K']] = r['ell']
        self.la = np.array([r['ell'][r['k']] for r in recs])
        self.fam = np.array([r['family'] for r in recs])

    def weights(self, mask):
        cnt = collections.Counter(self.fam[mask])
        return np.array([1.0 / cnt[f] if m else 0.0 for f, m in zip(self.fam, mask)])

    def nll(self, beta, w):
        with np.errstate(invalid='ignore'):
            x = np.where(np.isfinite(self.L), beta * self.L, -np.inf)
        return float((w * (np.log(np.exp(x).sum(axis=1)) - beta * self.la)).sum())

    def temperature(self, mask=None):
        from scipy.optimize import minimize_scalar
        if mask is None:
            mask = np.ones(len(self.recs), bool)
        w = self.weights(mask)
        res = minimize_scalar(lambda u: self.nll(math.exp(u), w), bounds=(math.log(1e-3), math.log(50.0)), method='bounded', options={'xatol': 1e-4})
        return 1.0 / math.exp(res.x)


def unit_means(recs, values):
    """values: one tuple per record -> {(workflow, family): mean tuple}"""
    by = collections.defaultdict(list)
    for r, v in zip(recs, values):
        by[(r['wf'], r['family'])].append(v)
    return {u: tuple(sum(x[i] for x in vs) / len(vs) for i in range(len(vs[0]))) for u, vs in by.items()}


def summarize(recs, cov, size, boot=0, seed=0):
    """coverage, set size and the kinds of sets, averaged per family within a workflow and then over these units.
    With boot > 0, 95% intervals from resampling families (all workflows of a family together)."""
    if not recs:
        return None
    vals = [(float(c), float(s), float(s == 1), float(s == 1 and r['first']), float(s == r['K'] and r['K'] > 1), s / r['K'], float(r['K']), float(r['first']))
            for r, c, s in zip(recs, cov, size)]
    um = unit_means(recs, vals)
    U = len(um)
    m = [sum(v[i] for v in um.values()) / U for i in range(8)]
    out = dict(coverage=m[0], size=m[1], single=m[2], single_right=(m[3] / m[2] if m[2] > 0 else float('nan')), every=m[4], share=m[5], K=m[6], top=m[7],
               units=U, families=len({u[1] for u in um}), records=len(recs))
    if boot:
        fam = collections.defaultdict(lambda: np.zeros(4))
        for (wf, f), v in um.items():
            fam[f] += np.array([v[0], v[1], v[2], 1.0])
        A = np.array(list(fam.values()))
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, len(A), size=(boot, len(A)))
        S = A[idx].sum(axis=1)
        for j, name in enumerate(('coverage', 'size', 'single')):
            q = np.quantile(S[:, j] / S[:, 3], [0.025, 0.975])
            out[name + '_lo'], out[name + '_hi'] = float(q[0]), float(q[1])
    return out


def pct(x, d=1):
    return '%.*f%%' % (d, 100 * x)
