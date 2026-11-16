"""Does the overconfidence grow with the size of the molecule?
Usage: temperature_size.py <pooled_candidates.csv> <out folder> [bootstrap=2000] [seed=20261114]

Model: the log-likelihoods of a record with N shifts are divided by a temperature T = T30 * (N / 30) ** b, with one
T30 per workflow and one exponent b for a group of workflows. b = 0: the same temperature for every molecule;
b = 1/2 or 1: the temperature grows with the square root of N or with N, as it would if the errors of the shifts of
one molecule were not independent. T30 and b maximise the family-weighted likelihood of the actual structures
(grid search); the interval of b comes from resampling families.
"""
import sys, os, csv, math, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cslib as cs

table, out = sys.argv[1:3]
BOOT = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 20261114
recs = cs.load(table)
BENCH = list(collections.OrderedDict((r['bench'], 1) for r in recs))
by_wf = collections.OrderedDict()
for r in recs:
    by_wf.setdefault(r['wf'], []).append(r)
B = np.round(np.arange(-1.0, 3.0001, 0.05), 2)                    # exponents
LT = np.linspace(math.log(0.05), math.log(1000.0), 160)           # log of T30
# NLL[workflow] = (array [b, T30, family of the workflow], its families): family-weighted negative log-likelihood of the actual structures
NLL = {}
for wf, rr in by_wf.items():
    R, Kmax = len(rr), max(r['K'] for r in rr)
    L = np.full((R, Kmax), -np.inf)
    for i, r in enumerate(rr):
        L[i, :r['K']] = r['ell']
    la = np.array([r['ell'][r['k']] for r in rr])
    n = np.array([r['N'] for r in rr], float)
    fl = sorted({r['family'] for r in rr})
    fp = {f: i for i, f in enumerate(fl)}
    cnt = collections.Counter(r['family'] for r in rr)
    w = np.array([1.0 / cnt[r['family']] for r in rr])
    fi = np.array([fp[r['family']] for r in rr])
    A = np.zeros((len(B), len(LT), len(fl)))
    for ib, b in enumerate(B):
        scale = (n / 30.0) ** b                                     # T = T30 * scale
        for it, lt in enumerate(LT):
            beta = 1.0 / (math.exp(lt) * scale)
            x = beta[:, None] * L
            np.add.at(A[ib, it], fi, w * (np.log(np.exp(x).sum(axis=1)) - beta * la))
    NLL[wf] = (A, fl)


def estimate(wfs, fset=None, M=None):
    """profile over T30 per workflow, then the exponent.
    Without M: -> (b, total nll at every b). With M (resamples x families of fset, multiplicities): -> b of every resample."""
    if M is None:
        tot = np.zeros(len(B))
        for wf in wfs:
            tot += NLL[wf][0].sum(axis=2).min(axis=1)
        return float(B[int(np.argmin(tot))]), tot
    col = {f: i for i, f in enumerate(fset)}
    tot = np.zeros((len(B), M.shape[0]))
    for wf in wfs:
        A, fl = NLL[wf]
        Mw = M[:, [col[f] for f in fl]]                             # resamples x families of the workflow
        tot += np.tensordot(A, Mw, axes=([2], [1])).min(axis=1)     # [b, resample]
    return B[np.argmin(tot, axis=0)]


rng = np.random.default_rng(SEED)
w = csv.writer(open(os.path.join(out, 'final_temperature_size.csv'), 'w', newline=''))
w.writerow(['workflows', 'number', 'families', 'exponent', 'lo', 'hi', 'share_of_resamples_above_zero', 'gain_in_log_likelihood_over_b0', 'gain_of_b_half', 'gain_of_b_one'])
print('exponent b of T = T30 * (N/30)^b: estimate, 95%% interval from resampling families (%d resamples), log-likelihood gained over b = 0 (family-weighted)' % BOOT)
groups = [(b, [wf for wf in by_wf if wf[0] == b]) for b in BENCH]
groups += [(g, [wf for wf in by_wf if cs.group_of(wf[1], wf[2]) == g]) for g in collections.OrderedDict((r['group'], 1) for r in recs)]
groups += [('all workflows', list(by_wf))]
i0, ih, i1 = int(np.argmin(abs(B))), int(np.argmin(abs(B - 0.5))), int(np.argmin(abs(B - 1.0)))
for name, wfs in groups:
    fset = sorted({r['family'] for wf in wfs for r in by_wf[wf]})
    b, tot = estimate(wfs)
    M = rng.multinomial(len(fset), np.full(len(fset), 1.0 / len(fset)), size=BOOT).astype(float)
    bs = estimate(wfs, fset, M)
    lo, hi = np.quantile(bs, [0.025, 0.975])
    w.writerow([name, len(wfs), len(fset), '%.8f' % b, '%.8f' % lo, '%.8f' % hi, '%.8f' % (bs > 0).mean(), '%.8f' % (tot[i0] - tot.min()), '%.8f' % (tot[i0] - tot[ih]), '%.8f' % (tot[i0] - tot[i1])])
    print('   %-32s workflows %2d families %3d | b = %+.2f (%+.2f to %+.2f), above zero in %5.1f%% | gain over b=0: best %6.2f, b=1/2 %6.2f, b=1 %6.2f' % (
        name, len(wfs), len(fset), b, lo, hi, 100 * (bs > 0).mean(), tot[i0] - tot.min(), tot[i0] - tot[ih], tot[i0] - tot[i1]))
