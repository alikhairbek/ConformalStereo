"""The statistical model behind the coupling term of the deposit outputs that the Gaussian of the PyDP4 source does not
reproduce. Usage: jmodel.py <jrepro_couplings.json> <out folder>

Each printed term is a product over the couplings of one candidate. The model assumed here is the general form of
PyDP4's 'k' model: term = product of f(|exp - (calc - intercept) / slope|), with one slope, one intercept and one
density f for all outputs. f is left free (a cubic spline for log f); slope and intercept are found by minimising the
residual of log(term). If the form is right the residual falls to the rounding of the printed terms (three digits).
The check that does not use an output's own terms: leave one output out, fit on the others, predict it."""
import sys, os, json, math, collections
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize
src, out = sys.argv[1:3]
B = json.load(open(src))
gauss = set(f for f in set(b['file'] for b in B) if all(abs(b['log_gauss'] - math.log(b['printed'])) < 0.006 for b in B if b['file'] == f and b['printed'] > 0))
use = [b for b in B if b['file'] not in gauss and b['all_as_printed'] and b['printed'] > 0]
files = sorted(set(b['file'] for b in use))
print('outputs reproduced by the Gaussian of the source: %d (%s)' % (len(gauss), ', '.join(sorted(gauss))))
print('other outputs: %d, with %d candidate blocks and %d couplings used (all couplings as printed)' % (len(files), len(use), sum(len(b['pairs']) for b in use)))
K = 3


def setup(blocks):
    return (np.array([math.log(b['printed']) for b in blocks]), [np.array([p[0] for p in b['pairs']]) for b in blocks],
            [np.array([p[1] for p in b['pairs']]) for b in blocks])


def errors(E, C, s, i):
    return [np.abs(e - (c - i) / s) for e, c in zip(E, C)]


def spline_fit(X, y, hi, nk):
    t = np.concatenate(([0.0] * K, np.linspace(0.0, hi, nk), [hi] * K))
    A = np.vstack([BSpline.design_matrix(np.clip(x, 0.0, hi), t, K).toarray().sum(0) for x in X])
    sol = np.linalg.solve(A.T @ A + 1e-8 * np.eye(A.shape[1]), A.T @ y)
    return t, sol, A @ sol - y


def fit(blocks, nk=30, start=(1.05, -0.65)):
    y, E, C = setup(blocks)

    def loss(p):
        X = errors(E, C, p[0], p[1])
        hi = max(x.max() for x in X) + 1e-6
        return float(np.sqrt(np.mean(spline_fit(X, y, hi, nk)[2] ** 2)))
    best = None
    for s0, i0 in (start, (1.16, -0.37), (1.0, 0.0)):
        r = minimize(loss, [s0, i0], method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-9, maxiter=3000))
        if best is None or r.fun < best.fun:
            best = r
    s, i = best.x
    X = errors(E, C, s, i)
    hi = max(x.max() for x in X) + 1e-6
    t, sol, res = spline_fit(X, y, hi, nk)
    return dict(slope=float(s), intercept=float(i), hi=float(hi), knots=t, coef=sol, res=res, X=X)


full = fit(use)
r = full['res']
print('all outputs together: slope %.4f, intercept %.4f Hz | residual of log(term): rms %.4f, largest %.4f | blocks within 0.006: %d of %d' % (
    full['slope'], full['intercept'], np.sqrt(np.mean(r ** 2)), np.abs(r).max(), int((np.abs(r) < 0.006).sum()), len(r)))
for nk in (20, 45):
    f2 = fit(use, nk)
    print('   with %d knots: slope %.4f, intercept %.4f | rms %.4f' % (nk, f2['slope'], f2['intercept'], np.sqrt(np.mean(f2['res'] ** 2))))
# leave one output out
cv = []
for f in files:
    train = [b for b in use if b['file'] != f]
    test = [b for b in use if b['file'] == f]
    m = fit(train, start=(full['slope'], full['intercept']))
    g = BSpline(m['knots'], m['coef'], K)
    y, E, C = setup(test)
    X = errors(E, C, m['slope'], m['intercept'])
    for b, x, yy in zip(test, X, y):
        inside = bool(x.max() <= m['hi'])
        cv.append((f, b['isomer'], float(g(np.clip(x, 0, m['hi'])).sum() - yy), inside, m['slope'], m['intercept']))
ins = np.array([c[2] for c in cv if c[3]])
print('leave one output out: %d blocks predicted from the other outputs (%d more have an error outside the range seen in training)' % (len(ins), len(cv) - len(ins)))
print('   residual of log(term): rms %.4f, median absolute %.4f, largest %.4f | within 0.006: %d, within 0.02: %d, within 0.1: %d' % (
    np.sqrt(np.mean(ins ** 2)), np.median(np.abs(ins)), np.abs(ins).max(), int((np.abs(ins) < 0.006).sum()), int((np.abs(ins) < 0.02).sum()), int((np.abs(ins) < 0.1).sum())))
print('   slope between %.4f and %.4f, intercept between %.4f and %.4f over the %d fits' % (
    min(c[4] for c in cv), max(c[4] for c in cv), min(c[5] for c in cv), max(c[5] for c in cv), len(files)))
worst = sorted(cv, key=lambda c: -abs(c[2]))[:4]
print('   largest:', '; '.join('%s isomer %d %.3f%s' % (c[0], c[1], c[2], '' if c[3] else ' (outside)') for c in worst))
# the density
g = BSpline(full['knots'], full['coef'], K)
allx = np.concatenate(full['X'])
grid = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10]
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'jmodel_density.csv'), 'w') as fh:
    fh.write('absolute_error_Hz,density,couplings_within_0.25_Hz\n')
    for x in np.arange(0, full['hi'], 0.05):
        fh.write('%.2f,%.6g,%d\n' % (x, math.exp(float(g(x))), int(((allx > x - 0.25) & (allx <= x + 0.25)).sum())))
print('recovered density f (per Hz) at absolute corrected errors of', ', '.join('%g' % x for x in grid), 'Hz:')
print('   ' + ', '.join('%.4f' % math.exp(float(g(x))) for x in grid))
print('   couplings behind it: %d, largest error %.2f Hz, %d above 8 Hz' % (len(allx), allx.max(), int((allx > 8).sum())))
json.dump(dict(slope=full['slope'], intercept=full['intercept'], hi=full['hi'], degree=K, knots=[float(v) for v in full['knots']], coef=[float(v) for v in full['coef']],
               gaussian_outputs=sorted(gauss), outputs=files), open(os.path.join(out, 'jmodel_recovered.json'), 'w'), indent=1)
