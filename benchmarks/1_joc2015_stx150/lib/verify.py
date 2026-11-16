"""Is a printed probability reproduced by the printed data?

The SI files print shieldings with two or three decimals. For a probability that is not close to 0 or 100% this
rounding matters: half a unit of the last decimal on every proton shielding moves a probability near 50% by 5 to 10
points. A record is therefore called
  reproduced   when the recomputed probability is within <tol> points of the printed one, and
  consistent   when the printed one lies inside the range that the recomputed probability takes when every printed
               shielding is moved at random within its rounding interval (n draws; the range is widened by the
               rounding of the printed probability itself).
The sampled range is narrower than the full range, so 'not reproduced' can be too severe but 'consistent' cannot be
too lenient."""
import numpy as np
import mmdp4 as M

def rounding_range(D, shieldings, exp, nuc, sp2, par, key, k, half_unit, pairs=None, n=200, seed=0):
    rng = np.random.default_rng(seed); vals = []
    for _ in range(n):
        sh = [np.asarray(s, float) + rng.uniform(-half_unit, half_unit, len(s)) for s in shieldings]
        if pairs: sh = [M.exchange(s, exp, pairs) for s in sh]
        vals.append(100 * float(M.probabilities(D, sh, exp, nuc, sp2, par)[key][k]))
    return min(vals), max(vals)

def status(recomputed, printed, rng_fn, tol=1.0, printed_half_unit=0.5):
    """rng_fn: function returning (lo, hi), called only when needed"""
    if abs(recomputed - printed) <= tol: return 'reproduced', None
    lo, hi = rng_fn()
    if lo - printed_half_unit <= printed <= hi + printed_half_unit: return 'consistent', (lo, hi)
    return 'not reproduced', (lo, hi)
