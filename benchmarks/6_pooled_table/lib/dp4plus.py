"""DP4+ as implemented in the authors' Excel sheet (J. Org. Chem. 2015, 80, 12526).

For one experimental data set and K candidate isomers, per nucleus k:
    unscaled shift  du = TMS - shielding
    scaled shift    ds = (du - b) / m,  m and b from the least-squares line du = m*dexp + b (C and H fitted separately)
    scaled term     T(|dexp - ds| / sigma_s ; nu_s)
    unscaled term   T(|du - dexp - mu| / sigma ; nu), with sp2 or sp3 parameters
    T(x; nu) is the one-tailed Student t probability. Excel's TDIST truncates nu to an integer; `truncate_nu`
    reproduces that (the published 2015 tables were made with the Excel sheet).
sDP4+ uses the scaled terms, uDP4+ the unscaled terms, DP4+ both; each for H, C or all data.
"""
import math
import numpy as np
from scipy import stats


def tail(x, nu, truncate_nu=True):
    nu = math.floor(nu) if truncate_nu else nu
    return stats.t.sf(np.asarray(x, float), nu)


def scale(du, dexp):
    """Regression of calculated (unscaled) on experimental shifts; returns scaled shifts."""
    du, dexp = np.asarray(du, float), np.asarray(dexp, float)
    if len(du) < 2 or np.ptp(dexp) == 0:
        return du.copy()
    m, b = np.polyfit(dexp, du, 1)
    return (du - b) / m


def log_terms(shield, dexp, nuc, sp2, par, truncate_nu=True):
    """shield, dexp: arrays for one candidate (rows with experimental data only); nuc: array of 'C'/'H';
    sp2: boolean array. Returns dict of per-row log terms: scaled ('s') and unscaled ('u')."""
    shield, dexp = np.asarray(shield, float), np.asarray(dexp, float)
    nuc, sp2 = np.asarray(nuc), np.asarray(sp2, bool)
    ls, lu = np.zeros(len(dexp)), np.zeros(len(dexp))
    es, eu = np.zeros(len(dexp)), np.zeros(len(dexp))
    for n in ('C', 'H'):
        idx = np.where(nuc == n)[0]
        if len(idx) == 0:
            continue
        du = par[f'TMS_{n}'] - shield[idx]
        ds = scale(du, dexp[idx])
        e_s = dexp[idx] - ds
        e_u = du - dexp[idx]
        ls[idx] = np.log(np.clip(tail(np.abs(e_s) / par[f'{n}_sc_sigma'], par[f'{n}_sc_nu'], truncate_nu), 1e-300, None))
        for flag, tag in ((True, 'sp2'), (False, 'sp3')):
            sel = sp2[idx] == flag
            if sel.any():
                x = np.abs(e_u[sel] - par[f'{n}_{tag}_mu']) / par[f'{n}_{tag}_sigma']
                lu[idx[sel]] = np.log(np.clip(tail(x, par[f'{n}_{tag}_nu'], truncate_nu), 1e-300, None))
        es[idx], eu[idx] = e_s, e_u
    return dict(s=ls, u=lu, e_s=es, e_u=eu)


def normalise(logp):
    logp = np.asarray(logp, float)
    w = np.exp(logp - logp.max())
    return w / w.sum()


def probabilities(terms_per_candidate, nuc):
    """terms_per_candidate: list of dicts from log_terms. Returns {(method, data): array of K probabilities}."""
    nuc = np.asarray(nuc)
    out = {}
    masks = {'H': nuc == 'H', 'C': nuc == 'C', 'All': np.ones(len(nuc), bool)}
    for data, mk in masks.items():
        if not mk.any():
            continue
        s = np.array([t['s'][mk].sum() for t in terms_per_candidate])
        u = np.array([t['u'][mk].sum() for t in terms_per_candidate])
        out[('sDP4+', data)] = normalise(s)
        out[('uDP4+', data)] = normalise(u)
        out[('DP4+', data)] = normalise(s + u)
    return out
