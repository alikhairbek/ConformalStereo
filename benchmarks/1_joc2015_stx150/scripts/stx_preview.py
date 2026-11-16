"""DP4+ confidence and split-conformal sets on the STX150 records, from the recomputed probabilities of every printed
candidate (stx_candidates.csv). Primary set: verified, close and corrected records. Sensitivity: the strict set
(verified and corrected only) and the set with the unverified records added.
Conformal: score = 1 - probability; calibration and test halves drawn at random by candidate set (cases that share
their candidates stay together); threshold = the ceil((n+1)(1-alpha))-th smallest calibration score of the correct
isomers; a candidate is kept when its score is at or below the threshold.
Usage: stx_preview.py <out dir with stx_records.csv and stx_candidates.csv> [n_splits] [seed]"""
import sys, os, csv, math, random, collections
out = sys.argv[1]; n_splits = int(sys.argv[2]) if len(sys.argv) > 2 else 2000; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261114
rec = {(r['case'], r['level']): r for r in csv.DictReader(open(os.path.join(out, 'stx_records.csv')))}
cand = collections.defaultdict(list)
for r in csv.DictReader(open(os.path.join(out, 'stx_candidates.csv'))): cand[(r['case'], r['level'])].append((int(r['is_correct']), float(r['probability']) / 100))
LEVELS = ['MMFF_PBE', 'xTB_PBE', 'xTB_mPW', 'xTB_wB9', 'B3L_PBE', 'B3L_wB9']
def select(statuses, lev): return [k for k, r in rec.items() if r['level'] == lev and r['status'] in statuses and k in cand]
def confidence(keys):
    n = len(keys); top = call = call_ok = low = 0
    for k in keys:
        ps = cand[k]; pt = [p for t, p in ps if t][0]; best = max(p for _, p in ps)
        top += pt >= best - 1e-12; low += pt <= 0.01
        if best >= 0.99: call += 1; call_ok += pt >= 0.99
    return n, top, call_ok, call, low
def conformal(keys, alpha, rng):
    units = sorted({rec[k]['candidate_set'] for k in keys}); member = collections.defaultdict(list)
    for k in keys: member[rec[k]['candidate_set']].append(k)
    cov = kept = total = single = single_ok = n = 0
    for _ in range(n_splits):
        rng.shuffle(units); h = len(units) // 2
        cal = [k for u in units[:h] for k in member[u]]; test = [k for u in units[h:] for k in member[u]]
        s = sorted(1 - [p for t, p in cand[k] if t][0] for k in cal); j = math.ceil((len(s) + 1) * (1 - alpha))
        q = s[j - 1] if j <= len(s) else float('inf')
        for k in test:
            keep = [(t, p) for t, p in cand[k] if 1 - p <= q + 1e-12]
            cov += any(t for t, _ in keep); kept += len(keep); total += len(cand[k]); n += 1
            if len(keep) == 1: single += 1; single_ok += keep[0][0]
    return 100 * cov / n, 100 * kept / total, 100 * single / n, 100 * single_ok / max(single, 1)
for name, statuses in (('verified, close and corrected records', ('verified', 'close', 'corrected')), ('verified and corrected records only', ('verified', 'corrected')), ('with the unverified records added', ('verified', 'close', 'corrected', 'unverified'))):
    print(name)
    print('  level     records  sets | first  | calls >=99%%: right/all | correct <=1%% || 90%% target: coverage, candidates kept, one left, of which right || 80%% target: same four')
    tot = [0] * 5
    for lev in LEVELS:
        keys = select(statuses, lev); n, top, cok, call, low = confidence(keys)
        rng = random.Random(seed); a = conformal(keys, 0.10, rng); rng = random.Random(seed); b = conformal(keys, 0.20, rng)
        for i, v in enumerate((n, top, cok, call, low)): tot[i] += v
        print('  %-9s %5d %5d | %3d    | %3d/%-3d                 | %3d          || %5.1f %5.1f %5.1f %5.1f || %5.1f %5.1f %5.1f %5.1f' % ((lev, n, len({rec[k]['candidate_set'] for k in keys}), top, cok, call, low) + a + b))
    print('  %-9s %5d       | %3d (%.1f%%) | %3d/%-3d (%.1f%%)       | %3d (%.1f%%)' % ('pooled', tot[0], tot[1], 100 * tot[1] / tot[0], tot[2], tot[3], 100 * tot[2] / tot[3], tot[4], 100 * tot[4] / tot[0]))
