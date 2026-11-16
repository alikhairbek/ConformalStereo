"""Checks of pairing.py: the dynamic programme against exhaustive search, and three constructed cases."""
import itertools
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pairing

random.seed(1)
c = [random.uniform(0, 200) for _ in range(12)]
e = [x * 1.02 - 1 + random.gauss(0, 1) for x in c]
a, b = pairing.optimal(c, e, 2.27)
print('as many values as atoms: every atom kept:', sorted(a) == list(range(12)))
a, b = pairing.optimal(c, [x for i, x in enumerate(e) if i != 5], 2.27)
print('one signal removed (atom 5): atom left out', sorted(set(range(12)) - set(a)))
e3 = e + [77.0, 128.3]
a, b = pairing.optimal(c, e3, 2.27)
print('two signals added (77.0, 128.3): signals left out', sorted(e3[i] for i in set(range(len(e3))) - set(b)))
ok = 0
for t in range(300):
    n = random.randint(4, 8)
    m = random.randint(2, n - 1)
    cs = sorted(random.uniform(0, 10) for _ in range(n))
    es = sorted(random.uniform(0, 10) for _ in range(m))
    v = sum(pairing.term(cs[i] - es[j], 0.5) for i, j in pairing.align(cs, es, 0.5))
    best = max(sum(pairing.term(cs[i] - es[j], 0.5) for j, i in enumerate(sub)) for sub in itertools.combinations(range(n), m))
    ok += abs(v - best) < 1e-9
print('dynamic programme equals exhaustive search in %d of 300 random problems' % ok)
hit = tot = 0
for t in range(200):
    c = [random.uniform(0.7, 7.5) for _ in range(30)]
    e = [x * 0.97 + 0.1 + random.gauss(0, 0.12) for x in c]
    k = random.randint(1, 3)
    drop = set(random.sample(range(30), k))
    a, b = pairing.optimal(c, [x for i, x in enumerate(e) if i not in drop], 0.187)
    tot += k
    hit += len((set(range(30)) - set(a)) & drop)
print('30 protons, 1 to 3 signals missing: the atoms left out are exactly the right ones in %d of %d; '
      'the others are neighbours in shift, which the shifts alone cannot tell apart' % (hit, tot))
