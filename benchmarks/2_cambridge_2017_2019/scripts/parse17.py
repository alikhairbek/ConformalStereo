"""Reads the PyDP4 outputs (*.dp4) stored in the 2017 deposit.

parse(text) -> dict
  desc      the lines of the experimental description echoed at the top
  isomers   {n: {'C': [(label, calc, scaled, exp, error)],
                 'H': [(label, calc, scaled, exp, error, [J exp], [J calc])],
                 'p': {'C': .., 'H': .., 'J': ..}}}      values of "DP4 based on ..."
  results   {'all' | 'CH' | 'C' | 'H' | 'J': [percent per isomer]}
"""
import re

HEAD = re.compile(r'^Assigned shifts for isomer (\d+):')
ROW = re.compile(r'^([A-Z][a-z]?\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)(.*)$')
BASED = re.compile(r'^DP4 based on (C|H|J):\s*(\S+)')
RES = re.compile(r'^Results of DP4 using (.+?):')
ISO = re.compile(r'^Isomer (\d+):\s*(-?[\d.]+(?:e-?\d+)?|nan)%')
KEYS = {'all available data': 'all', 'both carbon and proton data': 'CH', 'carbon data only': 'C',
        'proton data only': 'H', 'J value data only': 'J'}


def floats(s):
    return [float(x) for x in re.findall(r'-?\d+\.\d+', s)]


def parse(text):
    lines = text.splitlines()
    desc, isomers, results = [], {}, {}
    cur = sect = res = None
    seen = False
    for line in lines:
        m = HEAD.match(line)
        if m:
            cur = isomers.setdefault(int(m.group(1)), {'C': [], 'H': [], 'p': {}})
            sect = res = None
            seen = True
            continue
        if not seen:
            desc.append(line)
            continue
        if line.startswith('Assigned C shifts'):
            sect = 'C'
            continue
        if line.startswith('Assigned H shifts'):
            sect = 'H'
            continue
        m = RES.match(line)
        if m:
            res = results.setdefault(KEYS.get(m.group(1), m.group(1)), [])
            cur = sect = None
            continue
        if res is not None:
            m = ISO.match(line)
            if m:
                v = m.group(2)
                res.append(float('nan') if v == 'nan' else float(v))
            elif line.strip():
                res = None
            continue
        if cur is None:
            continue
        m = BASED.match(line)
        if m:
            try:
                cur['p'][m.group(1)] = float(m.group(2))
            except ValueError:
                cur['p'][m.group(1)] = float('nan')
            sect = None
            continue
        if line.startswith(('Max ', 'Direct J', 'Mean absolute')):
            sect = None
            continue
        m = ROW.match(line)
        if m and sect:
            row = (m.group(1),) + tuple(float(m.group(i)) for i in (2, 3, 4, 5))
            if sect == 'H':
                tail = m.group(6)
                if '||' in tail:
                    a, b = tail.split('||', 1)
                    row = row + (floats(a), floats(b))
                else:
                    row = row + ([], [])
            cur[sect].append(row)
    return {'desc': desc, 'isomers': isomers, 'results': results}
