"""Experimental descriptions of the 2017 deposit (PyDP4 format of 2015-2017).

line 1   carbon shifts,   e.g.  172.7(),106.0(C4),...       a label in brackets = assigned atom
line 3   proton shifts,   e.g.  5.84(),3.38(;7.4,5.5),...   after ';' the coupling constants
then     groups of equivalent atoms, one per line (H17,H19,H18), and a line 'OMIT H23,H33' with the atoms to leave out
"""
import re

ITEM = re.compile(r'(-?\d+(?:\.\d+)?)\(([^)]*)\)')


def shifts(line):
    out = []
    for value, inside in ITEM.findall(line):
        label, _, js = inside.partition(';')
        out.append((float(value), label.strip(), [float(x) for x in re.findall(r'-?\d+(?:\.\d+)?', js)]))
    return out


def parse(lines):
    """lines of the description -> dict(C=[(shift, label, J)], H=[...], groups=[[labels]], omit=[labels])"""
    blocks, cur = [], []
    for line in lines:
        if line.strip():
            cur.append(line.strip())
        elif cur:
            blocks.append(cur)
            cur = []
    if cur:
        blocks.append(cur)
    d = {'C': [], 'H': [], 'groups': [], 'omit': []}
    if blocks:
        d['C'] = shifts(' '.join(blocks[0]))
    if len(blocks) > 1:
        d['H'] = shifts(' '.join(blocks[1]))
    for block in blocks[2:]:
        for line in block:
            if line.upper().startswith('OMIT'):
                d['omit'] += [x.strip() for x in line[4:].split(',') if x.strip()]
            else:
                d['groups'].append([x.strip() for x in line.split(',') if x.strip()])
    return d
