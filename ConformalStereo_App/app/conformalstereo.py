#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ConformalStereo: calibrated candidate sets for NMR stereochemical assignments by DP4-type methods.

A DP4, DP4+, MM-DP4+ or xTB-DP4+ analysis gives every candidate structure a likelihood L. This module turns
these likelihoods into a set of candidates that contains the actual structure at a chosen rate:

    score of candidate i     s_i = (ln L_max - ln L_i) / N              N = number of shifts
    calibrated set           every candidate with s_i <= c              c = cut-off of the workflow

The cut-off c is a family-weighted conformal quantile of the scores that the actual structures of molecules of
known structure receive with the same workflow. Because ln L_max - ln L_i = ln(P_best / P_i), the unrounded
probabilities of a DP4-type output are all that is needed.

Only the Python standard library is used (Python 3.8 or newer). Everything here can be imported, and the file
is also a command line:  python conformalstereo.py --help
"""
import sys, os, re, io, csv, json, math, zipfile, argparse, collections
import xml.etree.ElementTree as ET

__version__ = '1.0.0'
INF = float('inf')
LEVELS = (0.80, 0.90, 0.95)
HERE = os.path.dirname(os.path.abspath(__file__))
# Single Gaussian of the PyDP4 program for internally scaled shifts (ppm); the DP4 workflows of the paper use it.
SIGMA = {'C': 2.269372270818724, 'H': 0.18731058105269952}


class InputError(ValueError):
    """the input cannot be used; the message is written for the person who supplied it"""


# ---------------------------------------------------------------------------------------------------------------
# numbers, with the precision they were written with
# ---------------------------------------------------------------------------------------------------------------
_NUMBER = re.compile(r'^([+-]?)(\d*)(?:\.(\d*))?(?:[eE]([+-]?\d+))?$')
_LAST_NUMBER = re.compile(r'([+\-\u2212]?(?:\d+[.,]?\d*|[.,]\d+)(?:[eE][+\-]?\d+)?|[+\-\u2212]?inf(?:inity)?|nan)\s*%?\s*$', re.I)


def parse_number(text):
    """'99.87' -> (99.87, 0.005): the value and half a unit of its last written digit"""
    t = text.strip().rstrip('%').strip().replace('\u2212', '-').replace(',', '.')
    low = t.lower()
    if low in ('inf', '+inf', 'infinity', '+infinity'):
        return INF, 0.0
    if low in ('-inf', '-infinity'):
        return -INF, 0.0
    m = _NUMBER.match(t)
    if not m or (not m.group(2) and not m.group(3)):
        raise InputError('"%s" is not a number' % text.strip())
    sign, whole, frac, expo = m.groups()
    frac = frac or ''
    expo = int(expo) if expo else 0
    value = float('%s%s.%se%d' % (sign, whole or '0', frac or '0', expo))
    return value, 0.5 * 10.0 ** (expo - len(frac))


def _ln(x):
    return math.log(x) if x > 0 else -INF


class Molecule(object):
    """The candidates of one molecule with bounds on their log-likelihoods.

    names   candidate names
    lo, hi  bounds of ln L for every candidate, up to a constant shared by all candidates. They differ when the
            input was rounded (printed percentages, printed errors); for unrounded input they coincide.
    mid     ln L from the values as written (for display)
    shifts  number of shifts N in the likelihood, or None if the source does not tell
    """

    def __init__(self, names, lo, hi, mid=None, shifts=None, source='', shown=None, info=None, notes=None):
        if len(names) < 1:
            raise InputError('No candidate was found.')
        if len(set(names)) != len(names):
            raise InputError('Two candidates have the same name: %s.' % ', '.join(sorted(set(n for n in names if names.count(n) > 1))))
        self.names = list(names)
        self.lo, self.hi = [float(x) for x in lo], [float(x) for x in hi]
        self.mid = [float(x) for x in (mid if mid is not None else [(a + b) / 2.0 if a > -INF else b for a, b in zip(lo, hi)])]
        for a, b in zip(self.lo, self.hi):
            if a != a or b != b or a > b:
                raise InputError('A value could not be read as a number.')
        if max(self.hi) == -INF or max(self.mid) == -INF:
            raise InputError('Every candidate has a probability of zero. At least one positive value is needed.')
        self.shifts = shifts
        self.source = source
        self.shown = list(shown) if shown is not None else ['' for _ in names]
        self.info = dict(info or {})
        self.notes = list(notes or [])

    @property
    def exact(self):
        """True if the values carry enough digits that rounding cannot matter (ln L known to 0.001 or better)"""
        return all((b - a <= 1e-3) if a > -INF else (b == -INF) for a, b in zip(self.lo, self.hi))

    def probabilities(self):
        """probabilities from the values as written (they add up to 1)"""
        m = max(self.mid)
        w = [math.exp(x - m) if x > -INF else 0.0 for x in self.mid]
        s = sum(w)
        return [x / s for x in w]

    def best(self):
        return max(range(len(self.names)), key=lambda i: (self.mid[i], -i))


def molecule_from_values(names, texts, kind='probability', shifts=None, source='typed values'):
    """kind 'probability': numbers on any common scale (per cent or fractions); kind 'loglik': natural logarithms"""
    lo, hi, mid = [], [], []
    for name, text in zip(names, texts):
        v, h = parse_number(text)
        if v != v:
            raise InputError('The value of %s is not a number.' % name)
        if kind == 'loglik':
            lo.append(v - h); hi.append(v + h); mid.append(v)
        else:
            if v < 0:
                raise InputError('The probability of %s is negative.' % name)
            if v == INF:
                raise InputError('The probability of %s is infinite.' % name)
            lo.append(_ln(v - h)); hi.append(_ln(v + h)); mid.append(_ln(v))
    return Molecule(names, lo, hi, mid, shifts=shifts, source=source, shown=[t.strip() for t in texts])


def read_typed(text, kind='probability', shifts=None):
    """One candidate per line: 'name value' or 'name, value' or just the value. Header lines without a number are skipped."""
    names, values = [], []
    for raw in text.replace('\r', '').split('\n'):
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        m = _LAST_NUMBER.search(line)
        if not m:
            if not names:
                continue                                    # a header line
            raise InputError('No number was found in the line "%s".' % line)
        name = line[:m.start()].strip().rstrip(',;:=\t').strip()
        names.append(name or 'candidate %d' % (len(names) + 1))
        values.append(m.group(1))
    if len(names) < 2:
        raise InputError('At least two candidates are needed, one per line: a name followed by its probability.')
    return molecule_from_values(names, values, kind=kind, shifts=shifts)


# ---------------------------------------------------------------------------------------------------------------
# the method
# ---------------------------------------------------------------------------------------------------------------
def score_bounds(mol, shifts=None):
    """-> [(lowest, as written, highest)] of s_i = (ln L_max - ln L_i) / N for every candidate"""
    n = shifts if shifts is not None else mol.shifts
    if not n or n <= 0:
        raise InputError('The number of shifts is needed (a whole number above zero).')
    k = len(mol.names)
    top = max(mol.mid)
    out = []
    for i in range(k):
        o_lo = max([mol.lo[j] for j in range(k) if j != i] or [-INF])
        o_hi = max([mol.hi[j] for j in range(k) if j != i] or [-INF])
        s_lo = 0.0 if o_lo == -INF else (INF if mol.hi[i] == -INF else max(0.0, o_lo - mol.hi[i]) / n)
        s_hi = 0.0 if o_hi == -INF else (INF if mol.lo[i] == -INF else max(0.0, o_hi - mol.lo[i]) / n)
        s = INF if mol.mid[i] == -INF else (top - mol.mid[i]) / n
        out.append((min(s_lo, s), s, max(s_hi, s)))
    return out


def cutoff(items, level):
    """Family-weighted conformal cut-off. items: [(score of the actual structure, family)].

    A score weighs one over the number of molecules of its family, every family weighs one over the number of
    families plus one, and the remaining weight stands for the family of the new molecule. The cut-off is the
    smallest score at which the weight reaches the level; it is infinite if the families are too few."""
    if not items:
        return INF
    count = collections.Counter(f for _, f in items)
    g = len(count)
    acc = 0.0
    for s, f in sorted(items, key=lambda x: x[0]):
        acc += 1.0 / count[f] / (g + 1)
        if acc >= level - 1e-12:
            return s
    return INF


def families_needed(level):
    """smallest number of families with which a level can be reached: 4, 9 and 19 for 80, 90 and 95%"""
    return int(math.ceil(level / (1.0 - level) - 1e-9))


IN, OUT, OPEN = 'in', 'out', 'open'


def conformal_set(bounds, c):
    """-> status of every candidate at cut-off c: 'in', 'out', or 'open' when the digits given do not decide"""
    out = []
    for s_lo, s, s_hi in bounds:
        if s_hi <= c:
            out.append(IN)
        elif s_lo > c:
            out.append(OUT)
        else:
            out.append(OPEN)
    return out


def level_text(level):
    return ('%.1f' % (100 * level)).rstrip('0').rstrip('.') + '%'


def join_names(names):
    names = list(names)
    return names[0] if len(names) == 1 else ', '.join(names[:-1]) + ' and ' + names[-1]


def assign(mol, cutoffs, shifts=None, workflow='', calibration=None, origin=None):
    """The calibrated sets of one molecule.

    cutoffs      {level: cut-off per shift}, e.g. {0.8: 0.208, 0.9: 0.461}; None or inf = level not supported
    workflow     name of the workflow of the paper, if the cut-offs are taken from it
    origin       where the cut-offs come from, as words for the report (default: the workflow, or "the user")
    calibration  optional dict with 'shifts_max' and 'candidates_max' of the calibration molecules, for warnings
    -> dict with the scores, the sets and sentences that can be quoted"""
    n = shifts if shifts is not None else mol.shifts
    bounds = score_bounds(mol, n)
    k = len(mol.names)
    if origin is None:
        origin = ('the workflow ' + workflow) if workflow else 'the user'
    probs = mol.probabilities()
    best = mol.best()
    res = dict(shifts=n, candidates=k, workflow=workflow, origin=origin, source=mol.source, exact=mol.exact, best=mol.names[best],
               rows=[dict(name=mol.names[i], shown=mol.shown[i], probability=probs[i], score=bounds[i][1], score_lo=bounds[i][0],
                          score_hi=bounds[i][2], log_ratio=(mol.mid[i] - max(mol.mid)) if mol.mid[i] > -INF else -INF) for i in range(k)],
               levels=[], warnings=list(mol.notes), info=mol.info)
    for level in sorted(cutoffs):
        c = cutoffs[level]
        if c is None or c == INF:
            res['levels'].append(dict(level=level, label=level_text(level), supported=False, cutoff=None, status=[IN] * k, members=list(mol.names), open=[],
                                      size=k, sentence='The %s level cannot be reached with this calibration (at least %d families are needed); every candidate remains.'
                                      % (level_text(level), families_needed(level))))
            continue
        status = conformal_set(bounds, c)
        members = [mol.names[i] for i in range(k) if status[i] != OUT]
        still_open = [mol.names[i] for i in range(k) if status[i] == OPEN]
        where = ' (cut-off %.3f per shift from %s; %d shifts, %d candidates)' % (c, origin, n, k)
        if len(members) == 1:
            text = 'At the %s level the calibrated set is one structure, %s%s.' % (level_text(level), members[0], where)
        elif len(members) == k:
            text = 'At the %s level the calibrated set keeps %s candidates%s: these data do not decide between them.' % (
                level_text(level), 'both' if k == 2 else 'all %d' % k, where)
        else:
            text = 'At the %s level the calibrated set holds %d of the %d candidates, %s%s.' % (level_text(level), len(members), k, join_names(members), where)
        if still_open:
            text += ' For %s the digits given do not decide; %s kept.' % (join_names(still_open), 'it is' if len(still_open) == 1 else 'they are')
        res['levels'].append(dict(level=level, label=level_text(level), supported=True, cutoff=c, limit=c * n, status=status, members=members,
                                  open=still_open, size=len(members), sentence=text))
    if any(l['open'] for l in res['levels']):
        res['warnings'].append('Some values were rounded too far to decide whether their candidate belongs to a set. Such candidates are kept, '
                               'so the sets may be larger than necessary. Use unrounded probabilities: read the DP4+ spreadsheet or the '
                               'PyDP4 output directly, or copy the values with all their digits.')
    if calibration:
        if calibration.get('shifts_max') and n > calibration['shifts_max']:
            res['warnings'].append('This molecule has more shifts (%d) than any calibration molecule of the workflow (at most %d). Report the 90%% set.'
                                   % (n, calibration['shifts_max']))
        if calibration.get('candidates_max') and k > calibration['candidates_max']:
            res['warnings'].append('This molecule has more candidates (%d) than any calibration molecule of the workflow (at most %d). Report the 90%% set.'
                                   % (k, calibration['candidates_max']))
    if k > 4:
        res['warnings'].append('With more than four candidates the 80% sets of the benchmarks contained the actual structure less often than stated '
                               '(66 to 74%). Report the 90% set.')
    res['caveat'] = ('A set of this kind contains the actual structure, on average, for at least the stated share of molecules that are like the '
                     'calibration molecules and are computed with the same workflow. It is conditional on the list of candidates: if the actual '
                     'structure was not proposed, the set cannot contain it.')
    return res


# ---------------------------------------------------------------------------------------------------------------
# calibration with molecules of known structure
# ---------------------------------------------------------------------------------------------------------------
_COLUMNS = dict(
    molecule=('molecule', 'case', 'compound', 'name', 'record', 'experiment'),
    family=('family', 'group', 'series'),
    shifts=('shifts', 'n_nuclei', 'n', 'nuclei', 'number_of_shifts', 'n_shifts'),
    candidate=('candidate', 'isomer', 'structure', 'diastereomer'),
    actual=('actual', 'is_true', 'true', 'correct', 'is_actual', 'is_correct', 'right'),
    probability=('probability', 'prob', 'p', 'dp4', 'dp4+', 'dp4plus', 'percent', 'probability_percent'),
    loglik=('log_likelihood', 'loglik', 'log_ratio', 'lnl', 'ln_l', 'logl', 'log_l'))
_TRUE = ('1', 'true', 'yes', 'y', 'x', 't', 'actual', 'right', 'correct')
_FALSE = ('0', 'false', 'no', 'n', '', 'f', '-')


def read_table(text):
    """A table with one row per candidate -> list of (molecule name, family, Molecule, index of the actual structure or None).

    Columns (comma, semicolon or tab separated; any order): molecule, family (optional), shifts, candidate,
    probability or log_likelihood, actual (1 for the actual structure). The pooled table of the paper has this
    form (case, family, n_nuclei, candidate, log_ratio, is_true)."""
    text = text.replace('\r\n', '\n').replace('\r', '\n').lstrip('\ufeff')
    lines = [l for l in text.split('\n') if l.strip() and not l.lstrip().startswith('#')]
    if len(lines) < 2:
        raise InputError('The table is empty. It needs a header line and one line per candidate.')
    head = lines[0]
    delim = max(('\t', ';', ','), key=lambda d: head.count(d))
    if head.count(delim) == 0:
        raise InputError('The columns of the table could not be separated. Use commas, semicolons or tabs between the columns.')
    rows = list(csv.reader(lines, delimiter=delim))
    header = [h.strip().lower().replace(' ', '_') for h in rows[0]]
    col = {}
    for key, names in _COLUMNS.items():
        for nm in names:
            if nm in header:
                col[key] = header.index(nm)
                break
    for need in ('molecule', 'shifts', 'candidate'):
        if need not in col:
            raise InputError('The table has no column "%s". Columns found: %s.' % (need, ', '.join(rows[0])))
    if 'probability' not in col and 'loglik' not in col:
        raise InputError('The table needs a column "probability" or "log_likelihood". Columns found: %s.' % ', '.join(rows[0]))
    kind = 'loglik' if 'loglik' in col else 'probability'
    vcol = col[kind]
    groups = collections.OrderedDict()
    for ln, r in enumerate(rows[1:], 2):
        if len(r) < len(header):
            r = r + [''] * (len(header) - len(r))
        name = r[col['molecule']].strip()
        if not name:
            raise InputError('Line %d of the table has no molecule name.' % ln)
        groups.setdefault(name, []).append((ln, r))
    out = []
    for name, rr in groups.items():
        fams = set(r[col['family']].strip() for _, r in rr) if 'family' in col else set([''])
        if len(fams) > 1:
            raise InputError('Molecule %s is given more than one family: %s.' % (name, ', '.join(sorted(fams))))
        family = fams.pop() or name
        ns = set(r[col['shifts']].strip() for _, r in rr)
        if len(ns) != 1:
            raise InputError('Molecule %s has different numbers of shifts in its lines: %s.' % (name, ', '.join(sorted(ns))))
        try:
            n = int(float(ns.pop()))
        except ValueError:
            raise InputError('The number of shifts of molecule %s is not a number.' % name)
        if n <= 0:
            raise InputError('The number of shifts of molecule %s must be above zero.' % name)
        actual = None
        if 'actual' in col:
            flags = []
            for ln, r in rr:
                f = r[col['actual']].strip().lower()
                if f in _TRUE:
                    flags.append(True)
                elif f in _FALSE:
                    flags.append(False)
                else:
                    raise InputError('Line %d: the column "actual" holds "%s"; use 1 for the actual structure and 0 otherwise.' % (ln, r[col['actual']]))
            if sum(flags) > 1:
                raise InputError('Molecule %s has more than one candidate marked as the actual structure.' % name)
            actual = flags.index(True) if any(flags) else None
        mol = molecule_from_values([r[col['candidate']].strip() or 'candidate %d' % (i + 1) for i, (_, r) in enumerate(rr)],
                                   [r[vcol] for _, r in rr], kind=kind, shifts=n, source='table')
        out.append((name, family, mol, actual))
    return out


def calibrate(entries, levels=LEVELS):
    """entries: list of (molecule name, family, Molecule, index of the actual structure) -> cut-offs and a summary.

    The score of an actual structure whose value was rounded is taken at its upper bound, which can only raise
    the cut-off."""
    items, rows = [], []
    for name, family, mol, actual in entries:
        if actual is None:
            raise InputError('Molecule %s has no candidate marked as the actual structure.' % name)
        b = score_bounds(mol)
        order = sorted(range(len(mol.names)), key=lambda i: (-mol.mid[i], i))
        items.append((b[actual][2], family))
        rows.append(dict(molecule=name, family=family, shifts=mol.shifts, candidates=len(mol.names), actual=mol.names[actual],
                         score=b[actual][1], score_used=b[actual][2], rank=order.index(actual) + 1, bounds=b,
                         exact=(b[actual][2] - b[actual][0]) <= 1e-4))
    fams = collections.Counter(f for _, f in items)
    g = len(fams)
    res = dict(molecules=len(rows), families=g, rows=rows, levels=[],
               shifts_min=min(r['shifts'] for r in rows), shifts_max=max(r['shifts'] for r in rows),
               candidates_min=min(r['candidates'] for r in rows), candidates_max=max(r['candidates'] for r in rows),
               top_right=sum(r['rank'] == 1 for r in rows), warnings=[])
    for level in levels:
        need = families_needed(level)
        c = cutoff(items, level)
        entry = dict(level=level, label=level_text(level), families_needed=need, supported=c < INF, cutoff=c if c < INF else None)
        if c < INF and g >= 2:                              # what the sets look like on these molecules, their own family left out
            cov = size = single = 0.0
            for (name, family, mol, actual), r in zip(entries, rows):
                c1 = cutoff([it for it in items if it[1] != family], level)
                st = conformal_set(r['bounds'], c1)
                w = 1.0 / fams[family] / g
                cov += w * (st[actual] != OUT)
                size += w * sum(s != OUT for s in st)
                single += w * (sum(s != OUT for s in st) == 1)
            entry.update(coverage_left_out=cov, mean_set_size=size, single_structure_sets=single)
        res['levels'].append(entry)
    if g < families_needed(min(levels)):
        res['warnings'].append('%d %s: at least %d are needed for the %s level, %d for %s and %d for %s.' % (
            g, 'family is too few' if g == 1 else 'families are too few', families_needed(0.8), level_text(0.8), families_needed(0.9),
            level_text(0.9), families_needed(0.95), level_text(0.95)))
    if any(not r['exact'] for r in rows):
        res['warnings'].append('Some calibration molecules were given with rounded values. Their scores are taken at the upper bound, which can '
                               'make the cut-offs larger than necessary.')
    if any(r['score_used'] == INF for r in rows):
        res['warnings'].append('For %d calibration molecule(s) the actual structure has a probability that was rounded to zero. Such a molecule '
                               'counts as missed at any cut-off; supply unrounded values.' % sum(r['score_used'] == INF for r in rows))
    for r in rows:
        del r['bounds']
    return res


# ---------------------------------------------------------------------------------------------------------------
# the cut-offs of the 76 workflows of the paper
# ---------------------------------------------------------------------------------------------------------------
def load_workflows(path=None):
    path = path or os.path.join(HERE, 'cutoffs.csv')
    out = []
    with io.open(path, encoding='utf-8', newline='') as f:
        for r in csv.DictReader(f):
            w = dict(id=r['id'], benchmark=r['benchmark'], method=r['method'], geometry=r['geometry'], level=r['level'],
                     records=int(r['records']), families=int(r['families']), shifts_min=int(r['shifts_min']), shifts_max=int(r['shifts_max']),
                     candidates_max=int(r['candidates_max']), top_right=float(r['top_right']), cutoffs={})
            for lv, key in ((0.8, 'cutoff_80'), (0.9, 'cutoff_90'), (0.95, 'cutoff_95')):
                w['cutoffs'][lv] = float(r[key]) if r[key] != '' else None
            w['label'] = '%s: %s on %s geometries, %s' % (w['benchmark'], w['method'], w['geometry'], w['level'])
            out.append(w)
    return out


def find_workflow(key, workflows=None):
    workflows = workflows if workflows is not None else load_workflows()
    exact = [w for w in workflows if w['id'] == key]
    if exact:
        return exact[0]
    words = key.lower().split()
    hits = [w for w in workflows if all(x in (w['id'] + ' ' + w['label']).lower() for x in words)]
    if len(hits) == 1:
        return hits[0]
    if not hits:
        raise InputError('No workflow matches "%s". List them with: python conformalstereo.py workflows' % key)
    raise InputError('"%s" matches %d workflows: %s' % (key, len(hits), '; '.join(w['id'] for w in hits[:8]) + (' ...' if len(hits) > 8 else '')))


# ---------------------------------------------------------------------------------------------------------------
# reading workbooks (.xlsx) with the standard library
# ---------------------------------------------------------------------------------------------------------------
_NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
_RNS = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
_MAX_MEMBER = 80 * 1024 * 1024


def _col_index(letters):
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n


class Workbook(object):
    """sheet names and, on demand, the stored values of their cells: {(row, column): float or str}; both start at 1"""

    def __init__(self, data):
        try:
            self.zip = zipfile.ZipFile(io.BytesIO(data))
            names = set(self.zip.namelist())
            if 'xl/workbook.xml' not in names:
                raise InputError('This file is not an Excel workbook (.xlsx).')
            self.shared = []
            if 'xl/sharedStrings.xml' in names:
                for si in ET.fromstring(self._read('xl/sharedStrings.xml')).findall(_NS + 'si'):
                    self.shared.append(''.join(t.text or '' for t in si.iter(_NS + 't')))
            rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(self._read('xl/_rels/workbook.xml.rels'))}
            self.sheets = collections.OrderedDict()
            for s in ET.fromstring(self._read('xl/workbook.xml')).find(_NS + 'sheets'):
                target = rels.get(s.get(_RNS + 'id'), '')
                self.sheets[s.get('name')] = target[1:] if target.startswith('/') else 'xl/' + target
        except (zipfile.BadZipFile, ET.ParseError, KeyError):
            raise InputError('This file could not be read as an Excel workbook (.xlsx). Files in the old .xls format have to be saved as .xlsx first.')
        self._cells = {}

    def _read(self, name):
        if self.zip.getinfo(name).file_size > _MAX_MEMBER:
            raise InputError('The workbook is too large to be read here.')
        return self.zip.read(name)

    def cells(self, sheet):
        if sheet not in self._cells:
            out = {}
            root = ET.fromstring(self._read(self.sheets[sheet]))
            data = root.find(_NS + 'sheetData')
            for row in (data if data is not None else []):
                for c in row.findall(_NS + 'c'):
                    ref = c.get('r') or ''
                    m = re.match(r'^([A-Z]+)(\d+)$', ref)
                    if not m:
                        continue
                    t = c.get('t')
                    v = c.find(_NS + 'v')
                    if t == 'inlineStr':
                        val = ''.join(x.text or '' for x in c.iter(_NS + 't'))
                    elif v is None or v.text is None:
                        continue
                    elif t == 's':
                        val = self.shared[int(v.text)]
                    elif t in ('str', 'e'):
                        val = v.text
                    elif t == 'b':
                        val = float(v.text)
                    else:
                        try:
                            val = float(v.text)
                        except ValueError:
                            val = v.text
                    out[(int(m.group(2)), _col_index(m.group(1)))] = val
            self._cells[sheet] = out
        return self._cells[sheet]


def _is_num(x):
    return isinstance(x, float) and x == x


def read_dp4plus_spreadsheet(wb):
    """The Excel sheet of DP4+ (Grimblat, Zanardi and Sarotti, J. Org. Chem. 2015): the stored DP4+ probabilities of the isomers"""
    for sheet in wb.sheets:
        cells = wb.cells(sheet)
        head = [(r, c) for (r, c), v in cells.items() if isinstance(v, str) and v.strip().lower() == 'experimental']
        rows = [(r, c) for (r, c), v in cells.items() if isinstance(v, str) and v.strip() == 'DP4+']
        if not head or not rows:
            continue
        hr, hc = sorted(head)[0]
        isomers = sorted((c, v.strip()) for (r, c), v in cells.items() if r == hr and c > hc and isinstance(v, str) and v.strip())
        pr = sorted(rows)[0][0]
        names, values = [], []
        for c, name in isomers:
            v = cells.get((pr, c))
            if _is_num(v):
                names.append(name); values.append(v)
        n = sum(1 for (r, c), v in cells.items() if c == hc and r > hr and _is_num(v))
        if not names:
            raise InputError('The DP4+ sheet "%s" holds no probabilities yet (the row "DP4+" is empty). Fill in the data in Excel, save the file '
                             'and read it again. If only scaled shifts were entered, the sheet gives sDP4+ alone, for which the paper has no cut-offs.' % sheet)
        if len(names) < 2:
            raise InputError('The DP4+ sheet "%s" holds one isomer only; at least two candidates are needed.' % sheet)
        if any(v < 0 for v in values):
            raise InputError('The DP4+ sheet "%s" holds a negative probability.' % sheet)
        info = {}
        for label in ('Functional', 'Solvent?', 'Basis Set', 'Type of Data'):
            for (r, c), v in cells.items():
                if isinstance(v, str) and v.strip() == label and isinstance(cells.get((r + 1, c)), str):
                    info[label.rstrip('?')] = cells[(r + 1, c)].strip()
        lo = [_ln(v) for v in values]
        mol = Molecule(names, lo, lo, lo, shifts=n or None, source='DP4+ spreadsheet, sheet "%s"' % sheet,
                       shown=['%.6g%%' % (100 * v) for v in values], info=info)
        if any(v == 0 for v in values):
            mol.notes.append('The spreadsheet stores a probability of exactly zero for %s; such a candidate is outside every set.'
                             % join_names(nm for nm, v in zip(names, values) if v == 0))
        return mol
    return None


def read_dp4plus_app(wb):
    """Results workbook of DP4+App (Franco, Zanardi and Sarotti, J. Nat. Prod. 2023), written for version 2.2.

    The sheet 'results' prints the probabilities with four decimals. The sheets 'p_sca' and 'p_uns' hold the
    probability of every nucleus with ten decimals; their product is the unrounded likelihood."""
    if 'results' not in wb.sheets:
        return None
    res = wb.cells('results')
    full = [(r, c) for (r, c), v in res.items() if isinstance(v, str) and v.strip() == 'Full']
    if not full:
        return None
    fr, fc = sorted(full)[0]
    labels = ('H_sca', 'C_sca', 'Sca', 'H_uns', 'C_uns', 'Uns', 'H_full', 'C_full', 'Full')
    hr = min(r for (r, c), v in res.items() if c == fc and isinstance(v, str) and v.strip() in labels) - 1

    def text(v):
        return str(int(v)) if _is_num(v) and v == int(v) else str(v).strip()
    isomers = sorted((c, text(v)) for (r, c), v in res.items() if r == hr and c > fc and text(v))
    printed = [(name, res.get((fr, c))) for c, name in isomers if _is_num(res.get((fr, c)))]
    if len(printed) < 2:
        raise InputError('The sheet "results" of this workbook holds fewer than two isomers.')
    names = [nm for nm, _ in printed]
    info = {}
    for (r, c), v in res.items():
        if isinstance(v, str) and v.startswith('Theory Level selected'):
            info['Theory level'] = v.split(':', 1)[1].strip()
        if isinstance(v, str) and c == 1 and r == 1:
            info['Mode'] = v.strip()
    mats = {}
    for sheet in ('p_sca', 'p_uns'):
        if sheet in wb.sheets:
            cells = wb.cells(sheet)
            cols = {text(v): c for (r, c), v in cells.items() if r == 1 and text(v) in names}
            if set(cols) == set(names):
                nrows = sorted(set(r for (r, c) in cells if r > 1 and c == cols[names[0]]))
                mats[sheet] = [[cells.get((r, cols[nm])) for nm in names] for r in nrows]
    if len(mats) == 2 and len(mats['p_sca']) == len(mats['p_uns']) and mats['p_sca'] and all(_is_num(x) for m in mats.values() for row in m for x in row):
        h = 0.5e-10                                         # the terms are written with ten decimals
        lo, hi, mid = [0.0] * len(names), [0.0] * len(names), [0.0] * len(names)
        for m in mats.values():
            for row in m:
                for i, x in enumerate(row):
                    lo[i] += _ln(x - h); hi[i] += _ln(min(x + h, 1.0)); mid[i] += _ln(x)
        mol = Molecule(names, lo, hi, mid, shifts=len(mats['p_sca']), source='DP4+App results, sheets "p_sca" and "p_uns"', info=info)
        p = mol.probabilities()
        mol.shown = ['%.6g%%' % (100 * x) for x in p]
        worst = max(abs(a - b) for a, (_, b) in zip(p, printed))
        if worst > 1.5e-4:
            mol.notes.append('The probabilities recomputed from the sheets "p_sca" and "p_uns" differ from the sheet "results" by up to %.2f '
                             'percentage points. Check that the workbook was not edited.' % (100 * worst))
        return mol
    lo = [_ln(v - 0.5e-4) for _, v in printed]
    hi = [_ln(v + 0.5e-4) for _, v in printed]
    mol = Molecule(names, lo, hi, [_ln(v) for _, v in printed], shifts=None, source='DP4+App results, sheet "results" (four decimals)',
                   shown=['%.2f%%' % (100 * v) for _, v in printed], info=info)
    mol.notes.append('The sheets "p_sca" and "p_uns" were not found, so the probabilities come from the sheet "results", which keeps four decimals. '
                     'Give the number of shifts by hand.')
    return mol


def read_workbook(data):
    wb = Workbook(data)
    mol = read_dp4plus_app(wb) or read_dp4plus_spreadsheet(wb)
    if mol is None:
        raise InputError('This workbook is neither the DP4+ spreadsheet nor a results file of DP4+App: no row "DP4+" or "Full" was found. '
                         'Sheets found: %s. Type the probabilities by hand instead.' % ', '.join(wb.sheets))
    return mol


# ---------------------------------------------------------------------------------------------------------------
# reading PyDP4 and DP4-AI output
# ---------------------------------------------------------------------------------------------------------------
_SECTION = re.compile(r'^Assigned (C|H) shifts(?: for isomer (\d+))?:', re.M)
_ISOMER = re.compile(r'^Assigned shifts for isomer (\d+):', re.M)
_ROW = re.compile(r'^([A-Z][a-z]?\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$', re.M)


def _ln_erfc(x):
    return math.log(max(math.erfc(x), 1e-300))


def read_pydp4(text):
    """Output of PyDP4 or DP4-AI. The printed probabilities keep one decimal, so the likelihood is computed again from
    the listed errors of the scaled shifts, with the single Gaussian of PyDP4 (the model of the DP4 workflows of the paper)."""
    marks = [(m.start(), m.group(1), m.group(2)) for m in _SECTION.finditer(text)]
    if not marks:
        return None
    iso_marks = [(m.start(), int(m.group(1))) for m in _ISOMER.finditer(text)]
    isomers = collections.OrderedDict()
    for k, (pos, nuc, iso) in enumerate(marks):
        end = marks[k + 1][0] if k + 1 < len(marks) else len(text)
        stop = text.find('Results of DP4', pos)
        if 0 <= stop < end:
            end = stop
        nxt = [p for p, _ in iso_marks if pos < p < end]
        if nxt:
            end = nxt[0]
        if iso is None:                                     # older layout: "Assigned shifts for isomer N:" stands above both tables
            prev = [n for p, n in iso_marks if p < pos]
            iso = prev[-1] if prev else 1
        rows = [(lab, float(e)) for lab, _, _, _, e in _ROW.findall(text[pos:end])]
        isomers.setdefault(int(iso), {})[nuc] = rows
    if len(isomers) < 2:
        raise InputError('The PyDP4 output lists the shifts of one isomer only; at least two candidates are needed.')
    counts = set((len(d.get('C', [])), len(d.get('H', []))) for d in isomers.values())
    if len(counts) != 1:
        raise InputError('The isomers of the PyDP4 output have different numbers of assigned shifts, so their likelihoods cannot be compared.')
    n_c, n_h = counts.pop()
    if n_c + n_h == 0:
        raise InputError('No assigned shifts were found in the PyDP4 output.')
    h = 0.005                                               # the errors are printed with two decimals
    names, lo, hi, mid = [], [], [], []
    for iso, d in isomers.items():
        a = b = c = 0.0
        for nuc in 'CH':
            z = SIGMA[nuc] * math.sqrt(2.0)
            for _, e in d.get(nuc, []):
                e = abs(e)
                a += _ln_erfc((e + h) / z); b += _ln_erfc(max(e - h, 0.0) / z); c += _ln_erfc(e / z)
        names.append('Isomer %d' % iso); lo.append(a); hi.append(b); mid.append(c)
    info = {}
    for key, val in re.findall(r'^(Solvent|Force Field|DFT NMR Functional|DFT NMR Basis|DFT Opt Functional|DFT Opt Basis|DFT Single Point Functional|'
                               r'DFT Single Point Basis|Stats model)\s*=\s*(.*)$', text, re.M):
        info[key] = val.strip()
    mol = Molecule(names, lo, hi, mid, shifts=n_c + n_h, source='PyDP4 output (%d carbon and %d proton shifts), likelihood recomputed' % (n_c, n_h), info=info)
    p = mol.probabilities()
    mol.shown = ['%.4g%%' % (100 * x) for x in p]
    printed = None
    for title in ('Results of DP4 using all available data', 'Results of DP4'):
        m = re.search(r'^%s:\s*\n((?:Isomer \d+:\s*\S+%%\s*\n?)+)' % re.escape(title), text, re.M)
        if m:
            vals = re.findall(r'Isomer \d+:\s*(\S+)%', m.group(1))
            try:
                printed = [float(v) for v in vals]
            except ValueError:
                printed = None
            if printed is None or any(v != v for v in printed):
                raise InputError('This PyDP4 output prints "nan" for its probabilities: the run did not give a valid result, usually because the '
                                 'experimental shifts were not read or assigned properly. No set is computed from it.')
            break
    if printed and len(printed) == len(p) and all(v == v for v in printed):
        worst = max(abs(100 * a - b) for a, b in zip(p, printed))
        mol.info['Printed probabilities'] = ', '.join('%.1f%%' % v for v in printed)
        if worst > 0.25:
            mol.notes.append('The probabilities printed in the output differ from those recomputed here by up to %.1f percentage points. The output was '
                             'probably made with another error model (for example a multi-Gaussian statistics file) or another rounding. The sets below '
                             'use the single Gaussian of PyDP4, with which the DP4 workflows of the paper were calibrated.' % worst)
    return mol


def read_file(name, data):
    """Any supported file -> Molecule. A table with several molecules is refused here (use it for calibration)."""
    if data[:2] == b'PK':
        return read_workbook(data)
    if data[:8] == b'\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1':
        raise InputError('This is a workbook in the old .xls format. Open it in Excel, save it as .xlsx and read it again.')
    try:
        text = data.decode('utf-8-sig')
    except UnicodeDecodeError:
        text = data.decode('latin-1')
    if '\x00' in text[:2000]:
        raise InputError('This file is not a text file, an .xlsx workbook or a PyDP4 output.')
    mol = read_pydp4(text)
    if mol is not None:
        return mol
    first = text.strip().split('\n')[0].lower() if text.strip() else ''
    if any(k in first for k in ('candidate', 'isomer')) and any(k in first for k in ('molecule', 'case', 'compound')):
        entries = read_table(text)
        if len(entries) != 1:
            raise InputError('This table holds %d molecules. A single molecule can be read here; use the table on the page "Calibrate".' % len(entries))
        entries[0][2].source = 'table, molecule %s' % entries[0][0]
        return entries[0][2]
    return read_typed(text)


# ---------------------------------------------------------------------------------------------------------------
# plain-text report and command line
# ---------------------------------------------------------------------------------------------------------------
def fmt_score(x):
    return 'inf' if x == INF else '%.3f' % x


def report_assign(res):
    out = ['ConformalStereo %s' % __version__, 'Source: %s' % res['source']]
    out.append('Cut-offs from %s' % res['origin'])
    out.append('%d candidates, %d shifts' % (res['candidates'], res['shifts']))
    out.append('')
    levels = res['levels']
    head = '%-22s %14s %10s' % ('candidate', 'probability', 'score') + ''.join('  %6s' % l['label'] for l in levels)
    out += [head, '-' * len(head)]
    mark = {IN: 'in', OUT: '-', OPEN: 'in?'}
    for i, r in enumerate(res['rows']):
        out.append('%-22s %13.6g%% %10s' % (r['name'][:22], 100 * r['probability'], fmt_score(r['score'])) + ''.join('  %6s' % mark[l['status'][i]] for l in levels))
    out.append('')
    for l in levels:
        out.append(l['sentence'])
    for w in res['warnings']:
        out.append('Note: ' + w)
    out += ['', res['caveat']]
    return '\n'.join(out)


def own_origin(cal):
    return 'an own calibration with %d molecules in %d families' % (cal['molecules'], cal['families'])


_BASIS = {'6-31g(d)': '6-31G*', '6-31g(d,p)': '6-31G**', '6-31+g(d,p)': '6-31+G**', '6-311g(d)': '6-311G*', '6-311g(d,p)': '6-311G**',
          '6-311+g(d,p)': '6-311+G**'}


def suggest_workflow(mol, workflows=None):
    """the workflow of the paper that carries the level named in a DP4+ spreadsheet, or None"""
    i = mol.info
    if not ('Functional' in i and 'Basis Set' in i and 'Solvent' in i):
        return None
    basis = _BASIS.get(i['Basis Set'].replace(' ', '').lower())
    solvent = {'pcm': 'PCM', 'gas': 'gas'}.get(i['Solvent'].strip().lower())
    if not basis or not solvent:
        return None
    level = '%s/%s/%s' % (i['Functional'].strip(), basis, solvent)
    for w in (workflows if workflows is not None else load_workflows()):
        if w['benchmark'] == 'JOC2015' and w['level'].lower() == level.lower():
            return w['id']
    return None


def report_calibrate(res):
    out = ['ConformalStereo %s' % __version__, '%d molecules of known structure in %d families; %d to %d shifts, %d to %d candidates; top candidate right in %d'
           % (res['molecules'], res['families'], res['shifts_min'], res['shifts_max'], res['candidates_min'], res['candidates_max'], res['top_right']), '']
    for l in res['levels']:
        if l['supported']:
            extra = ''
            if 'mean_set_size' in l:
                extra = '   (own family left out: actual structure in the set %.0f%%, %.2f structures per set, one structure in %.0f%%)' % (
                    100 * l['coverage_left_out'], l['mean_set_size'], 100 * l['single_structure_sets'])
            out.append('%-4s cut-off per shift %.6f%s' % (l['label'], l['cutoff'], extra))
        else:
            out.append('%-4s not supported: at least %d families are needed' % (l['label'], l['families_needed']))
    for w in res['warnings']:
        out.append('Note: ' + w)
    return '\n'.join(out)


def _cutoffs_from_args(a):
    if a.cutoff:
        cut = {}
        for item in a.cutoff:
            lv, _, val = item.partition('=')
            try:
                cut[float(lv) / 100.0 if float(lv) > 1 else float(lv)] = float(val)
            except ValueError:
                raise InputError('Write a cut-off as level=value, for example 90=0.461.')
        return cut, '', 'the user', None
    if a.calibration:
        with io.open(a.calibration, encoding='utf-8-sig') as f:
            cal = calibrate(read_table(f.read()))
        return {l['level']: l['cutoff'] for l in cal['levels']}, '', own_origin(cal), cal
    if a.workflow:
        w = find_workflow(a.workflow)
        return dict(w['cutoffs']), w['label'], None, w
    raise InputError('Give the cut-off: --workflow (one of the paper), --calibration (your own table) or --cutoff 90=0.46.')


def main(argv=None):
    p = argparse.ArgumentParser(prog='conformalstereo', description='Calibrated candidate sets for DP4-type NMR assignments.')
    sub = p.add_subparsers(dest='command')
    a = sub.add_parser('assign', help='the calibrated sets of one molecule')
    a.add_argument('file', nargs='?', help='DP4+ spreadsheet (.xlsx), DP4+App results (.xlsx), PyDP4 output, or a text file with one "name value" per line')
    a.add_argument('--values', help='instead of a file: "A=99.2,B=0.7,C=0.1"')
    a.add_argument('--loglik', action='store_true', help='the values are natural logarithms of likelihoods, not probabilities')
    a.add_argument('--shifts', type=int, help='number of shifts N (read from the file where possible)')
    a.add_argument('--workflow', help='a workflow of the paper, by words of its name, e.g. "JOC2015 mPW1PW91/6-31+G**/PCM"')
    a.add_argument('--calibration', help='a table of your own molecules of known structure (see "calibrate")')
    a.add_argument('--cutoff', action='append', help='a cut-off per shift as level=value, e.g. --cutoff 90=0.461 (repeatable)')
    a.add_argument('--json', action='store_true', help='print the result as JSON')
    c = sub.add_parser('calibrate', help='cut-offs from your own molecules of known structure')
    c.add_argument('table', help='table with one row per candidate: molecule, family, shifts, candidate, probability or log_likelihood, actual')
    c.add_argument('--json', action='store_true')
    sub.add_parser('workflows', help='list the 76 workflows of the paper with their cut-offs')
    sub.add_parser('selftest', help='check that the program works on this computer')
    p.add_argument('--version', action='version', version=__version__)
    args = p.parse_args(argv)
    try:
        if args.command == 'assign':
            if args.values:
                pairs = [x.partition('=') for x in args.values.split(',') if x.strip()]
                mol = molecule_from_values([(k.strip() if v else 'candidate %d' % (i + 1)) for i, (k, _, v) in enumerate(pairs)],
                                           [(v if v else k) for k, _, v in pairs], kind='loglik' if args.loglik else 'probability')
            elif args.file:
                with open(args.file, 'rb') as f:
                    data = f.read()
                mol = read_file(os.path.basename(args.file), data)
                if args.loglik and mol.source == 'typed values':
                    mol = read_typed(data.decode('utf-8-sig'), kind='loglik')
            else:
                raise InputError('Give a file or --values.')
            cut, label, origin, cal = _cutoffs_from_args(args)
            res = assign(mol, cut, shifts=args.shifts, workflow=label, calibration=cal, origin=origin)
            print(json.dumps(clean(res), indent=1) if args.json else report_assign(res))
        elif args.command == 'calibrate':
            with io.open(args.table, encoding='utf-8-sig') as f:
                res = calibrate(read_table(f.read()))
            print(json.dumps(clean(res), indent=1) if args.json else report_calibrate(res))
        elif args.command == 'workflows':
            print('%-58s %8s %8s %8s %8s' % ('workflow', 'families', 'c 80%', 'c 90%', 'c 95%'))
            for w in load_workflows():
                print('%-58s %8d %8s %8s %8s' % (w['id'][:58], w['families'], *[('%.4f' % w['cutoffs'][l]) if w['cutoffs'][l] is not None else '-' for l in LEVELS]))
        elif args.command == 'selftest':
            return selftest()
        else:
            p.print_help()
    except InputError as e:
        print('ConformalStereo: %s' % e, file=sys.stderr)
        return 2
    return 0


def clean(x):
    """infinities and missing numbers in a form that JSON can carry"""
    if isinstance(x, float):
        if x != x:
            return None
        if x == INF:
            return 'inf'
        if x == -INF:
            return '-inf'
        return x
    if isinstance(x, dict):
        return {(level_text(k) if isinstance(k, float) else k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    return x


def selftest():
    """a few fixed cases; returns 0 if all agree"""
    ok = True

    def check(what, cond):
        nonlocal ok
        print('  %s  %s' % ('ok  ' if cond else 'FAIL', what))
        ok = ok and cond
    m = read_typed('A 99.2\nB 0.7\nC 0.1', shifts=20)
    r = assign(m, {0.9: 0.25})
    check('typed percentages give the scores ln(P_best/P_i)/N', abs(r['rows'][1]['score'] - math.log(99.2 / 0.7) / 20) < 1e-12)
    check('the set at a cut-off of 0.25 holds A and B', r['levels'][0]['members'] == ['A', 'B'])
    check('families needed: 4, 9 and 19', [families_needed(l) for l in LEVELS] == [4, 9, 19])
    items = [(0.1 * i, 'f%d' % i) for i in range(1, 10)]
    check('cut-off of nine families at 90% is the largest score', abs(cutoff(items, 0.9) - 0.9) < 1e-12 and cutoff(items, 0.95) == INF)
    check('rounded values leave a candidate open', assign(read_typed('A 100.0\nB 0.0', shifts=30), {0.9: 0.4})['levels'][0]['open'] == ['B'])
    try:
        w = load_workflows()
        check('76 workflows with cut-offs are on file', len(w) == 76 and all(x['cutoffs'][0.8] is not None for x in w))
    except Exception as e:                                  # noqa
        check('76 workflows with cut-offs are on file (%s)' % e, False)
    print('ConformalStereo %s on Python %s: %s' % (__version__, sys.version.split()[0], 'all checks passed' if ok else 'A CHECK FAILED'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
