"""Reference tools for the coupling-constant tests.

read_cmatrix: the algorithm of PyDP4's ReadCMatrix (nmrPredictGaus.py, 2017 version, the reader the depositors used),
              carried over to Python 3. It is independent of the notebook's reader: it counts the rows of the first
              block and fills a square matrix by position, where the notebook keys every value by its row and column
              numbers.
write_triangle: a symmetric matrix as Gaussian prints it (lower triangle, five columns at a time, D exponents).
mixed_output: a two-step output laid out like the coupling jobs of the 2017 deposit (NMR=(SpinSpin,Mixed,ReadAtoms)).
"""
import random


def read_cmatrix(text, title):
    lines = [ln + '\n' for ln in text.split('\n')]
    toks = lambda ln: [t for t in ln.split(' ') if t]
    index = 0
    while title not in lines[index]:
        index += 1
        if index == len(lines):
            return [0], ['']
    start, end, started = index, 0, False
    for i in range(index, len(lines)):
        data = toks(lines[i])
        if 'D' in data[-1]:
            if not started:
                started, start = True, i
        elif started:
            end = i
            break
    natoms = end - start
    labels = []
    matrix = [[0.0] * natoms for _ in range(natoms)]
    index = start - 1
    coln = rown = 0
    while toks(lines[index])[-2].isnumeric():
        index += 1
        data = toks(lines[index])
        while 'D' in data[-1]:
            if len(labels) < natoms:
                labels.append(data[0])
            vals = [float(x.replace('D', 'E')) for x in data[1:]]
            matrix[rown][coln:coln + len(vals)] = vals
            rown += 1
            index += 1
            data = toks(lines[index])
        coln += 5
        rown = coln
    for x in range(natoms):
        matrix[x][x] = 0.0
    for x in range(natoms):
        for y in range(x, natoms):
            matrix[x][y] = matrix[y][x]
    return matrix, labels


def dfmt(x):
    """Fortran D14.6 without the padding: 0.926067D+01"""
    if x == 0:
        return '0.000000D+00'
    mant, exp = ('%.5E' % x).split('E')
    return '%s0.%sD%+03d' % ('-' if mant.startswith('-') else '', mant.lstrip('-').replace('.', ''), int(exp) + 1)


def printed(x):
    return float(dfmt(x).replace('D', 'E'))


def write_triangle(labels, m):
    out = []
    n = len(labels)
    for c0 in range(0, n, 5):
        cols = list(range(c0, min(c0 + 5, n)))
        out.append('%17d' % labels[cols[0]] + ''.join('%14d' % labels[c] for c in cols[1:]) + ' ')
        for r in range(c0, n):
            out.append('%7d' % labels[r] + ''.join('%14s' % dfmt(m[r][c]) for c in cols if c <= r))
    return '\n'.join(out) + '\n'


def random_matrix(n, rng, scale=12.0):
    m = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i):
            m[i][j] = m[j][i] = printed(rng.choice((1, 1, 1, 0.01, 30)) * rng.uniform(-scale, scale))
    return m


TITLES = (('FC', 'Fermi Contact (FC) contribution to %s (Hz):'), ('SD', 'Spin-dipolar (SD) contribution to %s (Hz):'),
          ('PSO', 'Paramagnetic spin-orbit (PSO) contribution to %s (Hz):'), ('DSO', 'Diamagnetic spin-orbit (DSO) contribution to %s (Hz):'),
          ('total', 'Total nuclear spin-spin coupling %s (Hz): '))
ROUTE0 = '# b3lyp/6-31g(d,p) nmr=(giao,spinspin,mixed,readatoms) scrf=(solvent=chloroform)'
ROUTE1 = '#N Geom=AllCheck Guess=TCheck SCRF=Check GenChk RB3LYP/ChkBas NMR=ReadFC IOp(3/6=22)'
RULE = ' ' + '-' * 70


def mixed_output(real_log, labels, seed, scf=-382.06, crlf=False, damage=None, cut=False):
    """A two-step coupling job built around the text of a real Gaussian NMR log.
    Step 1 (uncontracted basis): Fermi-contact matrices only. Step 2: every term and the total, the Fermi-contact part repeated.
    Returns (text, {(step, name): matrix}) for the J matrices; `damage` removes one row of the total matrix of step 2;
    `cut` ends the file before any matrix is printed."""
    rng = random.Random(seed)
    n = len(labels)
    truth = {}
    head, _, rest = real_log.partition(' Fermi Contact (FC) contribution to K (Hz):')
    assert rest, 'the template log prints no coupling matrices'
    a = head.index('\n #p ')
    b = head.index('\n', a + 1)
    step0 = head[:a] + '\n ' + ROUTE0 + head[b:]
    if cut:
        return step0 + ' Error termination via Lnk1e in /apps/gaussian/g09/l1002.exe\n', truth
    fc = random_matrix(n, rng)
    text = step0
    for kind in 'KJ':
        text += ' ' + TITLES[0][1] % kind + '\n' + write_triangle(labels, fc if kind == 'J' else random_matrix(n, rng, 2.0))
    truth[(0, 'FC')] = fc
    text += ' End of Minotr F.D. properties file   721 does not exist.\n Job cpu time:       0 days  1 hours  1 minutes  1.0 seconds.\n'
    text += ' Normal termination of Gaussian 09 at Thu Dec 10 08:19:44 2015.\n'
    shield = real_log[real_log.index(' SCF GIAO Magnetic shielding tensor (ppm):'):real_log.index(' Fermi Contact (FC) contribution to K (Hz):')]
    text += RULE + '\n ' + ROUTE1 + '\n' + RULE + '\n SCF Done:  E(RB3LYP) =  %.9f     A.U. after   10 cycles\n' % scf + shield
    parts = {name: random_matrix(n, rng, 0.6) for name, _ in TITLES[1:4]}
    total = [[printed(fc[i][j] + sum(parts[k][i][j] for k in parts)) for j in range(n)] for i in range(n)]
    for name, title in TITLES:
        m = fc if name == 'FC' else total if name == 'total' else parts[name]
        for kind in 'KJ':
            block = write_triangle(labels, m if kind == 'J' else random_matrix(n, rng, 2.0))
            if damage and name == 'total' and kind == 'J':
                rows = block.split('\n')
                del rows[damage]
                block = '\n'.join(rows)
            text += ' ' + title % kind + '\n' + block
    truth[(1, 'FC')], truth[(1, 'total')] = fc, total
    text += ' End of Minotr F.D. properties file   721 does not exist.\n End of Minotr F.D. properties file   722 does not exist.\n'
    text += ' Job cpu time:       0 days  0 hours 30 minutes  2.2 seconds.\n Normal termination of Gaussian 09 at Thu Dec 10 08:19:44 2015.\n'
    if crlf:
        text = text.replace('\n', '\r\n')
    return text, truth
