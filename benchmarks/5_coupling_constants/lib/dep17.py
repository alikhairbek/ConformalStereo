"""Access to the parsed Cambridge deposits of the Goodman group, and the layout of each.

CS_DEPOSIT selects the deposit (environment variable):
  dp4x2     doi 10.17863/CAM.13222, archive BasisFunc: twelve NMR levels on MMFF geometries (default)
  conf2018  doi 10.17863/CAM.32463, archives <Molecule>Geom.zip: three NMR levels on six kinds of geometry

Small files (descriptions, structures, DP4 outputs) are stored once per content; an index maps every path of the
deposit to the stored copy.
"""
import gzip
import json
import os
import re
import zipfile

KEY = os.environ.get('CS_DEPOSIT', 'dp4x2')
PREFIX = {'dp4x2': 'dep17', 'conf2018': 'dep19'}.get(KEY, KEY)   # prefix of the table names

if KEY == 'dp4x2':
    LEVELS = ['BF%d' % i for i in range(1, 13)]
    LEVEL_NAMES = {'BF1': 'B3LYP/6-31G(d,p)', 'BF2': 'B3LYP/6-311G(d)', 'BF3': 'mPW1PW91/6-31G(d,p)',
                   'BF4': 'mPW1PW91/6-311G(d)', 'BF5': 'WP04/6-31G(d,p)', 'BF6': 'WP04/6-311G(d)',
                   'BF7': 'M06-2X/6-31G(d,p)', 'BF8': 'M06-2X/6-311G(d)', 'BF9': 'B3LYP/gen', 'BF10': 'mPW1PW91/gen',
                   'BF11': 'WP04/gen', 'BF12': 'M06-2X/gen'}
    GEOMETRY = {l: 'MMFF' for l in LEVELS}
    NMR = dict(LEVEL_NAMES)
elif KEY == 'conf2018':
    _NMR = {'BF1': 'B3LYP/6-31G(d,p)', 'BF4': 'mPW1PW91/6-311G(d)', 'BF7': 'M06-2X/6-31G(d,p)'}
    _GEO = [('', 'MMFF'), ('s', 'PM6'), ('pm7', 'PM7'), ('h', 'HF/6-31G(d,p)'), ('d', 'B3LYP/6-31G(d,p)'), ('m', 'M06-2X/6-31G(d,p)')]
    LEVELS = [n + g for g, _ in _GEO for n in ('BF1', 'BF4', 'BF7')]
    GEOMETRY = {n + g: name for g, name in _GEO for n in _NMR}
    NMR = {n + g: _NMR[n] for g, _ in _GEO for n in _NMR}
    LEVEL_NAMES = {l: '%s // %s' % (NMR[l], GEOMETRY[l]) for l in LEVELS}
else:
    raise SystemExit('unknown deposit %r' % KEY)


def level_index(level):
    return LEVELS.index(level) if level in LEVELS else len(LEVELS)


def gaussian_files(core):
    return [os.path.join(core, fn) for fn in sorted(os.listdir(core)) if re.match(r'gaussian_%s_\d+\.jsonl\.gz$' % KEY, fn)]


def level_of(folder):
    """level folder -> (name of the calculation folder without its level suffix, level)"""
    top = folder.split('/')[0]
    if KEY == 'dp4x2':
        m = re.search(r'_(BF\d+|J)$', top)
        if m:
            return top[:m.start()], m.group(1)
        return top, 'BF1'
    m = re.search(r'_BF(\d)([a-z0-9]*)$', top)
    if m:
        return top[:m.start()], 'BF' + m.group(1) + m.group(2)
    return top, 'BF1'


class Small(object):
    def __init__(self, core):
        self.zip = zipfile.ZipFile(core + '/small_%s_001.zip' % KEY)
        self.index = {}
        with gzip.open(core + '/small_index_%s.tsv.gz' % KEY, 'rt') as f:
            next(f)
            for line in f:
                path, size, md5, stored, same = line.rstrip('\n').split('\t')
                self.index[path] = (int(size), md5, same or path)

    def paths(self):
        return sorted(self.index)

    def read(self, path):
        src = self.index[path][2]
        return self.zip.read(src.replace('!/', '/')).decode('utf-8', 'replace')


def split(path):
    """deposit path -> (molecule archive, folder, file name) or None"""
    p = path.split('!/')
    if KEY == 'dp4x2':
        if len(p) != 3:
            return None
        arch, inner = p[1].split('/')[-1].replace('.zip', ''), p[2]
    else:
        if len(p) != 2:
            return None
        arch, inner = p[0].replace('Geom.zip', ''), p[1]
    parts = inner.split('/')
    if len(parts) < 2:
        return arch, '', parts[-1]
    return arch, '/'.join(parts[:-1]), parts[-1]


def jsonl(path):
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.strip():
                yield json.loads(line)
