"""Tests of ConformalStereo. They need Python only:

    python tests/test_conformalstereo.py

Tests that use the pooled table or the tables of the analysis are skipped when the program folder stands alone."""
import os, sys, io, csv, json, math, time, base64, unittest, threading, collections
import urllib.request, urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
sys.path.insert(0, APP)
import conformalstereo as cs
import run as server

FIX = os.path.join(HERE, 'fixtures')
REPO = os.path.dirname(os.path.dirname(APP))
POOLED = os.path.join(REPO, 'data', 'pooled_candidates.csv')
WORKFLOWS = os.path.join(REPO, 'analysis', 'tables', 'final_workflows.csv')
RECORDS = os.path.join(REPO, 'analysis', 'tables', 'final_records.csv')


def data(name):
    with open(os.path.join(FIX, name), 'rb') as f:
        return f.read()


def text_of(*parts):
    with io.open(os.path.join(*parts), encoding='utf-8') as f:
        return f.read()


def pooled():
    groups = collections.OrderedDict()
    with io.open(POOLED, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            groups.setdefault((r['benchmark'], r['method'], r['geometry'], r['level']), collections.OrderedDict()).setdefault(r['case'], []).append(r)
    return groups


class Numbers(unittest.TestCase):
    def test_precision_follows_the_digits(self):
        self.assertEqual(cs.parse_number('99.87'), (99.87, 0.005))
        self.assertEqual(cs.parse_number('99,87 %'), (99.87, 0.005))
        self.assertEqual(cs.parse_number('100'), (100.0, 0.5))
        v, h = cs.parse_number('1.2e-5'); self.assertAlmostEqual(v, 1.2e-5); self.assertAlmostEqual(h, 0.05e-5)
        v, h = cs.parse_number('−3.25'); self.assertEqual(v, -3.25)
        self.assertEqual(cs.parse_number('-inf')[0], -cs.INF)
        for bad in ('abc', '', '1.2.3', '--1'):
            with self.assertRaises(cs.InputError):
                cs.parse_number(bad)

    def test_typed_lines(self):
        m = cs.read_typed('candidate, probability\nisomer 1, 99.2\nisomer 2;0.7%\nC\t0,1\n')
        self.assertEqual(m.names, ['isomer 1', 'isomer 2', 'C'])
        self.assertAlmostEqual(m.probabilities()[0], 0.992)
        m = cs.read_typed('0.9\n0.1')
        self.assertEqual(m.names, ['candidate 1', 'candidate 2'])
        m = cs.read_typed('a -3.5\nb -10.25', kind='loglik')
        self.assertAlmostEqual(m.mid[1] - m.mid[0], -6.75)
        for bad in ('only one 5', 'a 1\nb x', 'a -1\nb 2', 'a 0\nb 0.0', 'a 1\na 2'):
            with self.assertRaises(cs.InputError):
                cs.read_typed(bad)


class Method(unittest.TestCase):
    def test_scores(self):
        m = cs.read_typed('A 0.9000000000\nB 0.0900000000\nC 0.0100000000', shifts=10)
        b = cs.score_bounds(m)
        self.assertAlmostEqual(b[0][1], 0.0)
        self.assertAlmostEqual(b[1][1], math.log(10) / 10, places=9)
        self.assertAlmostEqual(b[2][1], math.log(90) / 10, places=9)
        self.assertTrue(m.exact)
        with self.assertRaises(cs.InputError):
            cs.score_bounds(cs.read_typed('A 1\nB 2'))

    def test_scores_do_not_depend_on_the_scale(self):
        a = cs.score_bounds(cs.read_typed('A 99.20000000\nB 0.70000000\nC 0.10000000', shifts=7))
        b = cs.score_bounds(cs.read_typed('A 0.9920000000\nB 0.0070000000\nC 0.0010000000', shifts=7))
        for x, y in zip(a, b):
            self.assertAlmostEqual(x[1], y[1], places=8)

    def test_rounded_values(self):
        m = cs.read_typed('A 100.0\nB 0.0', shifts=30)
        b = cs.score_bounds(m)
        self.assertEqual(b[1][2], cs.INF)
        self.assertAlmostEqual(b[1][0], math.log(99.95 / 0.05) / 30)
        self.assertFalse(m.exact)
        self.assertEqual(cs.conformal_set(b, 0.2), [cs.IN, cs.OUT])          # the lower bound alone excludes B
        self.assertEqual(cs.conformal_set(b, 0.3), [cs.IN, cs.OPEN])         # here the digits do not decide
        r = cs.assign(m, {0.9: 0.3})
        self.assertEqual(r['levels'][0]['members'], ['A', 'B'])
        self.assertEqual(r['levels'][0]['open'], ['B'])
        self.assertTrue(any('rounded' in w for w in r['warnings']))

    def test_exact_zero_is_outside(self):
        m = cs.Molecule(['A', 'B'], [0.0, -cs.INF], [0.0, -cs.INF], shifts=5)
        self.assertEqual(cs.conformal_set(cs.score_bounds(m), 99.0), [cs.IN, cs.OUT])

    def test_cutoff(self):
        self.assertEqual([cs.families_needed(l) for l in cs.LEVELS], [4, 9, 19])
        items = [(0.1 * i, 'f%d' % i) for i in range(1, 10)]
        self.assertAlmostEqual(cs.cutoff(items, 0.9), 0.9)
        self.assertEqual(cs.cutoff(items, 0.95), cs.INF)
        self.assertAlmostEqual(cs.cutoff(items, 0.8), 0.8)
        self.assertEqual(cs.cutoff(items[:3], 0.8), cs.INF)
        two = [(0.5, 'a'), (0.1, 'a'), (0.2, 'b'), (0.3, 'c'), (0.4, 'd')]    # family a counts once, split over its two molecules
        self.assertAlmostEqual(cs.cutoff(two, 0.8), 0.5)
        self.assertAlmostEqual(cs.cutoff(two, 0.7), 0.4)

    def test_cutoff_is_the_weighted_quantile(self):
        import random
        rnd = random.Random(3)
        for _ in range(200):
            g = rnd.randint(1, 25)
            items = [(rnd.random(), 'f%d' % rnd.randint(1, g)) for _ in range(rnd.randint(1, 60))]
            fams = collections.Counter(f for _, f in items)
            for level in cs.LEVELS:
                c = cs.cutoff(items, level)
                def mass(x):
                    return sum(1.0 / fams[f] / (len(fams) + 1) for s, f in items if s <= x)
                if c == cs.INF:
                    self.assertLess(mass(1e9), level - 1e-12)
                else:
                    self.assertGreaterEqual(mass(c), level - 1e-9)
                    below = [s for s, _ in items if s < c]
                    if below:
                        self.assertLess(mass(max(below)), level - 1e-12)

    def test_sets_are_nested_and_keep_the_best(self):
        m = cs.read_typed('A 60.00000\nB 30.00000\nC 9.90000\nD 0.10000', shifts=12)
        r = cs.assign(m, {0.8: 0.05, 0.9: 0.16, 0.95: 0.6})
        sizes = [l['size'] for l in r['levels']]
        self.assertEqual(sizes, sorted(sizes))
        self.assertTrue(all('A' in l['members'] for l in r['levels']))
        self.assertEqual(r['levels'][0]['members'], ['A'])
        r0 = cs.assign(m, {0.9: 0.0})
        self.assertEqual(r0['levels'][0]['members'], ['A'])
        rn = cs.assign(m, {0.95: None})
        self.assertFalse(rn['levels'][0]['supported']); self.assertEqual(rn['levels'][0]['size'], 4)

    def test_warnings(self):
        m = cs.read_typed('\n'.join('c%d %s' % (i, v) for i, v in enumerate(['50.000', '20.000', '10.000', '10.000', '5.000', '5.000'])), shifts=40)
        r = cs.assign(m, {0.8: 0.1}, calibration=dict(shifts_max=30, candidates_max=4))
        text = ' '.join(r['warnings'])
        self.assertIn('more shifts', text); self.assertIn('more candidates', text); self.assertIn('more than four candidates', text)

    def test_json(self):
        m = cs.read_typed('A 100.0\nB 0.0', shifts=30)
        s = json.dumps(cs.clean(cs.assign(m, {0.9: 0.3, 0.95: None})), allow_nan=False)
        self.assertIn('"inf"', s)


class Tables(unittest.TestCase):
    TABLE = ('molecule,family,shifts,candidate,probability,actual\n' +
             ''.join('m%d,f%d,20,a,%s,1\nm%d,f%d,20,b,%s,0\n' % (i, i, p, i, i, q) for i, (p, q) in enumerate(
                 [('99.0000', '1.0000'), ('90.0000', '10.0000'), ('60.0000', '40.0000'), ('30.0000', '70.0000'), ('99.9000', '0.1000')])))

    def test_calibrate(self):
        res = cs.calibrate(cs.read_table(self.TABLE))
        self.assertEqual((res['molecules'], res['families'], res['top_right']), (5, 5, 4))
        lv = {l['label']: l for l in res['levels']}
        self.assertTrue(lv['80%']['supported']); self.assertFalse(lv['90%']['supported'])
        self.assertAlmostEqual(lv['80%']['cutoff'], math.log(70.0 / 30.0) / 20, places=4)   # 5/6 of the weight: the largest score
        self.assertEqual(lv['90%']['families_needed'], 9)

    def test_families_count_once(self):
        t = self.TABLE.replace('m1,f1', 'm1,f0').replace('m2,f2', 'm2,f0')
        res = cs.calibrate(cs.read_table(t))
        self.assertEqual(res['families'], 3)
        self.assertFalse(res['levels'][0]['supported'])

    def test_other_separators_and_names(self):
        t = 'case;n_nuclei;isomer;log_ratio;is_true\nx;10;1;0.0;1\nx;10;2;-4.5;0\ny;12;1;-1.25;1\ny;12;2;0.0;0\n'
        e = cs.read_table(t)
        self.assertEqual([(n, f, a) for n, f, _, a in e], [('x', 'x', 0), ('y', 'y', 0)])
        self.assertAlmostEqual(cs.score_bounds(e[1][2])[0][1], 1.25 / 12, places=6)

    def test_errors(self):
        for bad in ('', 'molecule,shifts,candidate\nm,5,a\n', 'molecule,shifts,candidate,probability\nm,5,a,1\nm,6,b,2\n',
                    'molecule,shifts,candidate,probability,actual\nm,5,a,1,1\nm,5,b,2,1\n', 'molecule,shifts,candidate,probability,actual\nm,5,a,1,maybe\n',
                    'just one column\nx\n'):
            with self.assertRaises(cs.InputError):
                cs.read_table(bad)
        with self.assertRaises(cs.InputError):                               # no actual structure marked
            cs.calibrate(cs.read_table('molecule,shifts,candidate,probability\nm,5,a,1\nm,5,b,2\n'))


class Readers(unittest.TestCase):
    def test_dp4plus_sheet(self):
        m = cs.read_file('sheet.xlsx', data('dp4plus_sheet.xlsx'))
        self.assertEqual(m.names, ['Isomer 1', 'Isomer 2', 'Isomer 3', 'Isomer 4'])
        self.assertEqual(m.shifts, 18)
        self.assertTrue(m.exact)
        self.assertAlmostEqual(m.mid[1] - m.mid[0], math.log(0.1379999 / 0.862), places=12)
        self.assertEqual(m.mid[3], -cs.INF)
        self.assertEqual(m.info['Functional'], 'B3LYP')
        self.assertEqual(cs.suggest_workflow(m), 'JOC2015 | DP4+ | B3LYP | B3LYP/6-31G*/gas')
        r = cs.assign(m, {0.9: 0.4})
        self.assertEqual(r['levels'][0]['members'], ['Isomer 1', 'Isomer 2'])
        self.assertEqual(r['levels'][0]['status'][3], cs.OUT)

    def test_dp4plus_app(self):
        m = cs.read_file('res.xlsx', data('dp4plus_app_results.xlsx'))
        self.assertEqual(m.shifts, 25)
        expected = [float(l.split()[1]) for l in text_of(FIX, 'dp4plus_app_results.expected.txt').split('\n') if l.strip()]
        for a, b in zip(m.probabilities(), expected):
            self.assertAlmostEqual(a / b, 1.0, places=5)                     # the unrounded values, far below the four decimals of the sheet "results"
        self.assertLess(expected[2], 1e-10)
        self.assertEqual(m.notes, [])
        self.assertEqual(m.info.get('Theory level'), 'B3LYP/6-31G(d) (GAS)')
        b = cs.score_bounds(m)
        self.assertLess(b[1][2] - b[1][0], 1e-6)
        self.assertLess(b[3][2] - b[3][0], 0.05)                             # terms near 1e-9 keep two digits at ten decimals: a small, known range

    def test_pydp4_both_layouts(self):
        expected = [float(x) for x in text_of(FIX, 'pydp4.expected.txt').split()]
        for name in ('pydp4_new.dp4', 'pydp4_old.dp4'):
            m = cs.read_file(name, data(name))
            self.assertEqual((m.names, m.shifts), (['Isomer 1', 'Isomer 2', 'Isomer 3'], 16))
            for i in range(3):
                self.assertAlmostEqual(m.mid[i], expected[i], places=9)
                self.assertLessEqual(m.lo[i], m.mid[i]); self.assertGreaterEqual(m.hi[i], m.mid[i])
            self.assertEqual(m.notes, [])
            self.assertIn('Printed probabilities', m.info)

    def test_pydp4_nan_is_refused(self):
        text = data('pydp4_old.dp4').decode().replace('Isomer 1: ', 'Isomer 1:  nan%\nIsomer 9: ', 1)
        text = text.replace('Results of DP4 using all available data:\nIsomer 1:  nan%\nIsomer 9: ', 'Results of DP4 using all available data:\nIsomer 1:  nan%\nIsomer 2:  nan%\n\nJunk: ', 1)
        if 'nan%' in text:
            with self.assertRaises(cs.InputError):
                cs.read_file('x.dp4', text.encode())

    def test_other_files(self):
        with self.assertRaises(cs.InputError):
            cs.read_file('x.xls', b'\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1' + b'0' * 100)
        with self.assertRaises(cs.InputError):
            cs.read_file('x.xlsx', b'PK\x03\x04 not a workbook')
        with self.assertRaises(cs.InputError):
            cs.read_file('x.bin', b'\x00\x01\x02\x03' * 50)
        m = cs.read_file('values.txt', b'A 70.00\nB 30.00\n')
        self.assertEqual(m.names, ['A', 'B'])
        one = b'molecule,family,shifts,candidate,probability,actual\nm,f,9,a,80.00,1\nm,f,9,b,20.00,0\n'
        self.assertEqual(cs.read_file('t.csv', one).shifts, 9)


class Paper(unittest.TestCase):
    """the numbers of the paper, from the pooled table of the repository"""

    @unittest.skipUnless(os.path.exists(POOLED), 'the pooled table is not here')
    def test_the_76_cutoffs(self):
        on_file = {w['id']: w for w in cs.load_workflows()}
        groups = pooled()
        self.assertEqual(len(groups), 76)
        self.assertEqual(sum(len(c) for c in groups.values()), 3100)
        n = 0
        for key, cases in groups.items():
            items = []
            for case, rr in cases.items():
                ell = [float(r['log_ratio']) for r in rr]
                k = [int(r['is_true']) for r in rr].index(1)
                items.append(((max(ell) - ell[k]) / int(rr[0]['n_nuclei']), rr[0]['family']))
            w = on_file[' | '.join(key)]
            for level in cs.LEVELS:
                c = cs.cutoff(items, level)
                if w['cutoffs'][level] is None:
                    self.assertEqual(c, cs.INF)
                else:
                    self.assertAlmostEqual(c, w['cutoffs'][level], places=8); n += 1
        self.assertEqual(n, 174)

    @unittest.skipUnless(os.path.exists(POOLED) and os.path.exists(WORKFLOWS), 'the tables of the analysis are not here')
    def test_calibration_as_in_the_analysis(self):
        ref = {(r['benchmark'], r['method'], r['geometry'], r['level']): r for r in csv.DictReader(io.StringIO(text_of(WORKFLOWS)))}
        for key, cases in pooled().items():
            lines = ['molecule,family,shifts,candidate,log_likelihood,actual']
            for case, rr in cases.items():
                lines += ['%s,%s,%s,%s,%s,%s' % (case.replace(',', ' '), r['family'], r['n_nuclei'], r['candidate'].replace(',', ' '), r['log_ratio'], r['is_true']) for r in rr]
            res = cs.calibrate(cs.read_table('\n'.join(lines)))
            q = ref[key]
            self.assertEqual((res['molecules'], res['families']), (int(q['records']), int(q['families'])))
            for l in res['levels']:
                tag = '%d' % round(100 * l['level'])
                if q['cutoff_per_shift_' + tag] == '':
                    self.assertFalse(l['supported'])
                    continue
                self.assertAlmostEqual(l['cutoff'], float(q['cutoff_per_shift_' + tag]), places=5)
                self.assertAlmostEqual(l['coverage_left_out'], float(q['coverage_' + tag]), places=6)
                self.assertAlmostEqual(l['mean_set_size'], float(q['mean_set_size_' + tag]), places=6)
                self.assertAlmostEqual(l['single_structure_sets'], float(q['single_structure_sets_' + tag]), places=6)

    @unittest.skipUnless(os.path.exists(POOLED), 'the pooled table is not here')
    def test_worked_example_of_figure_1(self):
        key = ('JOC2015', 'DP4+', 'B3LYP', 'B3LYP/6-31G*/gas')
        cases = pooled()[key]
        fam = cases['87a'][0]['family']
        items = []
        for case, rr in cases.items():
            if rr[0]['family'] == fam:
                continue                                                     # the family of the record is left out
            ell = [float(r['log_ratio']) for r in rr]
            k = [int(r['is_true']) for r in rr].index(1)
            items.append(((max(ell) - ell[k]) / int(rr[0]['n_nuclei']), rr[0]['family']))
        c90, c80 = cs.cutoff(items, 0.9), cs.cutoff(items, 0.8)
        self.assertEqual(('%.2f' % c90, '%.2f' % c80), ('0.51', '0.21'))
        m = cs.read_file('example.txt', text_of(APP, 'examples', 'typed_example.txt').encode('utf-8'))
        r = cs.assign(m, {0.8: c80, 0.9: c90}, shifts=23)
        self.assertEqual(sorted('%.2f' % x['score'] for x in r['rows']), ['0.00', '0.37', '1.70', '1.93'])
        self.assertEqual([l['members'] for l in r['levels']], [['87b'], ['87a', '87b']])
        self.assertEqual(int(cases['87a'][0]['n_nuclei']), 23)
        self.assertEqual([x['candidate'] for x in cases['87a'] if x['is_true'] == '1'], ['87a'])   # the actual structure is in the 90% set, not on top


class Page(unittest.TestCase):
    """the local page: the same answers through the address of this computer"""

    @classmethod
    def setUpClass(cls):
        server.WORKFLOWS = cs.load_workflows()
        cls.httpd = server.ThreadingHTTPServer(('127.0.0.1', 0), server.Handler)
        cls.url = 'http://127.0.0.1:%d' % cls.httpd.server_address[1]
        threading.Thread(target=cls.httpd.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown(); cls.httpd.server_close()

    def call(self, path, body=None, headers=None):
        req = urllib.request.Request(self.url + path, data=None if body is None else json.dumps(body).encode(), headers=dict({'Content-Type': 'application/json'}, **(headers or {})))
        try:
            with urllib.request.urlopen(req, timeout=20) as r:
                return r.status, json.loads(r.read().decode()) if 'json' in r.headers.get('Content-Type', '') else r.read()
        except urllib.error.HTTPError as e:
            return e.code, json.loads(e.read().decode())

    def test_page_and_lists(self):
        code, html = self.call('/')
        self.assertEqual(code, 200); self.assertIn(b'ConformalStereo', html)
        code, d = self.call('/api/workflows')
        self.assertEqual((code, len(d['workflows'])), (200, 76))
        self.assertIn('90%', d['workflows'][0]['cutoffs'])
        code, ex = self.call('/api/examples')
        self.assertTrue(ex['calibration'].startswith('molecule,family')); self.assertIn('87a', ex['typed'])

    def test_assign_typed_and_file(self):
        wid = 'JOC2015 | DP4+ | B3LYP | B3LYP/6-31G*/gas'
        code, r = self.call('/api/assign', dict(source=dict(type='typed', text='87a 0.02184270534\n87b 99.97815729\n87c 4.7e-18\n87d 1.0e-15', kind='probability'), shifts='23', cutoff=dict(type='workflow', id=wid)))
        self.assertEqual(code, 200)
        self.assertEqual([l['members'] for l in r['levels'][:2]], [['87b'], ['87a', '87b']])
        self.assertFalse(r['levels'][2]['supported'])
        self.assertIn('At the 90% level', r['report'])
        b64 = base64.b64encode(data('dp4plus_sheet.xlsx')).decode()
        code, p = self.call('/api/read', dict(source=dict(type='file', name='sheet.xlsx', data=b64)))
        self.assertEqual((code, p['shifts'], p['suggestion']), (200, 18, wid))
        self.assertEqual(p['log_ratio'][3], '-inf')
        code, r = self.call('/api/assign', dict(source=dict(type='file', name='sheet.xlsx', data=b64), shifts='', cutoff=dict(type='manual', cutoffs={'0.9': '0,4', '0.8': ''})))
        self.assertEqual((code, r['levels'][0]['members']), (200, ['Isomer 1', 'Isomer 2']))

    def test_calibrate_then_assign(self):
        code, ex = self.call('/api/examples')
        code, cal = self.call('/api/calibrate', dict(text=ex['calibration']))
        self.assertEqual((code, cal['molecules'], cal['families']), (200, 48, 18))
        self.assertAlmostEqual(cal['levels'][1]['cutoff'], 0.46062125, places=5)
        own = dict(type='own', calibration=dict(molecules=cal['molecules'], families=cal['families'], shifts_max=cal['shifts_max'], candidates_max=cal['candidates_max'],
                                                levels=[dict(level=l['level'], cutoff=l['cutoff']) for l in cal['levels']]))
        code, r = self.call('/api/assign', dict(source=dict(type='typed', text='A 99.20000\nB 0.70000\nC 0.10000', kind='probability'), shifts=40, cutoff=own))
        self.assertEqual(code, 200)
        self.assertIn('own calibration with 48 molecules in 18 families', r['origin'])
        self.assertTrue(any('more shifts' in w for w in r['warnings']))

    def test_messages_for_the_user(self):
        for body in (dict(source=dict(type='typed', text='A 1', kind='probability'), shifts='5', cutoff=dict(type='manual', cutoffs={'0.9': '0.3'})),
                     dict(source=dict(type='typed', text='A 60\nB 40', kind='probability'), shifts='', cutoff=dict(type='manual', cutoffs={'0.9': '0.3'})),
                     dict(source=dict(type='typed', text='A 60\nB 40', kind='probability'), shifts='2.5', cutoff=dict(type='manual', cutoffs={'0.9': '0.3'})),
                     dict(source=dict(type='typed', text='A 60\nB 40', kind='probability'), shifts='5', cutoff=dict(type='manual', cutoffs={'0.9': 'x'})),
                     dict(source=dict(type='typed', text='A 60\nB 40', kind='probability'), shifts='5', cutoff=dict(type='own', calibration=None)),
                     dict(source=dict(type='file', name='a.xlsx', data=base64.b64encode(b'PK nothing').decode()), shifts='5', cutoff=dict(type='manual', cutoffs={'0.9': '0.3'}))):
            code, r = self.call('/api/assign', body)
            self.assertEqual(code, 400, r)
            self.assertTrue(r['error'] and 'Traceback' not in r['error'])
        self.assertEqual(self.call('/api/calibrate', dict(text='nonsense'))[0], 400)
        self.assertEqual(self.call('/api/nothing', {})[0], 404)

    def test_other_hosts_are_refused(self):
        code, r = self.call('/api/workflows', headers={'Host': 'example.com'})
        self.assertEqual(code, 403)


if __name__ == '__main__':
    unittest.main(verbosity=1)
