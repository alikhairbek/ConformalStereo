#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Starts ConformalStereo as a page in the web browser of this computer.

    python run.py            opens the page; keep the window open while you work
    python run.py --help     options

The page talks only to this program on this computer (address 127.0.0.1). Nothing is sent anywhere else."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:                    # the embeddable Python for Windows does not add the folder of the script itself
    sys.path.insert(0, HERE)
import io, json, base64, time, threading, argparse, traceback, webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import conformalstereo as cs

MAX_BODY = 40 * 1024 * 1024
STATE = dict(last=time.time(), page_seen=False, quiet_exit=180.0)


def read_source(src, shifts=None):
    """the 'source' object of a request -> Molecule"""
    kind = src.get('type')
    if kind == 'file':
        try:
            data = base64.b64decode(src.get('data', ''), validate=False)
        except Exception:
            raise cs.InputError('The file did not arrive in one piece. Open it again.')
        if not data:
            raise cs.InputError('The file is empty.')
        mol = cs.read_file(src.get('name', 'file'), data)
    elif kind == 'typed':
        mol = cs.read_typed(src.get('text', ''), kind='loglik' if src.get('kind') == 'loglik' else 'probability')
    else:
        raise cs.InputError('No DP4 result was given yet.')
    if shifts:
        mol.shifts = int(shifts)
    return mol


def preview(mol):
    return dict(names=mol.names, shown=mol.shown, probabilities=mol.probabilities(), shifts=mol.shifts, source=mol.source, info=mol.info,
                notes=mol.notes, exact=mol.exact, suggestion=cs.suggest_workflow(mol, WORKFLOWS), log_ratio=[x - max(mol.mid) if x > -cs.INF else -cs.INF for x in mol.mid])


def number(x, what):
    try:
        v = float(str(x).replace(',', '.'))
    except (TypeError, ValueError):
        raise cs.InputError('%s is not a number.' % what)
    if v != v or v < 0:
        raise cs.InputError('%s must be zero or more.' % what)
    return v


def api(path, body):
    if path == '/api/read':
        return preview(read_source(body.get('source') or {}))
    if path == '/api/assign':
        shifts = body.get('shifts')
        if shifts in ('', None):
            shifts = None
        else:
            shifts = number(shifts, 'The number of shifts')
            if shifts != int(shifts) or shifts < 1:
                raise cs.InputError('The number of shifts must be a whole number above zero.')
        mol = read_source(body.get('source') or {}, shifts)
        if not mol.shifts:
            raise cs.InputError('Give the number of shifts N (carbon plus proton shifts that entered the DP4 calculation).')
        cut = body.get('cutoff') or {}
        calibration, origin, label = None, None, ''
        if cut.get('type') == 'workflow':
            w = cs.find_workflow(cut.get('id', ''), WORKFLOWS)
            cutoffs, label, calibration = dict(w['cutoffs']), w['label'], w
        elif cut.get('type') == 'own':
            cal = cut.get('calibration') or {}
            cutoffs = {l['level']: (l['cutoff'] if l.get('cutoff') is not None else None) for l in cal.get('levels', [])}
            if not cutoffs:
                raise cs.InputError('There is no own calibration yet. Compute one on the page "Calibrate" first.')
            calibration, origin = cal, cs.own_origin(cal)
        elif cut.get('type') == 'manual':
            cutoffs = {}
            for lv, val in (cut.get('cutoffs') or {}).items():
                if str(val).strip() != '':
                    cutoffs[float(lv)] = number(val, 'The cut-off')
            if not cutoffs:
                raise cs.InputError('Type at least one cut-off.')
            origin = 'the user'
        else:
            raise cs.InputError('Choose where the cut-off comes from.')
        res = cs.assign(mol, cutoffs, workflow=label, calibration=calibration, origin=origin)
        res['report'] = cs.report_assign(res)
        res['preview'] = preview(mol)
        return res
    if path == '/api/calibrate':
        res = cs.calibrate(cs.read_table(body.get('text', '')))
        res['report'] = cs.report_calibrate(res)
        return res
    if path == '/api/quit':
        threading.Thread(target=lambda: (time.sleep(0.3), os._exit(0)), daemon=True).start()
        return dict(ok=True)
    raise KeyError(path)


class Handler(BaseHTTPRequestHandler):
    server_version = 'ConformalStereo/' + cs.__version__
    protocol_version = 'HTTP/1.1'

    def log_message(self, fmt, *args):                      # a quiet window
        pass

    def _send(self, code, payload, ctype='application/json; charset=utf-8'):
        data = payload if isinstance(payload, bytes) else json.dumps(cs.clean(payload), allow_nan=False).encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.end_headers()
        self.wfile.write(data)

    def _local(self):
        host = (self.headers.get('Host') or '').split(':')[0].strip('[]').lower()
        if host not in ('127.0.0.1', 'localhost', '::1'):
            self._send(403, dict(error='This page answers on this computer only.'))
            return False
        STATE['last'] = time.time()
        return True

    def do_GET(self):
        if not self._local():
            return
        path = self.path.split('?')[0]
        try:
            if path in ('/', '/index.html'):
                with open(os.path.join(HERE, 'ui.html'), 'rb') as f:
                    self._send(200, f.read(), 'text/html; charset=utf-8')
            elif path == '/api/ping':
                STATE['page_seen'] = True
                self._send(200, dict(ok=True, version=cs.__version__))
            elif path == '/api/workflows':
                self._send(200, dict(version=cs.__version__, python=sys.version.split()[0], workflows=WORKFLOWS))
            elif path == '/api/examples':
                out = {}
                for key, name in (('calibration', 'calibration_example.csv'), ('template', 'calibration_template.csv'), ('typed', 'typed_example.txt')):
                    fn = os.path.join(HERE, 'examples', name)
                    out[key] = io.open(fn, encoding='utf-8').read() if os.path.exists(fn) else ''
                self._send(200, out)
            else:
                self._send(404, dict(error='not found'))
        except Exception as e:                              # noqa
            traceback.print_exc()
            self._send(500, dict(error='The program met an error it did not expect: %s' % e))

    def do_POST(self):
        if not self._local():
            return
        try:
            n = int(self.headers.get('Content-Length') or 0)
            if n > MAX_BODY:
                self._send(413, dict(error='The file is too large (limit 30 MB).'))
                return
            body = json.loads(self.rfile.read(n).decode('utf-8') or '{}')
            self._send(200, api(self.path.split('?')[0], body))
        except cs.InputError as e:
            self._send(400, dict(error=str(e)))
        except KeyError:
            self._send(404, dict(error='not found'))
        except Exception as e:                              # noqa
            traceback.print_exc()
            self._send(500, dict(error='The program met an error it did not expect: %s. The window of the program shows the details.' % e))


WORKFLOWS = []


def serve(port=None, open_browser=True, quiet_exit=180.0):
    global WORKFLOWS
    WORKFLOWS = cs.load_workflows()
    STATE['quiet_exit'] = quiet_exit
    httpd = None
    for p in ([port] if port else list(range(8765, 8790)) + [0]):
        try:
            httpd = ThreadingHTTPServer(('127.0.0.1', p), Handler)
            break
        except OSError:
            continue
    if httpd is None:
        raise SystemExit('ConformalStereo could not open a local address. Is another copy running?')
    httpd.daemon_threads = True
    url = 'http://127.0.0.1:%d/' % httpd.server_address[1]
    print('')
    print('  ConformalStereo %s is running.' % cs.__version__)
    print('  Page:  %s' % url)
    print('  Keep this window open while you work. Close it, or press "Quit" on the page, to stop the program.')
    print('  Everything is computed on this computer; nothing is sent anywhere.')
    print('')
    sys.stdout.flush()

    def watchdog():                                         # stop some minutes after the page was closed
        while True:
            time.sleep(5)
            if STATE['page_seen'] and quiet_exit and time.time() - STATE['last'] > quiet_exit:
                print('  The page was closed. The program stops.')
                sys.stdout.flush()
                os._exit(0)
    threading.Thread(target=watchdog, daemon=True).start()
    if open_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print('  Stopped.')
    return 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='ConformalStereo as a page in the web browser of this computer.')
    ap.add_argument('--port', type=int, help='local port (default: the first free one from 8765)')
    ap.add_argument('--no-browser', action='store_true', help='do not open the browser; print the address only')
    ap.add_argument('--keep', action='store_true', help='keep running after the page was closed')
    a = ap.parse_args()
    sys.exit(serve(a.port, not a.no_browser, 0.0 if a.keep else 180.0))
