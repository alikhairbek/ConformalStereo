#!/usr/bin/env python3
"""Rebuilds every table of the paper from the pooled table of candidates (about two minutes).

    python analysis/run_all.py              tables and logs
    python analysis/run_all.py --figures    also draws the figures into analysis/figures (not part of the repository)

Needs Python 3.9 or newer with NumPy and SciPy:  pip install numpy scipy
The figures need Matplotlib 3.6 or newer, and the chemical structures in Figures 1 and 6 and in the graphic for the
table of contents need RDKit 2022.09 or newer (Pillow for the TIFF file):  pip install matplotlib rdkit pillow
The results are the files in analysis/tables and analysis/logs; a rerun gives the same files."""
import sys, os, subprocess, time, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
POOLED = os.path.join(REPO, 'data', 'pooled_candidates.csv')
COUPLINGS = os.path.join(REPO, 'data', 'coupling_candidates.csv')
TABLES = os.path.join(HERE, 'tables')
LOGS = os.path.join(HERE, 'logs')
APP = os.path.join(REPO, 'ConformalStereo_App', 'app')


def step(title, script, args, log=None):
    t0 = time.time()
    sys.stdout.write('%-58s' % title); sys.stdout.flush()
    cmd = [sys.executable, '-B', os.path.join(HERE, script)] + args
    if log:
        with open(os.path.join(LOGS, log), 'w', encoding='utf-8') as f:
            code = subprocess.call(cmd, stdout=f)
    else:
        code = subprocess.call(cmd, stdout=subprocess.DEVNULL)
    print('%s  %4.0f s' % ('done  ' if code == 0 else 'FAILED', time.time() - t0))
    if code:
        sys.exit(code)


def main():
    for d in (TABLES, LOGS):
        if not os.path.isdir(d):
            os.makedirs(d)
    step('1/7 the method on the five benchmarks', 'final_method.py', [POOLED, TABLES], '01_final_method.txt')
    step('2/7 temperature and size of the molecule', 'temperature_size.py', [POOLED, TABLES], '02_temperature_size.txt')
    step('3/7 checks of the procedure', 'check_final.py', [POOLED, TABLES], '03_checks.txt')
    step('4/7 the number of candidates in the score', 'explore_candidates.py', [POOLED], '04_explore_candidates.txt')
    step('5/7 analyses behind the figures; coupling constants', 'exhibits.py', [POOLED, TABLES, COUPLINGS, TABLES], '05_exhibits.txt')
    step('6/7 the tables of the paper', 'paper_tables.py', [TABLES, LOGS, TABLES, os.path.join(TABLES, 'paper', 'tables.json'), os.path.join(TABLES, 'paper')])
    step('7/7 the cut-offs for the companion program', 'make_cutoffs.py', [POOLED, APP, os.path.join(APP, 'cutoffs.csv'), os.path.join(TABLES, 'final_workflows.csv')])
    if '--figures' in sys.argv:
        fig = os.path.join(HERE, 'figures')
        if not os.path.isdir(fig):
            os.makedirs(fig)
        rdkit = importlib.util.find_spec('rdkit') is not None
        if rdkit:
            step('    structures: checks and table of SMILES', 'structures.py', [os.path.join(REPO, 'data', 'structures.csv')], '06_structures.txt')
        else:
            print('    RDKit is not installed: Figures 1 and 6 are drawn without structures, no graphic for the table of contents')
        step('    figures 1 to 7', 'figures.py', [TABLES, TABLES, fig])
        step('    figures S1 to S4', 'si_figures.py', [TABLES, POOLED, fig])
        if rdkit:
            step('    graphic for the table of contents', 'toc_graphic.py', [TABLES, fig])
    print('Tables: %s' % TABLES)


if __name__ == '__main__':
    main()
