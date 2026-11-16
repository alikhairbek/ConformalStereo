"""ConformalStereo: draft figures 1 to 7 of the paper.
Usage: figures.py <folder with the final tables> <folder with the exhibit tables> <out folder>

Conventions: Wong palette, PNG and PDF at 300 dpi, no black bars, markers or curves, no error bars on bars and no
whiskers on points (intervals are in the tables), text in black. Width 7.0 in (two columns).
"""
import sys, os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

final, exhibit, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
SANS = 'Liberation Sans' if any(f.name == 'Liberation Sans' for f in font_manager.fontManager.ttflist) else 'DejaVu Sans'
plt.rcParams.update({'font.family': SANS, 'font.size': 8, 'axes.titlesize': 8.5, 'axes.labelsize': 8, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'legend.fontsize': 7.5,
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#8a8a8a', 'axes.linewidth': 0.7,
                     'xtick.color': '#8a8a8a', 'ytick.color': '#8a8a8a', 'xtick.labelcolor': 'black', 'ytick.labelcolor': 'black', 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                     'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'savefig.dpi': 300, 'figure.dpi': 110, 'figure.constrained_layout.use': True})
BLUE, ORANGE, GREEN, VERM, PURPLE, SKY = '#0072B2', '#E69F00', '#009E73', '#D55E00', '#CC79A7', '#56B4E9'
GREY, LIGHT = '#8a8a8a', '#e3e3e3'
DASH = (0, (4, 2))
BENCH = ['JOC2015', 'STX150', 'DP4-AI', 'Cam2017', 'Cam2019']
W = 7.0


def read(folder, name):
    return list(csv.DictReader(open(os.path.join(folder, name))))


def save(fig, name):
    fig.savefig(os.path.join(out, name + '.png'), metadata={'Software': None})
    fig.savefig(os.path.join(out, name + '.pdf'), metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)
    print('written', name)


def panel(ax, letter, title, **kw):
    ax.set_title('%s  %s' % (letter, title), **kw)


# ------------------------------------------------------------------ Figure 1: one molecule
ex = read(exhibit, 'exhibit_example.csv')
names = [r['candidate'] for r in ex]
p = [100 * float(r['stated_probability']) for r in ex]
s = [float(r['distance_per_shift']) for r in ex]
c90, c80, nsh = float(ex[0]['cutoff_90']), float(ex[0]['cutoff_80']), int(ex[0]['shifts'])
actual = [int(r['is_actual']) for r in ex].index(1)
fig, axs = plt.subplots(1, 3, figsize=(W, 2.45), gridspec_kw={'width_ratios': [1.0, 1.5, 1.05]})
y = np.arange(len(names))[::-1]
ax = axs[0]
ax.barh(y, p, height=0.55, color=SKY)
for yi, v in zip(y, p):
    ax.text(v - 3 if v >= 60 else v + 3, yi, ('%.2f%%' % v) if v >= 0.005 else 'below 0.01%', va='center', ha='right' if v >= 60 else 'left')
ax.set_yticks(y)
ax.set_yticklabels(['%s%s' % (n, '\n(actual)' if i == actual else '') for i, n in enumerate(names)])
ax.set_xlim(0, 100)
ax.set_ylim(-0.6, len(names) - 0.2)
ax.set_xlabel('probability stated by DP4+ (%)')
panel(ax, 'a', 'What DP4+ states')
ax = axs[1]
ax.axvspan(0, c90, color=GREEN, alpha=0.14, lw=0)
ax.axvspan(0, c80, color=GREEN, alpha=0.22, lw=0)
for c, lab in ((c80, '80%'), (c90, '90%')):
    ax.axvline(c, color=GREEN, lw=1.0, ls=DASH)
    ax.text(c, len(names) - 0.42, lab, ha='center', va='bottom', fontsize=7.3)
ax.text(c90 + 0.22, len(names) - 0.42, 'cut-offs: %.2f and %.2f' % (c80, c90), ha='left', va='bottom', fontsize=7.3)
ax.scatter(s, y, s=46, color=[VERM if i == actual else BLUE for i in range(len(names))], zorder=3, edgecolor='white', linewidth=1.0)
for i, (yi, v) in enumerate(zip(y, s)):
    ax.text(v + 0.09, yi, '%.2f%s' % (v, '  (actual)' if i == actual else ''), va='center', ha='left')
ax.set_yticks(y)
ax.set_yticklabels([])
ax.set_ylim(-0.6, len(names) + 0.1)
axs[0].set_ylim(-0.6, len(names) + 0.1)
ax.set_xlim(-0.08, 2.45)
ax.set_xlabel('log-likelihood below the best candidate, per shift')
panel(ax, 'b', 'Distance of every candidate (%d shifts)' % nsh)
ax = axs[2]
ax.axis('off')
panel(ax, 'c', 'Calibrated sets')
in90 = [r['candidate'] for r in ex if r['in_90_set'] == '1']
in80 = [r['candidate'] for r in ex if r['in_80_set'] == '1']
ax.text(0.0, 0.86, '80% set', fontweight='bold', transform=ax.transAxes)
ax.text(0.0, 0.72, ', '.join(in80), transform=ax.transAxes)
ax.text(0.0, 0.58, 'one structure;\nthe actual one is missed', transform=ax.transAxes, fontsize=7.3, va='center')
ax.text(0.0, 0.34, '90% set', fontweight='bold', transform=ax.transAxes)
ax.text(0.0, 0.20, ', '.join(in90), transform=ax.transAxes)
ax.text(0.0, 0.06, 'two structures;\nthe actual one is kept', transform=ax.transAxes, fontsize=7.3, va='center')
save(fig, 'Figure1_worked_example')

# ------------------------------------------------------------------ Figure 2: reliability
rel = read(final, 'final_reliability.csv')
bins = ['below 50%', '50 to 90%', '90 to 99%', '99 to 99.9%', '99.9 to 99.999%', '99.999% or more']
edges = ['0', '50', '90', '99', '99.9', '99.999']
kinds = ['DP4+, B3LYP geometries', 'MM-DP4+, MMFF geometries', 'xTB-DP4+, xTB geometries', 'DP4, MMFF geometries', 'DP4, HF or DFT geometries', 'DP4, semi-empirical geometries']
fig, axs = plt.subplots(2, 3, figsize=(W, 4.0), sharex=True, sharey=True)
MINN = 10
for ax, kind in zip(axs.ravel(), kinds):
    rows = {r['stated_probability_of_top_candidate']: r for r in rel if r['kind_of_calculation'] == kind}
    xs = [i for i, b in enumerate(bins) if b in rows and int(rows[b]['records']) >= MINN]
    stated = [100 * float(rows[bins[i]]['mean_stated']) for i in xs]
    right = [100 * float(rows[bins[i]]['share_right']) for i in xs]
    ax.plot(xs, stated, color=GREY, lw=1.2, ls=DASH, marker='o', ms=3.5, mfc='white', mec=GREY, label='stated')
    ax.plot(xs, right, color=BLUE, lw=1.8, marker='o', ms=5.5, mec='white', mew=1.0, label='right')
    ax.text(xs[-1] + 0.22, right[-1], '%.0f%%' % right[-1], va='center', ha='left', fontsize=7.5)
    n = sum(int(r['records']) for r in rows.values())
    ax.set_title(kind, fontsize=7.8)
    ax.text(5.9, 4, '%s records' % format(n, ','), ha='right', va='bottom', fontsize=7.3)
    ax.set_ylim(0, 105)
    ax.set_xlim(-0.4, 6.0)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.grid(axis='y', color=LIGHT, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xticks(range(6))
    ax.set_xticklabels(edges, fontsize=7)
for ax in axs[:, 0]:
    ax.set_ylabel('top candidate right (%)')
fig.supxlabel('stated probability of the top candidate (%), lower edge of the class', fontsize=8)
axs[0, 0].legend(loc='lower center', frameon=False, handlelength=2.2, ncol=2, columnspacing=1.2, bbox_to_anchor=(0.45, 0.10))
save(fig, 'Figure2_reliability')

# ------------------------------------------------------------------ Figure 3: temperature and size
tb = read(exhibit, 'exhibit_temperature_bins.csv')
ts = {r['workflows']: r for r in read(final, 'final_temperature_size.csv')}
fig = plt.figure(figsize=(W, 4.0))
gs = fig.add_gridspec(2, 5, height_ratios=[1.0, 1.05])
first = None
for k, b in enumerate(BENCH):
    ax = fig.add_subplot(gs[0, k], sharex=first, sharey=first)
    first = first or ax
    rows = [r for r in tb if r['benchmark'] == b]
    x = np.array([float(r['median_shifts']) for r in rows])
    yv = np.array([float(r['temperature']) for r in rows])
    gx, gy = np.exp(np.log(x).mean()), np.exp(np.log(yv).mean())
    xx = np.array([12.0, 80.0])
    ax.plot(xx, gy * xx / gx, color=GREY, lw=1.0, ls=DASH, label='slope 1')
    ax.plot(x, yv, color=BLUE, lw=1.8, marker='o', ms=5.5, mec='white', mew=1.0, label='three groups of molecules')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(12, 80)
    ax.set_ylim(1.5, 30)
    ax.set_xticks([15, 30, 60])
    ax.set_xticklabels(['15', '30', '60'])
    ax.set_yticks([2, 5, 10, 20])
    ax.set_yticklabels(['2', '5', '10', '20'])
    ax.minorticks_off()
    ax.grid(color=LIGHT, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(('a  ' if k == 0 else '') + b, fontsize=8.2)
    if k:
        plt.setp(ax.get_yticklabels(), visible=False)
    else:
        ax.set_ylabel('temperature')
    if k == 2:
        ax.set_xlabel('number of shifts (median of the group)')
    if k == 0:
        ax.text(46, gy * 46 / gx * 0.80, 'slope 1', fontsize=6.8, ha='left', va='top')
ax = fig.add_subplot(gs[1, :])
rowsb = BENCH + ['all workflows']
yy = np.arange(len(rowsb))[::-1]
for g in (0, 0.5, 1):
    ax.axvline(g, color=LIGHT, lw=0.9, zorder=0)
for yi, b in zip(yy, rowsb):
    r = ts[b]
    est, lo, hi = float(r['exponent']), float(r['lo']), float(r['hi'])
    ax.scatter([est], [yi], s=46, color=BLUE, edgecolor='white', linewidth=1.0, zorder=3)
    ax.text(3.6, yi, '%.2f (%.2f to %.2f)' % (est, lo, hi), va='center', ha='left', fontsize=7.5)
ax.text(3.6, len(rowsb) - 0.45, 'estimate (95% interval)', va='bottom', ha='left', fontsize=7.3)
ax.set_yticks(yy)
ax.set_yticklabels(rowsb)
ax.set_xlim(-1.2, 5.6)
ax.set_xticks([-1, 0, 0.5, 1, 2, 3])
ax.set_xticklabels(['-1', '0', '1/2', '1', '2', '3'])
ax.set_ylim(-0.6, len(rowsb) - 0.1)
ax.spines['bottom'].set_bounds(-1.2, 3.2)
ax.set_xlabel('exponent b of the temperature, T proportional to N$^b$', loc='left')
panel(ax, 'b', 'Fitted exponent')
save(fig, 'Figure3_temperature_size')

# ------------------------------------------------------------------ Figure 4: what the sets look like
sk = read(exhibit, 'exhibit_set_kinds.csv')
fig, ax = plt.subplots(figsize=(W, 3.0))
ypos, labels, yv = [], [], 0.0
segs = [('one structure', GREEN), ('several, not all', SKY), ('every candidate', ORANGE)]
for b in BENCH:
    for t in ('90', '80'):
        r = next(x for x in sk if x['benchmark'] == b and x['target'] == t)
        left = 0.0
        for (lab, col), v in zip(segs, [100 * float(r['one_structure']), 100 * float(r['several_not_all']), 100 * float(r['every_candidate'])]):
            ax.barh(yv, v, left=left, height=0.62, color=col, edgecolor='white', linewidth=1.2)
            if v >= 7:
                ax.text(left + v / 2, yv, '%.0f' % v, ha='center', va='center', fontsize=7.3)
            left += v
        ax.text(101.5, yv, '%.2f of %.1f' % (float(r['mean_set_size']), float(r['mean_candidates'])), va='center', ha='left', fontsize=7.5)
        ypos.append(yv)
        labels.append('%s   %s%%' % (b if t == '90' else '', t))
        yv -= 0.74
    yv -= 0.34
ax.text(101.5, ypos[0] + 0.50, 'mean size', va='bottom', ha='left', fontsize=7.3)
ax.set_yticks(ypos)
ax.set_yticklabels(labels)
ax.set_xlim(0, 113)
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.spines['bottom'].set_bounds(0, 100)
ax.set_xlabel('share of the sets (%)')
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in segs]
fig.legend(handles, [l for l, _ in segs], loc='outside upper center', ncol=3, frameon=False, handlelength=1.2, columnspacing=1.6)
save(fig, 'Figure4_sets')

# ------------------------------------------------------------------ Figure 5: a cut-off does not travel
wc = read(exhibit, 'exhibit_workflow_coverage.csv')
origins = ['own workflow', 'same kind of calculation', 'other benchmarks']
olab = ['own workflow', 'same kind of\ncalculation', 'other\nbenchmarks']
family = lambda m: 'DP4' if m == 'DP4' else 'DP4+ and its variants'      # noqa: E731
fig, axs = plt.subplots(1, 2, figsize=(W, 3.2), sharey=True)
rngj = np.random.default_rng(7)
for ax, t, letter in zip(axs, ('0.9', '0.8'), 'ab'):
    for i, o in enumerate(origins):
        rows = [r for r in wc if r['target'] == t and r['cutoff_from'] == o]
        for lab, col in (('DP4', ORANGE), ('DP4+ and its variants', BLUE)):
            v = [100 * float(r['coverage']) for r in rows if family(r['method']) == lab]
            x = i + rngj.uniform(-0.24, 0.24, size=len(v))
            ax.scatter(x, v, s=17, color=col, alpha=0.85, edgecolor='white', linewidth=0.5, zorder=3)
        below = sum(1 for r in rows if 100 * float(r['coverage']) < 100 * float(t) - 5)
        ax.text(i, 33.5, '%d of %d' % (below, len(rows)), ha='center', va='bottom', fontsize=7.5)
    ax.axhline(100 * float(t), color=GREY, lw=1.0, ls=DASH, zorder=1)
    ax.axhspan(100 * float(t) - 5, 100 * float(t), color=LIGHT, alpha=0.6, lw=0, zorder=0)
    ax.set_xticks(range(3))
    ax.set_xticklabels(olab)
    ax.set_xlim(-0.55, 2.5)
    ax.set_ylim(32, 102)
    ax.text(-0.5, 38.0, 'workflows more than 5 points below the target:', ha='left', va='bottom', fontsize=7.3)
    panel(ax, letter, 'Target %s%%, cut-off from' % ('90' if t == '0.9' else '80'))
axs[0].set_ylabel('coverage of a workflow (%)')
handles = [plt.Line2D([], [], marker='o', ls='', color=c, ms=5) for c in (ORANGE, BLUE)] + [plt.Line2D([], [], color=GREY, lw=1.0, ls=DASH)]
fig.legend(handles, ['DP4', 'DP4+ and its variants', 'target'], loc='outside upper center', ncol=3, frameon=False, handletextpad=0.4, columnspacing=1.6)
save(fig, 'Figure5_transfer')

# ------------------------------------------------------------------ Figure 6: stress tests and misses
st = [r for r in read(final, 'final_stress.csv') if r['test'] in ('size', 'candidates')]
fm = read(exhibit, 'exhibit_family_misses.csv')
fig = plt.figure(figsize=(W, 3.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1.0], height_ratios=[1.25, 1.0])
ax = fig.add_subplot(gs[:, 0])
pairs = [('the half with fewer shifts', 'fewer shifts →\nmore shifts'), ('the half with more shifts', 'more shifts →\nfewer shifts'),
         ('at most four candidates', 'at most 4 candidates →\n5 or more'), ('five or more candidates', '5 or more candidates →\nat most 4')]
yv, ypos, labels, done = 0.0, [], [], False
for cal, lab in pairs:
    ymid = []
    for t in ('0.9', '0.8'):
        d = {r['score']: 100 * (float(r['coverage']) - float(t)) for r in st if r['calibrated_on'] == cal and r['target'] == t}
        a, bb = d['likelihood ratio'], d['likelihood ratio per shift']
        ax.plot([a, bb], [yv, yv], color=GREY, lw=1.0, zorder=1)
        ax.scatter([a], [yv], s=40, color=SKY, edgecolor='white', linewidth=1.0, zorder=3, label=None if done else 'plain distance')
        ax.scatter([bb], [yv], s=40, color=BLUE, edgecolor='white', linewidth=1.0, zorder=3, label=None if done else 'distance per shift')
        done = True
        ax.text(-18.6, yv, '%s%%' % ('90' if t == '0.9' else '80'), va='center', ha='left', fontsize=7.3)
        ymid.append(yv)
        yv -= 0.72
    ypos.append(sum(ymid) / 2)
    labels.append(lab)
    yv -= 0.40
ax.axvline(0, color=GREY, lw=1.0, ls=DASH, zorder=0)
ax.set_yticks(ypos)
ax.set_yticklabels(labels, fontsize=7.3)
ax.set_xlim(-19, 19)
ax.set_xticks([-15, -10, -5, 0, 5, 10, 15])
ax.set_xlabel('coverage minus target (percentage points)')
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)
ax.legend(loc='upper right', frameon=False, handletextpad=0.2, borderaxespad=0.0)
panel(ax, 'a', 'Calibrated on → tested on')
for k, (b, letter) in enumerate((('Cam2017', 'b'), ('Cam2019', 'c'))):
    ax = fig.add_subplot(gs[k, 1])
    allr = [r for r in fm if r['benchmark'] == b]
    rows = sorted([r for r in allr if int(r['missed_90']) > 0], key=lambda r: float(r['share_missed_90']))
    v = [100 * float(r['share_missed_90']) for r in rows]
    yb = np.arange(len(rows))
    ax.barh(yb, v, height=0.62, color=BLUE)
    for yi, r, val in zip(yb, rows, v):
        ax.text(val + 2, yi, '%s of %s' % (r['missed_90'], r['records']), va='center', ha='left', fontsize=7)
    ax.set_yticks(yb)
    ax.set_yticklabels([r['family'].split(':')[-1] for r in rows], fontsize=7)
    ax.set_xlim(0, 125)
    ax.set_xticks([0, 50, 100])
    ax.spines['bottom'].set_bounds(0, 100)
    ax.tick_params(axis='y', length=0)
    if k == 1:
        ax.set_xlabel('records missed by the 90% set (%)')
    ax.text(124, -0.3, '%d of %d families\nnever missed' % (len(allr) - len(rows), len(allr)), ha='right', va='bottom', fontsize=7.3)
    panel(ax, letter, 'Misses in %s' % b)
save(fig, 'Figure6_stress_and_misses')

# ------------------------------------------------------------------ Figure 7: set size ranks workflows
wf = read(final, 'final_workflows.csv')
fig, axs = plt.subplots(1, 3, figsize=(W, 2.9))


def cloud(ax, rows, groups, ycol):
    for lab, col, test in groups:
        rr = [r for r in rows if test(r)]
        ax.scatter([100 * float(r['top_candidate_right']) for r in rr], [float(r[ycol]) for r in rr], s=30, color=col, edgecolor='white', linewidth=0.8, label=lab, zorder=3)
    ax.set_xlabel('top candidate right (%)')
    ax.grid(color=LIGHT, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc='upper right', frameon=False, handletextpad=0.1, borderaxespad=0.1)


ax = axs[0]
rows = [r for r in wf if r['benchmark'] == 'JOC2015']
cloud(ax, rows, (('gas phase', ORANGE, lambda r: r['level'].endswith('/gas')), ('PCM', BLUE, lambda r: r['level'].endswith('/PCM'))), 'mean_set_size_90')
best = next(r for r in rows if r['level'] == 'mPW1PW91/6-31+G**/PCM')
ax.text(100 * float(best['top_candidate_right']) + 0.4, float(best['mean_set_size_90']) - 0.055, 'PCM/mPW1PW91/6-31+G**', fontsize=6.8, ha='right', va='top')
ax.set_ylabel('mean size of the set')
ax.set_ylim(0.86, 2.1)
panel(ax, 'a', 'JOC2015, 90% sets')
ax = axs[1]
rows = [r for r in wf if r['benchmark'] == 'Cam2019']
cloud(ax, rows, (('MMFF geometries', ORANGE, lambda r: r['geometry'] == 'MMFF'), ('PM6 or PM7', PURPLE, lambda r: r['geometry'] in ('PM6', 'PM7')),
                 ('HF or DFT', BLUE, lambda r: r['geometry'] in ('HF', 'B3LYP', 'M06-2X'))), 'mean_set_size_80')
ax.set_ylim(0.9, 4.9)
panel(ax, 'b', 'Cam2019, 80% sets')
ax = axs[2]
rows = [r for r in wf if r['benchmark'] == 'STX150']
cloud(ax, rows, (('DP4+', BLUE, lambda r: r['method'] == 'DP4+'), ('xTB-DP4+', GREEN, lambda r: r['method'] == 'xTB-DP4+'), ('MM-DP4+', ORANGE, lambda r: r['method'] == 'MM-DP4+')), 'mean_set_size_95')
ax.set_ylim(0.97, 1.36)
panel(ax, 'c', 'STX150, 95% sets')
save(fig, 'Figure7_workflows')
