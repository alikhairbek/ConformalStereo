"""Supplementary Figures S1 to S4.
Usage: python3 -I si_figures.py <final tables dir> <pooled_candidates.csv> <out dir>"""
import sys, os, csv, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cslib as cs

final, pooled, out = sys.argv[1:4]
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
BCOL = dict(zip(BENCH, [BLUE, GREEN, ORANGE, VERM, PURPLE]))
W = 7.0
FINAL = 'likelihood ratio per shift'


def read(folder, name):
    return list(csv.DictReader(open(os.path.join(folder, name))))


def save(fig, name):
    fig.savefig(os.path.join(out, name + '.png'), metadata={'Software': None})
    fig.savefig(os.path.join(out, name + '.pdf'), metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)
    print('written', name)


rec = read(final, 'final_records.csv')
true = lambda v: v in ('1', 'True', 'true')

# ------------------------------------------------------------------ S1 reliability per benchmark
EDGES = [0.0, 0.5, 0.9, 0.99, 0.999, 0.99999]
TICKS = ['0', '50', '90', '99', '99.9', '99.999']
fig, axes = plt.subplots(1, 5, figsize=(W, 1.95), sharey=True)
print('S1: stated probability of the top candidate, classes with ten records or more: records, mean stated, share right')
for ax, b in zip(axes, BENCH):
    rr = [r for r in rec if r['benchmark'] == b]
    xs, st, ri = [], [], []
    for i, lo in enumerate(EDGES):
        hi = EDGES[i + 1] if i + 1 < len(EDGES) else 2.0
        cl = [r for r in rr if lo <= float(r['stated_probability_top']) < hi]
        if len(cl) >= 10:
            xs.append(i); st.append(100 * np.mean([float(r['stated_probability_top']) for r in cl])); ri.append(100 * np.mean([true(r['top_right']) for r in cl]))
            print('   %-8s %-7s %4d %6.2f %5.1f' % (b, TICKS[i], len(cl), st[-1], ri[-1]))
    ax.plot(xs, st, color=GREY, lw=1.2, ls=DASH, marker='o', ms=4, mfc='white', mec=GREY, mew=1.0, label='stated')
    ax.plot(xs, ri, color=BLUE, lw=1.8, marker='o', ms=5, mec='white', mew=1.0, label='right')
    ax.text(xs[-1], ri[-1] - 9, '%.0f%%' % ri[-1], ha='center', va='top', fontsize=7.5)
    ax.set_title(b)
    ax.set_xticks(range(6)); ax.set_xticklabels(TICKS, rotation=45, ha='right', rotation_mode='anchor')
    ax.set_xlim(-0.4, 5.4); ax.set_ylim(0, 105); ax.set_yticks([0, 25, 50, 75, 100])
    ax.grid(axis='y', color=LIGHT, lw=0.7); ax.set_axisbelow(True)
    ax.text(0.97, 0.04, '{:,} records'.format(len(rr)), transform=ax.transAxes, ha='right', va='bottom', fontsize=7)
axes[0].set_ylabel('top candidate right (%)')
axes[1].legend(frameon=False, loc='lower left', handlelength=2.2, borderaxespad=0.3)
fig.supxlabel('stated probability of the top candidate (%), lower edge of the class', fontsize=8)
save(fig, 'FigureS1_reliability_benchmarks')

# ------------------------------------------------------------------ S2 share of candidates kept against the number of families
wf = read(final, 'final_workflows.csv')
fig, axes = plt.subplots(1, 2, figsize=(W, 2.7), sharey=True)
for ax, t in zip(axes, ('90', '80')):
    for b in BENCH:
        rr = [r for r in wf if r['benchmark'] == b]
        x = [int(r['families']) for r in rr]
        y = [100 * float(r['mean_set_size_' + t]) / float(r['mean_candidates']) for r in rr]
        ax.scatter(x, y, s=22, color=BCOL[b], edgecolor='white', linewidth=0.6, label=b, zorder=3)
    ax.axhline(float(t), color=GREY, lw=1.0, ls=DASH)
    ax.text(79, float(t) + 1.5, 'a score without information', ha='right', va='bottom', fontsize=7, color='black')
    ax.set_title('%s%% sets' % t)
    ax.set_xlabel('families of the workflow')
    ax.set_xlim(10, 82); ax.set_ylim(0, 102)
    ax.grid(axis='y', color=LIGHT, lw=0.7); ax.set_axisbelow(True)
axes[0].set_ylabel('candidates kept in the set (%)')
axes[1].legend(frameon=False, loc='lower right', handletextpad=0.2, borderaxespad=0.2, ncol=2, columnspacing=1.0)
save(fig, 'FigureS2_families_and_set_size')

# ------------------------------------------------------------------ S3 distance of the actual structure, one curve per workflow
fig, axes = plt.subplots(1, 5, figsize=(W, 2.1), sharey=True)
XMAX = 1.5
for ax, b in zip(axes, BENCH):
    by = collections.defaultdict(list)
    for r in rec:
        if r['benchmark'] == b:
            by[(r['method'], r['geometry'], r['level'])].append((float(r['distance_of_actual_per_shift']), r['family']))
    for k, items in by.items():
        cnt = collections.Counter(f for _, f in items)
        items.sort()
        xs, ys, acc = [0.0], [0.0], 0.0
        for s, f in items:
            acc += 1.0 / cnt[f] / len(cnt)
            xs += [min(s, XMAX), min(s, XMAX)]; ys += [ys[-1], 100 * acc]
        xs.append(XMAX); ys.append(ys[-1])
        ax.plot(xs, ys, color=BCOL[b], lw=0.7, alpha=0.55)
    for lev in (80, 90):
        ax.axhline(lev, color=GREY, lw=0.9, ls=DASH)
    ax.set_title('%s (%d)' % (b, len(by)))
    ax.set_xlim(0, XMAX); ax.set_ylim(0, 101); ax.set_xticks([0, 0.5, 1.0, 1.5])
    ax.grid(axis='y', color=LIGHT, lw=0.7); ax.set_axisbelow(True)
axes[0].set_ylabel('families with the actual structure\nwithin the distance (%)')
axes[0].set_yticks([0, 20, 40, 60, 80, 90, 100])
fig.supxlabel('log-likelihood of the actual structure below the best candidate, per shift', fontsize=8)
save(fig, 'FigureS3_distance_of_the_actual_structure')

# ------------------------------------------------------------------ S4 energy averaging per family
recs_p = cs.load(pooled)
recs_e = cs.load(pooled, 'log_ratio_energy_4')
keys = {(r['wf'], r['case']) for r in recs_e}
recs_p = [r for r in recs_p if (r['wf'], r['case']) in keys]
assert [(r['wf'], r['case']) for r in recs_p] == [(r['wf'], r['case']) for r in recs_e]


def family_sizes(recs, target):
    by = collections.defaultdict(list)
    for r in recs:
        by[r['wf']].append(r)
    unit = collections.defaultdict(list)
    for w, rr in by.items():
        for f in sorted({r['family'] for r in rr}):
            thr = cs.threshold([(x['t'][FINAL], x['family']) for x in rr if x['family'] != f], target)
            for r in rr:
                if r['family'] == f:
                    unit[(w, f)].append(cs.apply_set(r, FINAL, thr)[1])
    fam = collections.defaultdict(list)
    for (w, f), v in unit.items():
        fam[f].append(float(np.mean(v)))
    return {f: float(np.mean(v)) for f, v in fam.items()}


fam_bench = {}
for r in recs_e:
    fam_bench.setdefault(r['family'], r['bench'])
fig, axes = plt.subplots(1, 2, figsize=(W * 0.82, 2.9))
print('S4: mean set size per family, plain against energy-averaged likelihood')
for ax, t in zip(axes, (0.90, 0.80)):
    a, e = family_sizes(recs_p, t), family_sizes(recs_e, t)
    fams = sorted(a)
    smaller = sum(1 for f in fams if e[f] < a[f] - 1e-9); larger = sum(1 for f in fams if e[f] > a[f] + 1e-9)
    print('   %d%%: %d families | smaller with averaged energies %d, larger %d, equal %d | mean %.3f -> %.3f' % (100 * t, len(fams), smaller, larger, len(fams) - smaller - larger, np.mean([a[f] for f in fams]), np.mean([e[f] for f in fams])))
    hi = max(max(a.values()), max(e.values())) * 1.05
    ax.plot([1, hi], [1, hi], color=GREY, lw=1.0, ls=DASH)
    for lab, bs, col in (('DP4-AI', ('DP4-AI',), ORANGE), ('Cam2017 and Cam2019', ('Cam2017', 'Cam2019'), VERM)):
        ff = [f for f in fams if fam_bench[f] in bs]
        ax.scatter([a[f] for f in ff], [e[f] for f in ff], s=20, color=col, edgecolor='white', linewidth=0.6, label=lab, zorder=3)
    ax.set_xscale('log'); ax.set_yscale('log')
    tk = [1, 2, 4, 8, 16, 32]
    ax.set_xticks(tk); ax.set_xticklabels([str(x) for x in tk]); ax.set_yticks(tk); ax.set_yticklabels([str(x) for x in tk])
    ax.minorticks_off()
    ax.set_xlim(0.9, hi); ax.set_ylim(0.9, hi)
    ax.set_title('%d%% sets' % round(100 * t))
    ax.set_xlabel('structures in the set, plain likelihood')
    ax.text(0.05, 0.95, 'smaller in %d families,\nlarger in %d' % (smaller, larger), transform=ax.transAxes, ha='left', va='top', fontsize=7.5)
axes[0].set_ylabel('structures in the set,\nenergies averaged')
axes[1].legend(frameon=False, loc='lower right', handletextpad=0.2, borderaxespad=0.2)
save(fig, 'FigureS4_energy_averaging_per_family')
