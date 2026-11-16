#!/usr/bin/env python3
"""ConformalStereo: graphic for the table of contents (3.25 in x 1.75 in, the size asked for by ACS journals).
Usage: toc_graphic.py <folder with the exhibit tables> <out folder>

The worked example of the paper in one picture: the probability that DP4+ states for the two endo esters of
compound 87, and the calibrated 90% set, which keeps both. Written as PNG and TIFF at 600 dpi and as PDF.
Needs Matplotlib, RDKit and, for the TIFF file, Pillow."""
import sys, os, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import structures, moldraw

exhibit, out = sys.argv[1:3]
os.makedirs(out, exist_ok=True)
_HAVE = {f.name for f in font_manager.fontManager.ttflist}
SANS = next(f for f in ('Liberation Sans', 'Arial', 'Helvetica', 'DejaVu Sans') if f in _HAVE)      # Arial has the same metrics as Liberation Sans
plt.rcParams.update({'font.family': SANS, 'font.size': 7.5, 'pdf.fonttype': 42, 'ps.fonttype': 42})
FONT_FILE = font_manager.findfont(SANS)
BLUE, GREEN, VERM, PURPLE, SKY = '#0072B2', '#009E73', '#D55E00', '#CC79A7', '#56B4E9'
W_PT, H_PT = 3.25 * 72, 1.75 * 72

ex = list(csv.DictReader(open(os.path.join(exhibit, 'exhibit_example.csv'))))
kept = [r for r in ex if r['in_90_set'] == '1']                     # the members of the 90% set, nearest first
assert len(kept) == 2 and kept[0]['is_actual'] == '0' and kept[1]['is_actual'] == '1'
assert all(structures.EXPECTED_87[r['candidate']][0] == 'endo' for r in kept)

fig = plt.figure(figsize=(3.25, 1.75))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W_PT)
ax.set_ylim(0, H_PT)
ax.axis('off')
scale = 0.78
draw = [moldraw.Drawing(structures.compound_87(r['candidate']), wedging='molblock', chiral_hs=False, font_file=FONT_FILE) for r in kept]
sw, sh = draw[0].size(scale)
xc = [W_PT / 2 - sw / 2 - 13.0, W_PT / 2 + sw / 2 + 13.0]
# row 1: what DP4+ states for the two structures
bar_y, bar_h = H_PT - 19.5, 13.0
for x, r in zip(xc, kept):
    p = 100 * float(r['stated_probability'])
    ax.add_patch(Rectangle((x - sw / 2, bar_y), sw, bar_h, facecolor='#e6ebef', edgecolor='none'))
    if sw * p / 100 >= 0.25:                                    # a bar thinner than a hairline is left out
        ax.add_patch(Rectangle((x - sw / 2, bar_y), sw * p / 100, bar_h, facecolor=SKY, edgecolor='none'))
    ax.text(x, bar_y + bar_h / 2 - 0.3, 'DP4+: %.2f%%' % p, fontweight='bold', ha='center', va='center', fontsize=8)
# row 2: the calibrated set keeps both
f0, f1 = 4.0, H_PT - 26.0
for face, edge, al, z in ((GREEN, 'none', 0.14, 0), ('none', GREEN, 1.0, 1)):
    ax.add_patch(FancyBboxPatch((4.0, f0), W_PT - 8.0, f1 - f0, boxstyle='round,pad=0,rounding_size=7', facecolor=face, edgecolor=edge, alpha=al, lw=0.9, zorder=z))
ax.text(11.0, f1 - 5.5, 'calibrated 90% set', fontweight='bold', ha='left', va='top', fontsize=8)
for x, r, d in zip(xc, kept, draw):
    mol = structures.compound_87(r['candidate'])
    d.draw(ax, x, f0 + 20.0, scale=scale, align=(0.5, 0.0), marks=structures.stereocentres(mol), mark_color=PURPLE)
    if r['is_actual'] == '1':
        t = ax.text(x + 4.5, f0 + 10.5, 'actual structure', fontweight='bold', ha='center', va='center', fontsize=7.5)
        fig.canvas.draw()
        bb = t.get_window_extent().transformed(ax.transData.inverted())
        ax.add_patch(Circle((bb.x0 - 5.5, f0 + 10.7), 2.6, facecolor=VERM, edgecolor='none'))
for ext, kw in (('png', {'dpi': 600, 'metadata': {'Software': None}}), ('pdf', {'metadata': {'CreationDate': None, 'ModDate': None}})):
    fig.savefig(os.path.join(out, 'TOC_graphic.' + ext), **kw)
plt.close(fig)
try:
    from PIL import Image
    im = Image.open(os.path.join(out, 'TOC_graphic.png')).convert('RGB')
    im.save(os.path.join(out, 'TOC_graphic.tif'), dpi=(600, 600), compression='tiff_lzw')
    print('written TOC_graphic (png, pdf, tif): %d x %d pixels at 600 dpi' % im.size)
except ImportError:
    print('written TOC_graphic (png, pdf); the TIFF file needs Pillow')
