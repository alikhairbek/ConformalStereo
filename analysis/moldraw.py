"""Structure drawings inside matplotlib figures.

RDKit lays the molecule out and draws it in its ACS 1996 mode (bond 14.4 pt, line 0.6 pt, labels 10 pt, the
settings of the ACS document style of the drawing programs). The drawing is read back as outlines and placed in the
figure as vector paths, so that structures, text and data marks share one file, one scale and one typeface.
Needs RDKit 2023.09 or newer."""
import re
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Circle
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

_NUM = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?')
_TOK = re.compile(r'([MLQCZz])([^MLQCZz]*)')
ACS_BOND = 14.4


def _path(d):
    verts, codes = [], []
    for cmd, args in _TOK.findall(d):
        v = [float(x) for x in _NUM.findall(args)]
        pts = list(zip(v[0::2], v[1::2]))
        if cmd == 'M':
            for k, pt in enumerate(pts):
                verts.append(pt); codes.append(Path.MOVETO if k == 0 else Path.LINETO)
        elif cmd == 'L':
            for pt in pts:
                verts.append(pt); codes.append(Path.LINETO)
        elif cmd == 'Q':
            for pt in pts:
                verts.append(pt); codes.append(Path.CURVE3)
        elif cmd == 'C':
            for pt in pts:
                verts.append(pt); codes.append(Path.CURVE4)
        else:
            verts.append(verts[-1] if verts else (0.0, 0.0)); codes.append(Path.CLOSEPOLY)
    return np.array(verts, float), codes


def _wedge_from(mol, centre, neighbour):
    """Put the wedge of a stereocentre on the bond to one chosen neighbour (the wedge has to start at the centre)."""
    for b in mol.GetAtomWithIdx(centre).GetBonds():
        if b.GetBeginAtomIdx() == centre:
            b.SetBondDir(Chem.BondDir.NONE)
    b = mol.GetBondBetweenAtoms(centre, neighbour)
    if b.GetBeginAtomIdx() != centre:
        rw = Chem.RWMol(mol)
        order = b.GetBondType()
        rw.RemoveBond(centre, neighbour)
        rw.AddBond(centre, neighbour, order)
        mol = rw.GetMol()
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        b = mol.GetBondBetweenAtoms(centre, neighbour)
    Chem.WedgeBond(b, centre, mol.GetConformer())
    return mol


class Drawing:
    """Outlines of one molecule in points at ACS size (y upwards).

    wedging: 'auto' (RDKit chooses the wedged bonds from the chirality) or 'molblock' (the wedges of the mol block are
    kept as drawn). chiral_hs: show the hydrogen at stereocentres whose three other bonds lie in rings.
    wedge_at: {stereocentre: neighbour} to move a wedge to a chosen bond. open_stereo: stereocentres without a stated
    configuration get a wavy bond."""

    def __init__(self, mol, wedging='auto', chiral_hs=True, wedge_at=None, open_stereo=False, font_file=None):
        mol = Chem.Mol(mol)
        if mol.GetNumConformers() == 0:
            rdDepictor.Compute2DCoords(mol)
        n_atoms = mol.GetNumAtoms()
        md = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True, addChiralHs=chiral_hs, wedgeBonds=(wedging == 'auto'))
        if wedging == 'molblock':
            Chem.ReapplyMolBlockWedging(md)
        for centre, neighbour in (wedge_at or {}).items():
            md = _wedge_from(md, centre, neighbour)
        d = rdMolDraw2D.MolDraw2DSVG(-1, -1)
        o = d.drawOptions()
        rdMolDraw2D.SetACS1996Mode(o, Draw.MeanBondLength(mol))
        o.useBWAtomPalette()
        o.clearBackground = False
        o.prepareMolsBeforeDrawing = False
        if open_stereo:
            o.unspecifiedStereoIsUnknown = True
        if font_file:
            o.fontFile = font_file
        d.DrawMolecule(md)
        d.FinishDrawing()
        svg = d.GetDrawingText()
        w, h = [float(x) for x in re.search(r"width='([\d.]+)px' height='([\d.]+)px'", svg).groups()]
        self.items = []
        for tag in re.findall(r'<path\b[^>]*>', svg):
            dm = re.search(r"\bd='([^']*)'", tag)
            if not dm:
                continue
            verts, codes = _path(dm.group(1))
            verts[:, 1] = h - verts[:, 1]
            style = (re.search(r"style='([^']*)'", tag) or [None, ''])[1]
            st = dict(x.split(':', 1) for x in style.split(';') if ':' in x)
            cls = (re.search(r"class='([^']*)'", tag) or [None, ''])[1]
            if re.search(r"\bfill='([^']*)'", tag):          # glyph of an atom label
                self.items.append(dict(kind='glyph', verts=verts, codes=codes, cls=cls))
            elif st.get('fill', 'none') != 'none':           # bold wedge
                self.items.append(dict(kind='wedge', verts=verts, codes=codes, cls=cls, lw=float(st.get('stroke-width', '0').replace('px', ''))))
            else:
                self.items.append(dict(kind='line', verts=verts, codes=codes, cls=cls, lw=float(st.get('stroke-width', '0.6').replace('px', ''))))
        self.atoms = np.array([[d.GetDrawCoords(i).x, h - d.GetDrawCoords(i).y] for i in range(n_atoms)])
        v = np.vstack([it['verts'] for it in self.items])
        self.x0, self.y0 = v[:, 0].min(), v[:, 1].min()        # crop to the ink
        for it in self.items:
            it['verts'] = it['verts'] - np.array([self.x0, self.y0])
        self.atoms = self.atoms - np.array([self.x0, self.y0])
        self.width, self.height = v[:, 0].max() - self.x0, v[:, 1].max() - self.y0
        self.mol = md

    def size(self, scale=1.0):
        return self.width * scale, self.height * scale

    def draw(self, ax, x, y, scale=1.0, align=(0.5, 0.5), color='black', min_lw=0.5, marks=None, mark_color='#CC79A7', mark_radius=4.4, mark_alpha=0.45,
             transform=None, zorder=3):
        """Draw into an axes whose data unit is the point. (x, y): where the point `align` of the drawing goes.
        marks: atoms that get a filled circle behind the drawing. Returns the box (x0, y0, x1, y1)."""
        tr = transform if transform is not None else ax.transData
        w, h = self.size(scale)
        ox, oy = x - align[0] * w, y - align[1] * h
        off = np.array([ox, oy])
        for i in (marks or []):
            ax.add_patch(Circle(tuple(self.atoms[i] * scale + off), mark_radius * scale, facecolor=mark_color, edgecolor='none', alpha=mark_alpha, zorder=zorder - 1, transform=tr, clip_on=False))
        for it in self.items:
            p = Path(it['verts'] * scale + off, it['codes'])
            if it['kind'] == 'glyph':
                ax.add_patch(PathPatch(p, facecolor=color, edgecolor='none', lw=0, zorder=zorder, transform=tr, clip_on=False))
            elif it['kind'] == 'wedge':
                ax.add_patch(PathPatch(p, facecolor=color, edgecolor=color, lw=max(it['lw'] * scale, 0.1), joinstyle='miter', zorder=zorder, transform=tr, clip_on=False))
            else:
                ax.add_patch(PathPatch(p, facecolor='none', edgecolor=color, lw=max(it['lw'] * scale, min_lw), capstyle='butt', joinstyle='miter', zorder=zorder, transform=tr, clip_on=False))
        return ox, oy, ox + w, oy + h


def point_axes(fig, ax, full_width=False, margin=4.0):
    """Freeze the layout and give an empty axes the point as its data unit (1 unit = 1/72 in). Returns width, height.
    full_width: stretch the axes over the whole width of the figure (less the margin in points) first."""
    fig.draw_without_rendering()
    fig.set_layout_engine('none')
    box = ax.get_position()
    fw, fh = fig.get_size_inches()
    if full_width:
        m = margin / 72.0 / fw
        ax.set_position([m, box.y0, 1.0 - 2 * m, box.height])
        box = ax.get_position()
    w, h = box.width * fw * 72.0, box.height * fh * 72.0
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    return w, h
