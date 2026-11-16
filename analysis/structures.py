#!/usr/bin/env python3
"""ConformalStereo: the chemical structures drawn in Figures 1 and 6 and in the graphic for the table of contents.

    python analysis/structures.py [out.csv]      checks every structure and, with a file name, writes the SMILES table

Compound 87 (worked example, Figure 1) is taken from the drawings in the Supporting Information of the DP4+ paper
(Grimblat, Zanardi and Sarotti, J. Org. Chem. 2015, 80, 12526; pages S218 to S222): the four candidates 87a to 87d
share one drawing and differ in its wedges. The molecules of the Cambridge benchmarks (Figure 6) are the actual
candidates of the deposits, read from the stored geometries (Lopinavir1, Zidovudine1, hydroxycurvularin1 and 2,
Similanamide3) and written here as isomeric SMILES.

Similanamide: the stored geometry is the mirror image of the natural product. NMR shifts cannot tell the two apart,
so the benchmark is not affected; the drawing shows the natural enantiomer,
cyclo(Ant-D-Val-L-Leu-D-Ala-N-Me-L-Leu-L-Pip), the structure of PF1171C to which similanamide was revised
(Masuda et al., J. Nat. Prod. 2015, 78, 2286).

Needs RDKit (pip install rdkit). The tables of the paper do not depend on this file."""
import os
import sys
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDepictor, rdCIPLabeler, rdMolDescriptors
from rdkit.Geometry import Point2D

RDLogger.DisableLog('rdApp.*')

# ---------------------------------------------------------------------------------------------- compound 87
# atoms: label, element, x, y in the pixels of the published drawing of 87a (y points down there)
ATOMS_87 = [('C1', 'C', 124, 236), ('C2', 'C', 178, 266), ('C3', 'C', 178, 322), ('C4', 'C', 124, 352), ('C5', 'C', 72, 322), ('C6', 'C', 72, 266),
            ('C7', 'C', 146, 293), ('O2', 'O', 226, 237), ('C8', 'C', 277, 266), ('O8', 'O', 277, 322), ('C9', 'C', 328, 237), ('O9', 'O', 328, 178),
            ('C10', 'C', 379, 148), ('O10', 'O', 379, 88), ('C11', 'C', 430, 178), ('C12', 'C', 379, 266), ('C13', 'C', 430, 237), ('C14', 'C', 481, 266),
            ('C15', 'C', 481, 322), ('C16', 'C', 430, 352), ('C17', 'C', 379, 322)]
# bonds: begin, end, order (a wedge starts at the stereocentre, which is the first atom)
BONDS_87 = [('C1', 'C2', 1), ('C2', 'C3', 1), ('C3', 'C4', 1), ('C4', 'C5', 1), ('C5', 'C6', 2), ('C6', 'C1', 1), ('C1', 'C7', 1), ('C4', 'C7', 1),
            ('C2', 'O2', 1), ('O2', 'C8', 1), ('C8', 'O8', 2), ('C8', 'C9', 1), ('C9', 'O9', 1), ('O9', 'C10', 1), ('C10', 'O10', 2), ('C10', 'C11', 1),
            ('C9', 'C12', 1), ('C12', 'C13', 2), ('C13', 'C14', 1), ('C14', 'C15', 2), ('C15', 'C16', 1), ('C16', 'C17', 2), ('C17', 'C12', 1)]
# wedges as drawn: 1 = bold (towards the viewer), 6 = hashed (away); 'bridge' = C1-C7 and C4-C7, 'ester' = C2-O
WEDGES_87 = {'87a': {'bridge': 1, 'ester': 6}, '87b': {'bridge': 6, 'ester': 1}, '87c': {'bridge': 6, 'ester': 6}, '87d': {'bridge': 1, 'ester': 1}}
# what the wedges mean (checked below with RDKit): position of the ester on the bicycle and descriptors of C1, C2, C4, C9
EXPECTED_87 = {'87a': ('endo', 'RRRS'), '87b': ('endo', 'SSSS'), '87c': ('exo', 'SRSS'), '87d': ('exo', 'RSRS')}


def molblock_87(name):
    w = WEDGES_87[name]
    idx = {a[0]: i + 1 for i, a in enumerate(ATOMS_87)}
    lines = [name, '  after the Supporting Information of the DP4+ paper', '', '%3d%3d  0  0  0  0  0  0  0  0999 V2000' % (len(ATOMS_87), len(BONDS_87))]
    for lab, el, x, y in ATOMS_87:
        lines.append('%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0' % (x / 56.0 * 1.5, -y / 56.0 * 1.5, 0.0, el))
    for a, b, o in BONDS_87:
        st = w['bridge'] if (a, b) in (('C1', 'C7'), ('C4', 'C7')) else w['ester'] if (a, b) == ('C2', 'O2') else 1 if (a, b) == ('C9', 'O9') else 0
        lines.append('%3d%3d%3d%3d' % (idx[a], idx[b], o, st))
    return '\n'.join(lines + ['M  END']) + '\n'


def compound_87(name):
    """One of 87a, 87b, 87c, 87d with the layout and the wedges of the published drawing."""
    return Chem.MolFromMolBlock(molblock_87(name))


# ---------------------------------------------------------------------------------------------- Cambridge molecules
# isomeric SMILES of the actual candidates, from the geometries of the deposits (canonical form written by RDKit)
FROM_DEPOSIT = {
    'Zidovudine': 'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O',
    'Similanamide': 'CC(C)C[C@@H]1C(=O)N[C@@H](C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)Nc2ccccc2C(=O)N2CCCC[C@@H]2C(=O)N1C',
    'Hydroxycurvularin (11R,15S)': 'C[C@H]1CCC[C@@H](O)CC(=O)c2c(O)cc(O)cc2CC(=O)O1',
    'Hydroxycurvularin (11S,15S)': 'C[C@H]1CCC[C@H](O)CC(=O)c2c(O)cc(O)cc2CC(=O)O1',
    'Lopinavir': 'Cc1cccc(C)c1OCC(=O)N[C@@H](Cc1ccccc1)[C@@H](O)C[C@H](Cc1ccccc1)NC(=O)[C@H](C(C)C)N1CCCNC1=O',
}
# the same drugs written independently of the deposits, as a check of the geometries
REFERENCE = {
    'Zidovudine': 'CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)CO)N=[N+]=[N-]',
    'Lopinavir': 'CC1=C(C(=CC=C1)C)OCC(=O)N[C@@H](CC2=CC=CC=C2)[C@H](C[C@H](CC3=CC=CC=C3)NC(=O)[C@H](C(C)C)N4CCCNC4=O)O',
}
STEP_X, STEP_Y = 1.5 * np.cos(np.pi / 6), 1.5 * np.sin(np.pi / 6)


def canonical(smiles, stereo=True):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=stereo)


def mirror(smiles):
    m = Chem.MolFromSmiles(smiles)
    for a in m.GetAtoms():
        if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            a.InvertChirality()
    return Chem.MolToSmiles(m)


def cip(mol):
    m = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(m)
    return {a.GetIdx(): a.GetProp('_CIPCode') for a in m.GetAtoms() if a.HasProp('_CIPCode')}


def stereocentres(mol):
    """Stereocentres of a structure, with or without a stated configuration (the atoms marked in the figures)."""
    return sorted(i for i, _ in Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))


def _zigzag(mol):
    """Lay the atoms with atom-map numbers 1, 2, 3 ... on a horizontal zigzag; RDKit places the other atoms."""
    chain = sorted((a.GetAtomMapNum(), a.GetIdx()) for a in mol.GetAtoms() if a.GetAtomMapNum())
    fixed = {idx: Point2D(k * STEP_X, STEP_Y if k % 2 else 0.0) for k, (_, idx) in enumerate(chain)}
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    rdDepictor.Compute2DCoords(mol, coordMap=fixed, canonOrient=False)
    return mol


def lopinavir():
    """Main chain drawn as a zigzag from the phenoxy group to the cyclic urea."""
    m = _zigzag(Chem.MolFromSmiles('Cc1cccc(C)[c:1]1[O:2][CH2:3][C:4](=O)[NH:5][C@@H:6](Cc1ccccc1)[C@@H:7](O)[CH2:8][C@H:9](Cc1ccccc1)'
                                   '[NH:10][C:11](=O)[C@H:12](C(C)C)[N:13]1CCCNC1=O'))
    assert Chem.MolToSmiles(m) == canonical(FROM_DEPOSIT['Lopinavir'])
    return m


LOPINAVIR_WEDGES = {24: 25}      # the wedge of this stereocentre goes to its benzyl group, not into the main chain


def zidovudine():
    """Azido group written as N3: its two outer nitrogen atoms are left out of the drawing."""
    m = Chem.MolFromSmiles(FROM_DEPOSIT['Zidovudine'])
    hit = m.GetSubstructMatch(Chem.MolFromSmarts('[#6]-[N;X2]=[N+]=[N-]'))
    rw = Chem.RWMol(m)
    for i in sorted(hit[2:], reverse=True):
        rw.RemoveAtom(i)
    n = rw.GetAtomWithIdx(hit[1] - sum(1 for i in hit[2:] if i < hit[1]))
    n.SetNoImplicit(True)
    n.SetNumExplicitHs(0)
    n.SetProp('_displayLabel', 'N<sub>3</sub>')
    n.SetProp('_displayLabelW', 'N<sub>3</sub>')
    out = rw.GetMol()
    out.UpdatePropertyCache(strict=False)
    rdDepictor.Compute2DCoords(out)
    return out


def similanamide_residues(mol):
    """Descriptor of the alpha carbon of each residue of the cyclic hexapeptide."""
    ring = mol.GetRingInfo()
    out = {}
    for i, label in cip(mol).items():
        a = mol.GetAtomWithIdx(i)
        n = [x for x in a.GetNeighbors() if x.GetSymbol() == 'N'][0]
        side = [x for x in a.GetNeighbors() if x.GetSymbol() == 'C' and not any(y.GetSymbol() == 'O' for y in x.GetNeighbors())][0]
        if any(len(r) == 6 and i in r and n.GetIdx() in r for r in ring.AtomRings()):
            out['Pip'] = label
        elif side.GetDegree() == 1:
            out['Ala'] = label
        elif side.GetDegree() == 3:
            out['Val'] = label
        elif n.GetDegree() == 3:
            out['N-Me-Leu'] = label
        else:
            out['Leu'] = label
    return out


def similanamide():
    """The natural enantiomer: D-Val, L-Leu, D-Ala, N-Me-L-Leu, L-Pip (the deposit holds its mirror image)."""
    m = Chem.MolFromSmiles(mirror(FROM_DEPOSIT['Similanamide']))
    assert similanamide_residues(m) == {'Val': 'R', 'Leu': 'S', 'Ala': 'R', 'N-Me-Leu': 'S', 'Pip': 'S'}
    rdDepictor.Compute2DCoords(m)
    return m


def hydroxycurvularin():
    """One drawing for the two epimers of the benchmark: the configuration of the carbinol carbon is left open."""
    m = Chem.MolFromSmiles('C[C@H]1CCCC(O)CC(=O)c2c(O)cc(O)cc2CC(=O)O1')
    rdDepictor.Compute2DCoords(m)
    # the ketone oxygen and the neighbouring phenol oxygen come out side by side: turn the C=O bond away by 30 degrees
    ketone_c, ketone_o, phenol_o = m.GetSubstructMatch(Chem.MolFromSmarts('[#6;R](=[O])(-[#6;R;!a])-c:c-[OX2H]'))[0:1] + m.GetSubstructMatch(Chem.MolFromSmarts('[O]=[#6;R](-[#6;R;!a])-c:c-[OX2H]'))[0:1] + m.GetSubstructMatch(Chem.MolFromSmarts('[OX2H]-c:c-[#6;R]=O'))[0:1]
    conf = m.GetConformer()
    c, o, far = (np.array(conf.GetAtomPosition(i))[:2] for i in (ketone_c, ketone_o, phenol_o))
    best = None
    for angle in (np.radians(30), -np.radians(30)):
        rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        new = c + rot @ (o - c)
        if best is None or np.linalg.norm(new - far) > np.linalg.norm(best - far):
            best = new
    conf.SetAtomPosition(ketone_o, (float(best[0]), float(best[1]), 0.0))
    return m


FIGURE_6 = {'Zidovudine': zidovudine, 'Similanamide': similanamide, 'Hydroxycurvularin': hydroxycurvularin, 'Lopinavir': lopinavir}


# ---------------------------------------------------------------------------------------------- checks
def face_87(mol):
    """endo or exo from a three-dimensional model: an exo substituent lies on the side of the one-carbon bridge."""
    mh = Chem.AddHs(mol)
    assert AllChem.EmbedMolecule(mh, randomSeed=7) == 0
    AllChem.MMFFOptimizeMolecule(mh, maxIters=2000)
    pos = mh.GetConformer().GetPositions()
    i = {a[0]: k for k, a in enumerate(ATOMS_87)}
    to_bridge = np.linalg.norm(pos[i['O2']] - pos[i['C7']])
    to_alkene = np.linalg.norm(pos[i['O2']] - (pos[i['C5']] + pos[i['C6']]) / 2)
    return 'exo' if to_bridge < to_alkene else 'endo'


def drawn_stereo(drawing_mol):
    """Read a finished drawing back: forget the chirality of every atom, work it out again from the positions and
    the wedges alone, and return the isomeric SMILES that the drawing shows."""
    m = Chem.Mol(drawing_mol)
    for a in m.GetAtoms():
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    Chem.AssignChiralTypesFromBondDirs(m)
    Chem.SanitizeMol(m)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return Chem.MolToSmiles(Chem.RemoveHs(m))


def _plain(mol):
    m = Chem.Mol(mol)
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(Chem.RemoveHs(m))


def checks(out=sys.stdout):
    ok = True

    def line(passed, text):
        nonlocal ok
        ok = ok and passed
        out.write('%s  %s\n' % ('ok    ' if passed else 'FAILED', text))

    smiles = {}
    for name, (face, labels) in EXPECTED_87.items():
        m = compound_87(name)
        lab = cip(m)
        got = ''.join(lab[k] for k in (0, 1, 3, 10))
        smiles[name] = Chem.MolToSmiles(m)
        line(rdMolDescriptors.CalcMolFormula(m) == 'C17H18O4' and face_87(m) == face and got == labels,
             '%s  C17H18O4, %s ester, C1 C2 C4 C9 = %s   %s' % (name, face_87(m), ' '.join(got), smiles[name]))
    line(len(set(smiles.values())) == 4 and len({canonical(s, False) for s in smiles.values()}) == 1, 'compound 87: four stereoisomers of one constitution')
    line(not any(mirror(s) in smiles.values() for s in smiles.values()), 'compound 87: no two candidates are mirror images (four diastereomers)')
    for name, ref in REFERENCE.items():
        line(canonical(FROM_DEPOSIT[name]) == canonical(ref), '%s: geometry of the deposit = known structure of the drug' % name)
    dep = Chem.MolFromSmiles(FROM_DEPOSIT['Similanamide'])
    line(similanamide_residues(dep) == {'Val': 'S', 'Leu': 'R', 'Ala': 'S', 'N-Me-Leu': 'R', 'Pip': 'R'},
         'Similanamide: the deposit holds L-Val, D-Leu, L-Ala, N-Me-D-Leu, D-Pip, the mirror image of PF1171C')
    line(similanamide_residues(similanamide()) == {'Val': 'R', 'Leu': 'S', 'Ala': 'R', 'N-Me-Leu': 'S', 'Pip': 'S'}, 'Similanamide: drawn as D-Val, L-Leu, D-Ala, N-Me-L-Leu, L-Pip')
    a, b = FROM_DEPOSIT['Hydroxycurvularin (11R,15S)'], FROM_DEPOSIT['Hydroxycurvularin (11S,15S)']
    open11 = Chem.MolToSmiles(hydroxycurvularin())
    line(canonical(a, False) == canonical(b, False) == canonical(open11, False) and a != b and list(cip(hydroxycurvularin()).values()) == ['S']
         and sorted(cip(Chem.MolFromSmiles(a)).values()) == ['R', 'S'] and sorted(cip(Chem.MolFromSmiles(b)).values()) == ['S', 'S'],
         'Hydroxycurvularin: two epimers at the carbinol carbon, (11R,15S) and (11S,15S); the drawing leaves this carbon open')
    for name, n in (('Zidovudine', 3), ('Similanamide', 5), ('Hydroxycurvularin', 2), ('Lopinavir', 4)):
        line(len(stereocentres(FIGURE_6[name]())) == n, '%s: %d stereocentres, %d diastereomers' % (name, n, 2 ** (n - 1)))
    try:
        import moldraw
        for name in EXPECTED_87:
            d = moldraw.Drawing(compound_87(name), wedging='molblock', chiral_hs=False)
            line(drawn_stereo(d.mol) == smiles[name], '%s: the drawing shows this stereoisomer' % name)
        for name, f in FIGURE_6.items():
            m = f()
            d = moldraw.Drawing(m, wedge_at=LOPINAVIR_WEDGES if name == 'Lopinavir' else None, open_stereo=(name == 'Hydroxycurvularin'))
            line(drawn_stereo(d.mol) == _plain(m), '%s: the drawing shows this stereoisomer' % name)
    except ImportError as e:
        out.write('skipped the checks of the drawings (%s)\n' % e)
    return ok


def table():
    rows = [('structure', 'figure', 'smiles', 'source')]
    for name in WEDGES_87:
        rows.append((name, 'Figure 1', Chem.MolToSmiles(compound_87(name)), 'drawing in the Supporting Information of the DP4+ paper'))
    rows.append(('Zidovudine', 'Figure 6', canonical(FROM_DEPOSIT['Zidovudine']), 'geometry Zidovudine1 of the deposit'))
    rows.append(('Similanamide', 'Figure 6', Chem.MolToSmiles(similanamide()), 'mirror image of geometry Similanamide3 of the deposit (natural enantiomer)'))
    rows.append(('Hydroxycurvularin (11R,15S)', 'Figure 6', FROM_DEPOSIT['Hydroxycurvularin (11R,15S)'], 'geometry hydroxycurvularin1 of the deposit'))
    rows.append(('Hydroxycurvularin (11S,15S)', 'Figure 6', FROM_DEPOSIT['Hydroxycurvularin (11S,15S)'], 'geometry hydroxycurvularin2 of the deposit'))
    rows.append(('Lopinavir', 'Figure 6', canonical(FROM_DEPOSIT['Lopinavir']), 'geometry Lopinavir1 of the deposit'))
    return rows


if __name__ == '__main__':
    good = checks()
    if len(sys.argv) > 1:
        import csv
        with open(sys.argv[1], 'w', newline='', encoding='utf-8') as f:
            csv.writer(f, lineterminator='\n').writerows(table())
        print('written', os.path.basename(sys.argv[1]))
    sys.exit(0 if good else 1)
