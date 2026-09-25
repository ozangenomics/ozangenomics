"""
Programmatic construction of Primula-type saponin CLASS MODELS.

Why models rather than named congeners
--------------------------------------
Primula root saponins (primulasaponins I-VI, priverosaponin B 22-acetate,
primulic acids) share one 13,28-epoxy-oleanane sapogenin core carrying a
branched oligosaccharide at C-3. Congeners differ in glycan composition and
linkage. Published structures for individual congeners are not consistently
reported in machine-readable form, so rather than assert a specific
stereo-defined structure we build a SERIES: the sapogenin bearing 1-6 sugar
residues. This brackets the property envelope of the whole class.

This is valid for the question being asked because MW, TPSA, HBD, HBA and
Fsp3 are composition-level properties: two glycans of identical composition
but different linkage regiochemistry give identical values.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors

# Sapogenin core: 13,28-epoxy-oleanane-3,16-diol (protoprimulagenin A),
# constitution only - stereocentres deliberately unassigned.
PROTOPRIMULAGENIN_A = "CC1(C)CCC2(CO6)C(O)CC3(C)C6(CCC4C5(C)CCC(O)C(C)(C)C5CCC34C)C2C1"

# Glycosyl donors written so that atom 0 is the anomeric carbon.
HEXOSE = "C1OC(CO)C(O)C(O)C1O"   # extension oxygen at index 6 (O-4)
PENTOSE = "C1OCC(O)C(O)C1O"      # extension oxygen at index 4 (O-4)
DEOXYHEXOSE = "C1OC(C)C(O)C(O)C1O"  # rhamnose-type; extension oxygen index 5

DONORS = {
    "hexose": (HEXOSE, 6),
    "pentose": (PENTOSE, 4),
    "deoxyhexose": (DEOXYHEXOSE, 5),
}


def _find_c3_hydroxyl(mol: Chem.Mol) -> int:
    """Index of the C-3 hydroxyl oxygen (the one flanking the gem-dimethyl C-4)."""
    patt = Chem.MolFromSmarts("[OX2;H1][CX4;H1][CX4]([CH3])([CH3])")
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        raise ValueError("C-3 hydroxyl not found on sapogenin core")
    return matches[0][0]


def attach(acceptor: Chem.Mol, acceptor_o_idx: int, donor_key: str):
    """Form a glycosidic bond: acceptor O-H + donor anomeric C -> O-glycoside."""
    donor_smiles, donor_ext_o = DONORS[donor_key]
    donor = Chem.MolFromSmiles(donor_smiles)
    offset = acceptor.GetNumAtoms()
    combo = Chem.RWMol(Chem.CombineMols(acceptor, donor))
    combo.AddBond(acceptor_o_idx, offset + 0, Chem.BondType.SINGLE)
    o = combo.GetAtomWithIdx(acceptor_o_idx)
    o.SetNoImplicit(False)
    o.SetNumExplicitHs(0)
    mol = combo.GetMol()
    Chem.SanitizeMol(mol)
    return mol, offset + donor_ext_o


def build_series(max_sugars: int = 6, pattern=("hexose", "pentose", "deoxyhexose")):
    """Sapogenin + n sugars, n = 0..max_sugars, as a linear model glycan."""
    core = Chem.MolFromSmiles(PROTOPRIMULAGENIN_A)
    out = [("protoprimulagenin A (aglycone, 0 sugars)", core)]

    mol = core
    next_o = _find_c3_hydroxyl(core)
    for n in range(1, max_sugars + 1):
        donor_key = pattern[(n - 1) % len(pattern)]
        mol, next_o = attach(mol, next_o, donor_key)
        label = f"Primula-type saponin model, {n} sugar residue{'s' if n > 1 else ''}"
        out.append((label, Chem.Mol(mol)))
    return out


if __name__ == "__main__":
    for name, mol in build_series():
        print(f"{name:58s} {rdMolDescriptors.CalcMolFormula(mol):>18s} "
              f"MW={Descriptors.MolWt(mol):8.2f}")
