"""
Physicochemical + similarity engine for comparing Primula root constituents
against PROTAC chemical space.

Design notes
------------
All descriptors used here (MW, TPSA, HBD, HBA, rotatable bonds, Fsp3, Crippen
logP) are constitution-level properties. For isomeric glycans they depend on
glycan *composition*, not on linkage regiochemistry or anomeric configuration.
That is why constitution-only models of the Primula saponins are adequate for
the property-envelope question being asked, and why we do not assert
stereochemistry we cannot verify against an authentic standard.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import MACCSkeys

RDLogger.DisableLog("rdApp.*")

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Derivatisable handles: where a linker could actually be attached.
EXIT_VECTOR_SMARTS = {
    "aliphatic_OH": "[CX4][OX2H]",
    "phenol_OH": "[c][OX2H]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "primary_aromatic_amine": "[c][NX3;H2]",
    "aliphatic_amine": "[NX3;H1,H2;!$(NC=O)]",
    "ester": "[CX3](=O)[OX2][#6]",
    "anomeric_glycosidic": "[OX2;$(O([CX4;R][OX2;R]))]",
}


@dataclass
class Descriptors_:
    name: str
    formula: str
    mw: float
    heavy_atoms: int
    clogp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    rings: int
    aromatic_rings: int
    fsp3: float
    stereocentres: int
    qed: float
    ro5_violations: int
    space: str


def classify_space(mw: float, clogp: float, hbd: int, hba: int, tpsa: float) -> str:
    """Lipinski / beyond-rule-of-5 bucketing.

    Ro5:      MW<=500, logP<=5, HBD<=5, HBA<=10
    bRo5:     the 'extended' druggable envelope oral macrocycles and PROTACs
              occupy: MW<=1000 and TPSA<=250
    beyond:   everything larger / more polar than that.
    """
    ro5_ok = mw <= 500 and clogp <= 5 and hbd <= 5 and hba <= 10
    if ro5_ok:
        return "Ro5"
    if mw <= 1000 and tpsa <= 250 and hbd <= 10:
        return "bRo5 (PROTAC-like envelope)"
    return "beyond bRo5"


def ro5_violations(mw: float, clogp: float, hbd: int, hba: int) -> int:
    return sum([mw > 500, clogp > 5, hbd > 5, hba > 10])


def describe(mol: Chem.Mol, name: str) -> Descriptors_:
    mw = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    try:
        qed = QED.qed(mol)
    except Exception:
        qed = float("nan")
    return Descriptors_(
        name=name,
        formula=rdMolDescriptors.CalcMolFormula(mol),
        mw=round(mw, 2),
        heavy_atoms=mol.GetNumHeavyAtoms(),
        clogp=round(clogp, 2),
        tpsa=round(tpsa, 1),
        hbd=hbd,
        hba=hba,
        rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol),
        rings=rdMolDescriptors.CalcNumRings(mol),
        aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(mol),
        fsp3=round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        stereocentres=len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)),
        qed=round(qed, 3) if not math.isnan(qed) else float("nan"),
        ro5_violations=ro5_violations(mw, clogp, hbd, hba),
        space=classify_space(mw, clogp, hbd, hba, tpsa),
    )


def morgan_fp(mol: Chem.Mol):
    return _MORGAN.GetFingerprint(mol)


def tanimoto(a: Chem.Mol, b: Chem.Mol) -> float:
    return round(DataStructs.TanimotoSimilarity(morgan_fp(a), morgan_fp(b)), 3)


def maccs_tanimoto(a: Chem.Mol, b: Chem.Mol) -> float:
    return round(
        DataStructs.TanimotoSimilarity(MACCSkeys.GenMACCSKeys(a), MACCSkeys.GenMACCSKeys(b)), 3
    )


def exit_vectors(mol: Chem.Mol) -> dict:
    out = {}
    for label, smarts in EXIT_VECTOR_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        out[label] = len(mol.GetSubstructMatches(patt)) if patt is not None else 0
    out["total_handles"] = (
        out["aliphatic_OH"]
        + out["phenol_OH"]
        + out["carboxylic_acid"]
        + out["primary_aromatic_amine"]
        + out["aliphatic_amine"]
    )
    return out


def load(smiles: str, name: str = "") -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Unparsable SMILES for {name!r}: {smiles}")
    return mol


def row(mol: Chem.Mol, name: str) -> dict:
    d = asdict(describe(mol, name))
    d.update(exit_vectors(mol))
    return d
