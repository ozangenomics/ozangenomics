"""Run the full Primula-vs-PROTAC chemical space comparison."""
from __future__ import annotations

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import chemspace as cs
import build_models
from rdkit import Chem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main():
    primula_rows = read_csv(os.path.join(DATA, "primula_constituents.csv"))
    protac_rows = read_csv(os.path.join(DATA, "protac_reference.csv"))

    primula = []
    formula_checks = []
    for r in primula_rows:
        mol = cs.load(r["smiles"], r["name"])
        d = cs.row(mol, r["name"])
        d["chemical_class"] = r["chemical_class"]
        d["plant_part"] = r["plant_part"]
        d["structure_provenance"] = r["structure_provenance"]
        primula.append((r["name"], mol, d))
        calc = d["formula"]
        lit = r["lit_formula"]
        formula_checks.append((r["name"], lit, calc, "MATCH" if calc == lit else "MISMATCH"))

    # Saponin class models
    for label, mol in build_models.build_series():
        if label.startswith("protoprimulagenin"):
            continue
        d = cs.row(mol, label)
        d["chemical_class"] = "triterpenoid saponin (class model)"
        d["plant_part"] = "root"
        d["structure_provenance"] = "constructed-class-model"
        primula.append((label, mol, d))

    protac = []
    for r in protac_rows:
        mol = cs.load(r["smiles"], r["name"])
        d = cs.row(mol, r["name"])
        d["chemical_class"] = r["protac_module"]
        d["plant_part"] = "-"
        d["structure_provenance"] = r["source"]
        protac.append((r["name"], mol, d))

    # ---- descriptors ----
    fields = list(primula[0][2].keys())
    with open(os.path.join(RESULTS, "descriptors.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["dataset"] + fields)
        w.writeheader()
        for _, _, d in primula:
            w.writerow({"dataset": "primula", **d})
        for _, _, d in protac:
            w.writerow({"dataset": "protac_reference", **d})

    # ---- formula validation ----
    with open(os.path.join(RESULTS, "formula_validation.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["compound", "literature_formula", "computed_formula", "status"])
        w.writerows(formula_checks)

    # ---- similarity ----
    with open(os.path.join(RESULTS, "similarity_matrix.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["primula_constituent", "protac_reference", "protac_module",
                    "tanimoto_ecfp4", "tanimoto_maccs"])
        for pname, pmol, _ in primula:
            for rname, rmol, rd in protac:
                w.writerow([pname, rname, rd["chemical_class"],
                            cs.tanimoto(pmol, rmol), cs.maccs_tanimoto(pmol, rmol)])

    # ---- console summary ----
    print("=" * 108)
    print("FORMULA VALIDATION (drawn structures vs literature-reported molecular formulas)")
    print("=" * 108)
    for name, lit, calc, status in formula_checks:
        flag = " " if status == "MATCH" else "<<<"
        print(f"  {status:9s} {flag:3s} {name:46s} lit={lit:12s} calc={calc}")

    print()
    print("=" * 108)
    print("PROPERTY SPACE")
    print("=" * 108)
    hdr = (f"{'compound':52s} {'MW':>8s} {'cLogP':>7s} {'TPSA':>7s} {'HBD':>4s} "
           f"{'HBA':>4s} {'RTB':>4s} {'Fsp3':>6s} {'QED':>6s}  space")
    print(hdr)
    print("-" * 108)
    for tag, group in (("PROTAC ref", protac), ("PRIMULA", primula)):
        print(f"--- {tag} " + "-" * (104 - len(tag)))
        for _, _, d in group:
            print(f"{d['name'][:52]:52s} {d['mw']:8.1f} {d['clogp']:7.2f} {d['tpsa']:7.1f} "
                  f"{d['hbd']:4d} {d['hba']:4d} {d['rotatable_bonds']:4d} {d['fsp3']:6.2f} "
                  f"{d['qed']:6.2f}  {d['space']}")

    print()
    print("=" * 108)
    print("MAXIMUM ECFP4 TANIMOTO OF EACH PRIMULA CONSTITUENT TO ANY PROTAC REFERENCE")
    print("=" * 108)
    for pname, pmol, _ in primula:
        best = max(((cs.tanimoto(pmol, rmol), rname) for rname, rmol, _ in protac))
        print(f"  {pname[:58]:58s} max Tc = {best[0]:.3f}  (vs {best[1]})")

    print()
    print("=" * 108)
    print("EXIT-VECTOR / LINKER-ATTACHMENT HANDLE COUNT")
    print("=" * 108)
    print(f"{'compound':52s} {'alkOH':>6s} {'ArOH':>5s} {'COOH':>5s} {'ArNH2':>6s} {'amine':>6s} {'TOTAL':>6s}")
    for _, _, d in primula + protac:
        print(f"{d['name'][:52]:52s} {d['aliphatic_OH']:6d} {d['phenol_OH']:5d} "
              f"{d['carboxylic_acid']:5d} {d['primary_aromatic_amine']:6d} "
              f"{d['aliphatic_amine']:6d} {d['total_handles']:6d}")

    print("\nWrote: results/descriptors.csv, results/similarity_matrix.csv, "
          "results/formula_validation.csv")


if __name__ == "__main__":
    main()
