# PROTAC architecture vs. *Primula* root constituents — a cheminformatics SAR assessment

Computational assessment of whether any constituent of primrose (*Primula*) root can
serve as a module of a PROTAC degrader — warhead, linker, or E3 ligand — and what the
structure–activity comparison actually supports.

**Read the analysis: [`docs/PROTAC_primula_SAR.md`](docs/PROTAC_primula_SAR.md)**

## Headline findings

- **No structural overlap.** Maximum ECFP4 Tanimoto of any *Primula* constituent to any
  PROTAC reference is **0.133**, well below the ~0.4 threshold implying shared
  bioactivity. There is no natural CRBN or VHL binder in this extract.
- **Two motifs are transferable.** The triterpenoid sapogenin (protoprimulagenin A) is a
  conformationally locked spacer — **0 rotatable bonds, Fsp3 = 1.00**, two orthogonal
  hydroxyl exit vectors — but at cLogP 6.35 it is too lipophilic to graft into an
  already-bRo5 degrader. The primeverose disaccharide delivers **2.9× the polar surface
  area of a PEG3 linker with fewer than half the rotatable bonds**, and its glycosidic
  bond is an enzymatically cleavable release trigger.
- **Intact saponins are outside every permeability envelope.** At the same mass as
  cyclosporin A (~1200 Da), the modelled pentaglycoside carries **14 hydrogen-bond donors
  against cyclosporin's 5**, with TPSA 385 vs 279. Passive permeability is structurally
  precluded, not merely poor.
- **The credible translational angle is delivery, not degradation** — and the SAR argues
  *Primula* is the *wrong genus* for it. Endosomal-escape-enhancing saponins (QS-21,
  SO1861) require a **C-28 carboxylic acid** to carry their second glycan chain. In the
  *Primula* sapogenin class C-28 is consumed by the 13,28-ether bridge, so these saponins
  are constitutionally barred from forming that chain and are locked as monodesmosides —
  the strongly hemolytic subtype. This predicts a narrower therapeutic window than
  *Saponaria*/*Quillaja* saponins. It is a falsifiable prediction; Section 7 of the report
  gives the experiment that would test it.

## Species caveat

The European Pharmacopoeia drug *Primulae radix* is ***P. veris*** and/or ***P. elatior***,
**not *P. vulgaris***. All quantitative constituent data used here derive from those
species, and published work shows saponin pattern is species-discriminating. The
constituent list is treated as genus-level throughout and flagged in the report.

## Layout

```
data/primula_constituents.csv   13 constituents: class, species evidence, SMILES,
                                structure provenance, literature molecular formula
data/protac_reference.csv       10 references: CRBN/VHL ligands, two clinical PROTACs,
                                cyclosporin A benchmark, three canonical linkers
src/chemspace.py                descriptor panel, Ro5/bRo5 classification, ECFP4 +
                                MACCS Tanimoto, exit-vector (linker handle) detection
src/build_models.py             programmatic glycosylation of the sapogenin, 1-6 sugars,
                                to bracket the saponin class property envelope
src/run_analysis.py             orchestration; writes results/ and the console report
results/                        descriptors, similarity matrix, formula validation
docs/PROTAC_primula_SAR.md      the full analysis
```

## Reproduce

```bash
pip install rdkit
python3 src/run_analysis.py
```

## Verification built into the run

Every drawn structure is checked against its literature-reported molecular formula before
any descriptor is reported. All 13 pass — see `results/formula_validation.csv`. Reference
PROTAC structures come from ChEMBL rather than being drawn by hand.

Stereochemistry is asserted only where it could be verified: the phenolic glycosides and
flavonoids carry it, the triterpenoid sapogenins are written constitution-only. Every
descriptor used is constitution-level and Morgan fingerprints are generated without
chirality, so no conclusion depends on unverified stereochemistry.

## Sources

Constituent identities and quantitation from PubMed-indexed literature — Trendafilova
*et al.* 2026 ([DOI](https://doi.org/10.3390/plants15152259)), Müller *et al.* 2005
([DOI](https://doi.org/10.1016/j.chroma.2005.10.067)), Lacchini *et al.* 2025
([DOI](https://doi.org/10.1111/pbi.70122)). Reference structures from ChEMBL v34.
Full citations in the report.
