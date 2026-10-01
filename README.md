# ozangenomics

Computational analyses in cheminformatics and proteomics.

## Analyses

### DNA damage, oxidation and repair panel for plasma LC-MS/MS

A 118-protein target list built against UniProtKB/Swiss-Prot, digested in silico,
checked for peptide uniqueness across all 20,525 human entries, and tiered by
whether each protein is realistically measurable in plasma.

**Read the analysis: [`docs/DNA_damage_plasma_panel.md`](docs/DNA_damage_plasma_panel.md)**

The central finding is negative and it matters for study design: 71 of the 118
proteins are nuclear or cytosolic enzymes that do not circulate at a level any
mass spectrometer reaches from plasma. Every protein in double-strand break
repair, nucleotide excision repair, mismatch repair and base excision repair
falls in that group, OGG1, APEX1 and PARP1 among them. Twenty-eight proteins are
realistically measurable, and they are the consequences of DNA damage rather than
the repair machinery: nucleosomal histones from cell death and NETosis, the redox
enzymes that create the damaging environment, and the acute-phase iron and heme
handlers that drive Fenton chemistry.

The exception worth pursuing is H2AX. Its C-terminal tryptic peptide ATQASQEY is
specific to H2AX across the human proteome and contains Ser140, the residue the
literature calls Ser139, so a phospho-PRM assay on it reads the canonical
double-strand-break marker from blood rather than from cells on a slide.

Thirteen of the 28 measurable proteins are erythrocyte-rich, so a hemolysis index
belongs on every sample as a covariate.

```
data/panel_targets.tsv          118 targets with pathway category
data/plasma_evidence.tsv        tier, expected level, assay and confounder per protein
data/human_swissprot.fasta.gz   the 20,525 human Swiss-Prot entries used
src/fetch_swissprot.py          rebuilds that file from the NCBI distribution
src/panel_build.py              digestion, uniqueness, scoring, transitions, 13 checks
src/validate.py                 re-derives every number by an independent route
results/                        proteins, peptides, PRM transitions, QC
```

```bash
python3 src/panel_build.py
python3 src/validate.py
```

### PROTAC architecture vs. *Primula* root constituents

Computational assessment of whether any constituent of primrose (*Primula*) root can
serve as a module of a PROTAC degrader — warhead, linker, or E3 ligand — and what the
structure–activity comparison actually supports.

**Read the analysis: [`docs/PROTAC_primula_SAR.md`](docs/PROTAC_primula_SAR.md)**

#### Headline findings

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

#### Species caveat

The European Pharmacopoeia drug *Primulae radix* is ***P. veris*** and/or ***P. elatior***,
**not *P. vulgaris***. All quantitative constituent data used here derive from those
species, and published work shows saponin pattern is species-discriminating. The
constituent list is treated as genus-level throughout and flagged in the report.

#### Layout

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

#### Reproduce

```bash
pip install rdkit
python3 src/run_analysis.py
```

#### Verification built into the run

Every drawn structure is checked against its literature-reported molecular formula before
any descriptor is reported. All 13 pass — see `results/formula_validation.csv`. Reference
PROTAC structures come from ChEMBL rather than being drawn by hand.

Stereochemistry is asserted only where it could be verified: the phenolic glycosides and
flavonoids carry it, the triterpenoid sapogenins are written constitution-only. Every
descriptor used is constitution-level and Morgan fingerprints are generated without
chirality, so no conclusion depends on unverified stereochemistry.

#### Sources

Constituent identities and quantitation from PubMed-indexed literature — Trendafilova
*et al.* 2026 ([DOI](https://doi.org/10.3390/plants15152259)), Müller *et al.* 2005
([DOI](https://doi.org/10.1016/j.chroma.2005.10.067)), Lacchini *et al.* 2025
([DOI](https://doi.org/10.1111/pbi.70122)). Reference structures from ChEMBL v34.
Full citations in the report.
