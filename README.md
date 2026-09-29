# ozangenomics — computational biology and cheminformatics analyses

Two self-contained studies. Each carries its own data, code, generated results and report.

| Study | Question | Report |
|---|---|---|
| **DNA stress and repair panel** | Can 8-OHdG be extended into a panel that quantifies DNA damage *and* the repair capacity protecting against it? | [markdown](docs/DNA_stress_repair_panel.md) · [PDF](docs/DNA_stress_repair_panel.pdf) |
| **PROTAC vs *Primula* SAR** | Can any constituent of primrose root serve as a PROTAC module? | [`docs/PROTAC_primula_SAR.md`](docs/PROTAC_primula_SAR.md) |

---

# Study 1 — A multiplex panel for DNA stress and DNA repair capacity

Built outward from **8-OHdG** (8-oxo-7,8-dihydro-2'-deoxyguanosine) to a three-module
panel that measures oxidative and alkylation DNA damage, the cellular response to it,
and the capacity of the pathways that protect against mutation — mismatch repair, base
excision repair, nucleotide excision repair, homologous recombination, end joining,
direct alkylation reversal, crosslink repair and translesion synthesis.

**Read the design: [`docs/DNA_stress_repair_panel.md`](docs/DNA_stress_repair_panel.md)** — also typeset as [`DNA_stress_repair_panel.pdf`](docs/DNA_stress_repair_panel.pdf) (13 pages)

## The design problem

A raised 8-OHdG cannot distinguish the only two states that matter clinically: damage
production is high and repair is keeping up, or damage production is high and repair is
failing. Both give the same number. Separating them requires measuring repair in the
same subject at the same time.

A single cartridge cannot do both — 8-OHdG is a small molecule measured by mass
spectrometry, MSH2 is a protein epitope measured by antibody. The panel is therefore
modular by analyte class and matrix, with each module genuinely multiplexed internally.

| Module | Measures | Platform | Matrix | Plex |
|---|---|---|---|---|
| **A** | oxidative, alkylation, halogenation, deamination adducts | LC-MS/MS, one MRM method | urine, DNA hydrolysate | 18 analytes |
| **A2** | markers whose chemistry bars them from that method | IHC, qPCR, flow, sequencing | tissue, blood, plasma | 11 markers |
| **B** | cellular damage response | flow cytometry, one stained tube | PBMC | 10 parameters |
| **C** | repair capacity | IHC + MSI/HRD + expression array | tissue, PBMC | 9 tests + 83 targets |

## Headline design decisions

- **ELISA is excluded, not de-emphasised.** A multi-laboratory consensus exercise found
  immunoassay returns systematically higher urinary 8-oxodG than chromatographic methods
  with poor agreement between them. Every analyte carries a stable-isotope-labelled
  internal standard and is quantified by LC-MS/MS.
- **The panel can detect its own artefact.** DNA extraction manufactures 8-oxodG through
  Fenton chemistry. The FPG/OGG1-modified comet assay measures the same lesion class
  inside intact cells, so divergence between the two flags an extraction artefact. No
  single-analyte assay can do this.
- **Nine damage chemistries, not one.** A panel of ten guanine-oxidation markers would be
  ten measurements of one thing. This one spans closed-ring, ring-opened and tandem-cyclised
  purine oxidation, pyrimidine oxidation, nitration, halogenation, alkylation,
  aldehyde-derived exocyclic adducts and deamination. The chemistries point at different
  causes: halogenation implicates neutrophils, nitration implicates chronic inflammation,
  deamination implicates folate status.
- **Cyclopurines let one oxidant interrogate two repair pathways.** 8-oxodG and 8,5'-cyclo-dA
  come from the same hydroxyl radical, but the cyclopurine's C5'-C8 bond makes it a
  nucleotide-excision-repair substrate rather than a base-excision one. Normal 8-oxodG with
  elevated cdA localises the fault to NER, which no guanine-oxidation marker can do alone.
- **Mismatch repair is not the only mutation-protection pathway, and testing it alone leaves
  a hole.** Polymerase proofreading failure (POLE/POLD1 exonuclease domain) produces
  ultramutated tumours that are nonetheless *microsatellite stable*, so an MMR and MSI
  workup calls them normal. The panel treats replication fidelity as two layers and tests
  both.
- **The panel measures mutation, not only damage and machinery.** PIG-A mutant frequency and
  error-corrected duplex sequencing report fixed mutations, which is the outcome everything
  else is a proxy for.
- **Every damage analyte is coupled to the enzyme that handles it** — 8-oxodG to
  OGG1/MUTYH, O6-methyl-dG to MGMT, thymidine glycol to NTHL1, cyclopurines to the NER
  core. An abnormal composite therefore names the enzyme to follow up, not just a number.
- **The headline index is a difference of z-scores, not a ratio.** A z-scored denominator
  crosses zero by construction, so damage/repair is undefined at the population median
  and explodes near it. On the log scale a difference is exactly what a ratio reaches
  for, and it is finite everywhere.
- **Mismatch repair is the anchor** because it is where clinical standing is strongest.
  The four-antibody IHC panel is enforced with pair-aware interpretation: MSH6 loss can
  be secondary to MSH2 loss, and PMS2 loss secondary to MLH1 loss, so scoring the four
  independently over-calls primary mutation.

## Verification built into the run

`src/dna_stress_panel.py` recomputes every mass transition from molecular formula,
validates every internal standard, checks the flow panel for spectral crowding, and
confirms damage-to-repair coupling. Current state: **0 failures across 11 computable
transitions, 0 transition collisions, 14 of 14 damage analytes coupled.**

Three findings from those checks changed the design rather than being written up
afterwards. The internal-standard check initially tested *exact* mass shift against
3.0 Da and wrongly failed four `[13C1,15N2]` and `[15N3]` standards at 2.991–2.997 Da —
these are +3 *nominal* shifts, fully resolvable on a unit-resolution quadrupole. The
check now works on nominal mass and issues a bleed-through advisory at exactly +3.
Separately, two internal standards were upgraded to `[15N5]` once the check made the
+3 crowding visible. Third, the neutral-loss rule correctly *refuses* to apply to the two
cyclopurines, because the C5'-C8 bond that defines those lesions is exactly what stops the
sugar leaving as a neutral fragment. Section 8 of the report documents all three.

## Reproduce

```bash
python3 src/dna_stress_panel.py     # no dependencies beyond the standard library
python3 src/build_report_pdf.py     # typeset the report to PDF (needs reportlab)
```

## Layout

```
data/dna_damage_adducts.csv        18 Module A analytes: formula, transition, internal standard
data/dna_stress_other_markers.csv  11 Module A2 markers on other platforms
data/ddr_flow_panel.csv            10 Module B flow parameters with laser/emission assignment
data/ddr_cytogenetic_assays.csv     5 cytogenetic and comet endpoints
data/repair_clinical_tests.csv      9 clinically established repair tests with pitfalls
data/repair_expression_panel.csv   83 repair genes across 11 pathways, 2 tiers
src/dna_stress_panel.py            builder, validator and scoring framework
src/build_report_pdf.py            markdown-to-PDF renderer for the report
docs/DNA_stress_repair_panel.md    the full design document
docs/DNA_stress_repair_panel.pdf   the same report, typeset
```

## Status

No composite index here is clinically validated. The components differ enormously in
maturity: mismatch repair IHC and MSI are routine diagnostics, the micronucleus assay
has prospective cancer-incidence data in a 6718-subject cohort, and most adducts beyond
8-oxodG are research-grade. Section 9 of the report states the limitations in full.

---

# Study 2 — PROTAC architecture vs. *Primula* root constituents

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
