# A DNA-damage, DNA-oxidation and DNA-repair protein panel for plasma LC-MS/MS

> **Looking for the lymphocyte panel?** See
> [`DNA_damage_lymphocyte_panel.md`](DNA_damage_lymphocyte_panel.md). For DNA
> repair enzymes the cell is the better matrix: 60 of the 128 proteins are
> measurable in a lymphocyte lysate against 28 here, and the repair machinery
> that this page reports as unreachable is abundant there. This page remains the
> right reference for what circulates.

A target list built against UniProtKB/Swiss-Prot, digested in silico, checked for
peptide uniqueness across the whole human proteome, and tiered by whether each
protein is realistically measurable in plasma. The counts below describe the
original 118-protein list; the panel has since grown to 128 with the addition of
lymphocyte-specific targets, all of which are intracellular and tier T4 here.

## The headline you need before designing the experiment

**Most DNA repair enzymes are not plasma analytes.** Of the 118 proteins, 71 are
nuclear or cytosolic enzymes with no credible evidence of circulating freely at a
concentration any mass spectrometer will reach from plasma. That includes almost
everything in the classical repair pathways:

| Pathway | Proteins | Realistically measurable in plasma |
|---|---|---|
| Double-strand break repair | 24 | 0 |
| Nucleotide excision repair | 11 | 0 |
| Mismatch repair | 5 | 0 |
| Base excision repair | 14 | 0 |
| Oxidative base excision repair, including OGG1 | 7 | 0 |
| DNA damage response signalling | 14 | 0 |
| Nuclear damage-associated molecular patterns | 10 | 6 |
| Oxidative stress and redox | 31 | 20 |
| Stress chaperones | 2 | 2 |

A plasma experiment built on OGG1, APEX1, PARP1, ATM or BRCA1 as direct analytes
will not work. The measurable signal in plasma is of three kinds, and the panel is
organised around them.

**What is measurable is the consequence of DNA damage, not the repair machinery.**
Nucleosomal histones released by cell death and NETosis, the redox enzymes that
define the oxidative environment causing the damage, and the acute-phase iron and
heme handlers that drive Fenton chemistry.

## Tiers

Every protein carries an evidence tier in `data/plasma_evidence.tsv`, with the
reason, the expected concentration, the recommended assay and the confounder.

| Tier | n | Meaning |
|---|---|---|
| T1 abundant plasma | 5 | Classical secreted plasma proteins, µg/mL to mg/mL. Neat plasma, DIA or PRM. |
| T2 detected plasma | 23 | Documented in plasma by MS or validated immunoassay, ng/mL. Needs depletion or targeted PRM. |
| T3 EV or enriched | 19 | Only realistic via extracellular vesicle or nucleosome enrichment, or immunocapture. |
| T4 intracellular | 71 | Nuclear or cytosolic. Measure in PBMC or tumour lysate, not plasma. |

The 28 proteins in T1 and T2 are the panel you would actually run on plasma. The
T3 set is the honest extension if you add an enrichment step. The T4 set is kept
because it defines the pathway context and because those proteins are the right
targets in a PBMC or tissue arm of the same study.

## The measurable panel

Twenty-eight proteins, every one with a lead peptide whose quantification scope is
stated. Five are abundant enough for neat plasma:

| Gene | Protein | Expected level | Lead peptide |
|---|---|---|---|
| CP | Ceruloplasmin | 200-400 µg/mL | GAYPLSIEPIGVR |
| HP | Haptoglobin | 0.3-2 mg/mL | DYAEVGR |
| HPX | Hemopexin | 0.5-1 mg/mL | LYLVQGTQVYVFLTK |
| GPX3 | Plasma glutathione peroxidase | 10-30 µg/mL | YVRPGGGFVPNFQLFEK |
| SELENOP | Selenoprotein P | 3-6 µg/mL | LPTDSELAPR |

The remaining 23 sit in the ng/mL range and need depletion or a targeted method:
HMGB1, histones H4, H2B, H3.1, H3.3 and H2A, SOD1, SOD2, SOD3, CAT, GPX1, PRDX1,
PRDX2, PRDX4, PRDX6, TXN, GSTP1, GLRX, HMOX1, PARK7, MPO, HSPA1A and HSP90AA1.

## The one DNA-damage marker that does survive into plasma

H2AX is the exception worth the effort. Its C-terminal tryptic peptide **ATQASQEY**
spans residues 136-143 and is specific to H2AX across the entire human proteome.
It contains **Ser140** in UniProt numbering, which is the residue the literature
calls Ser139: the ATM, ATR and DNA-PK phosphorylation site that defines γH2AX.

That makes a phospho-PRM assay on this peptide the one route by which the canonical
DNA double-strand-break marker can be read from a blood sample rather than from
cells on a slide. The peptide ends in tyrosine, not lysine or arginine, so a
conventional C-terminal heavy Lys or Arg standard is not available and a custom
labelled standard is required. H2AX reaches plasma as nucleosomal cargo, so pair it
with nucleosome enrichment and read it against total H2A family signal.

## What the peptides actually measure

Uniqueness was tested by substring search against all 20,525 human Swiss-Prot
entries, not merely against other tryptic products. That is the conservative test:
a peptide that is a tryptic product of the target but an internal stretch of an
unrelated protein is still ambiguous under a semi-tryptic search and under the
in-vivo proteolysis that plasma proteins undergo.

- 114 of 118 proteins have a **protein-specific** lead peptide.
- 4 have only a **family-level** lead: H2BC4, H3C1, H2AC11 and H2AZ1. The core
  histone variants are near-identical, so no tryptic peptide separates one histone
  gene product from its family. The panel names the family rather than pretending
  to gene-level specificity.
- 13 proteins are **sequence-identical to another accession** in this release and
  cannot be separated by any peptide at all. HSPA1A and HSPA1B are the case most
  likely to catch you out: they are one entry, so an HSPA1A result is an
  HSPA1A-plus-HSPA1B result. The others are listed in the results table.

## Confounders you must design around

**Hemolysis dominates this panel.** Thirteen of the 28 measurable proteins are
erythrocyte-rich: SOD1, CAT, GPX1, PRDX1, PRDX2, PRDX6, TXN, GSTP1, GLRX, PARK7,
and the histone H4 signal, while HP and HPX move in the opposite direction because
they are consumed by free hemoglobin. PRDX2 is the worst case, being one of the most
abundant erythrocyte proteins. A redox panel measured on visibly or subclinically
hemolysed plasma will produce a confident, reproducible and entirely artefactual
result. Measure a hemolysis index on every sample and treat it as a covariate, not
as a sample-rejection rule only.

Three further confounders are built into the evidence table: heparin tubes raise
SOD3 and MPO by releasing them from endothelium; neutrophil activation drives both
MPO and the histone signal through NETosis; and the acute-phase response moves CP
and HP independently of any DNA damage.

## Assay design that follows from this

1. **Plasma arm, neat.** The five T1 proteins by DIA or PRM, no depletion.
2. **Plasma arm, depleted or targeted.** The 23 T2 proteins by PRM with heavy SIS
   peptides. The transition list is in `results/prm_transition_list.tsv`: 468
   transitions over 55 peptides, each with light and heavy precursor m/z, the
   required label, and six fragment ions chosen above the precursor m/z.
3. **Nucleosome arm.** Enrich circulating nucleosomes, then read the histone family
   peptides and the H2AX phosphopeptide.
4. **Cellular arm.** The 71 T4 repair enzymes in PBMC lysate from the same draw.
   This is where OGG1, APEX1, PARP1 and the rest belong.
5. **Orthogonal damage readout.** Pair the protein panel with urinary or plasma
   8-oxo-dG by LC-MS/MS. The nucleoside adduct is the direct measure of oxidative
   DNA damage; the protein panel measures the response to it.

## Reproducing

```bash
python3 src/panel_build.py plasma     # builds the panel and runs 13 internal checks
python3 src/validate.py plasma        # re-derives everything by an independent route
```

Both must exit cleanly. `src/validate.py` recomputes every peptide mass and every
fragment ion with pyteomics, re-reads every peptide coordinate from the sequence
database, re-tests proteotypicity by brute force, and compares the digestion to the
ExPASy enzymology rule.

## Method notes

**Sequence source.** UniProtKB/Swiss-Prot as distributed in the NCBI BLAST database
`swissprot`, release dated 15 September 2026, 487,728 entries, of which 20,525 are
human. The UniProt REST API is not reachable from this environment, so sequences
were taken from this mirror of the same reviewed database and every target
accession was confirmed present with its Swiss-Prot description.

**Digestion.** Trypsin, fully specific, no cleavage before proline, zero missed
cleavages for the candidate list. This is the convention used by MaxQuant, Skyline
and Spectronaut. It differs from the ExPASy enzymology rule at exactly two contexts,
W-K-P and M-R-P, which affects 10 of the 118 proteins; the validator asserts that
these are the only differences.

**Masses.** Monoisotopic, cysteine fixed at carbamidomethyl (+57.021464), with
methionine oxidation treated as a liability flag rather than a separate precursor.

**Peptide scoring.** Favours protein-specific peptides of 7 to 25 residues near
13 residues, penalises each chemical liability, penalises selenocysteine heavily,
and penalises precursors outside 600-2600 Da. Liabilities flagged: methionine and
tryptophan oxidation, cysteine alkylation, N-glycosylation sequons, deamidation and
isomerisation motifs, acid-labile Asp-Pro, ragged termini and protein termini.

**What the tiers are not.** The tier assignments are a curated reading of protein
localisation, secretion and the plasma biomarker literature. They are a design
prior, not a measurement. Confirm the T2 and T3 calls in your own matrix with
heavy peptide standards before committing to a large cohort.

## Files

| File | Contents |
|---|---|
| `data/panel_targets.tsv` | 118 targets: accession, gene, Swiss-Prot name, pathway category, length |
| `data/plasma_evidence.tsv` | Tier, expected level, basis, recommended assay and confounder per protein |
| `results/plasma_panel_proteins.tsv` | One row per protein with digestion statistics and lead peptide |
| `results/plasma_panel_peptides.tsv` | Ranked candidate peptides with m/z, scope, sharing partners and liabilities |
| `results/plasma_prm_transitions.tsv` | PRM transitions for the measurable tiers |
| `results/plasma_panel_qc.tsv` | The 13 internal checks and their results |
