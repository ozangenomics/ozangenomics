# A DNA-damage, DNA-oxidation and DNA-repair protein panel for lymphocyte LC-MS/MS

A 128-protein target list built against UniProtKB/Swiss-Prot, digested in silico,
checked for peptide uniqueness across all 20,525 human entries, and tiered by
whether each protein is realistically measurable in isolated peripheral blood
lymphocytes.

## Why the cell is the right matrix

Sixty of the 128 proteins are measurable in a lymphocyte lysate against 28 in
plasma, and the proteins that change status are precisely the ones you asked
about. Forty-one proteins that plasma cannot reach become measurable in cells,
including the entire core of the repair machinery:

> H2AX, PARP1, PRKDC, XRCC5, XRCC6, MRE11, RAD50, NBN, RPA1, RPA2, RPA3, XRCC4,
> TOP1, PCNA, NUDT1, APEX1, XRCC1, POLB, LIG3, XPC, RAD23B, ERCC2, ERCC3, DDB1,
> MLH1, MSH2, MSH6, HMGB2, DEK, NPM1, H2AZ1, GPX4, PRDX3, PRDX5, TXNRD1, TXNIP,
> GSR, G6PD, TRIM28, APOBEC3G, TFAM

Several are not merely detectable but abundant. Ku70 and Ku80 run at roughly
4×10⁵ copies per cell, PARP1 and APE1 near 10⁶, and DNA-PKcs near 10⁵ and higher
in lymphoid cells than in most tissues. These are DIA-quantifiable from a plain
whole-cell lysate, no enrichment required.

The traffic runs both ways. Seven proteins that are genuine plasma analytes are
not lymphocyte proteins at all: SOD3, GPX3, SELENOP, CP, HP, HPX and MPO. In a
cell preparation they are not analytes but **quality indicators**, and that is how
the panel uses them.

## Tiers

| Tier | n | Meaning |
|---|---|---|
| L1 abundant | 24 | High copy number. Whole-cell lysate, DIA or PRM, no enrichment. |
| L2 deep proteome | 36 | Routinely seen in a deep or fractionated lymphocyte proteome, or by PRM. |
| L3 enrichment | 44 | Needs nuclear fractionation or immunoenrichment to be quantified reliably. |
| L4 stimulation required | 14 | Near-absent in resting cells. Only appears after mitogen stimulation. |
| L4 not a lymphocyte protein | 10 | Absent. Seven are contamination markers, three are developmental. |

Coverage by pathway, counting only L1 and L2:

| Pathway | Measurable | Total |
|---|---|---|
| Nuclear damage-associated molecular patterns | 10 | 10 |
| Oxidative stress and redox | 18 | 31 |
| Double-strand break repair | 11 | 24 |
| Base excision repair | 5 | 14 |
| Nucleotide excision repair | 5 | 11 |
| Mismatch repair | 3 | 5 |
| DNA damage response signalling | 3 | 15 |
| Oxidative base excision repair | 1 | 7 |
| Stress chaperones | 2 | 2 |
| Mitochondrial DNA maintenance | 1 | 2 |
| Lymphocyte deaminases | 1 | 4 |
| Lymphocyte programmed breaks | 0 | 3 |

## The three things that will decide whether this experiment works

### 1. Resting lymphocytes are in G0, and that silences a whole arm of repair

Fourteen proteins are effectively absent from quiescent cells and appear only on
mitogen stimulation: **RAD51, BRCA1, BRCA2, FANCD2, BLM, TOP2A, CHEK1, NEIL3,
RAD52, MDM2, GADD45A, SESN2, NFE2L2 and AICDA**.

This is not a detection limit, it is biology. Homologous recombination needs a
sister chromatid, so a G0 lymphocyte does not run it. If your design compares
donors or exposures and some samples have been stimulated or contain more
activated cells, you will measure proliferation and report it as repair capacity.

Decide explicitly which experiment you are running. Resting cells give you the
constitutive, non-homologous end joining and base excision repair machinery.
PHA or anti-CD3 stimulation opens up homologous recombination but changes the
redox proteome at the same time: SOD2 rises, TXNIP falls sharply. Both are
legitimate; measuring a mixture of the two is not.

PCNA and the H3.3 to H3.1 ratio both track proliferation and are already in the
panel. Use them as internal proliferation indices rather than assuming your
samples are comparable.

### 2. The isolation procedure itself damages the cells

Density-gradient separation, ambient temperature and time to processing all
induce oxidative stress and a DNA damage response before the sample ever reaches
the mass spectrometer. HSPA1A and HSP90AA1 are in the panel partly as handling
stress reporters.

The practical consequence is that time from venepuncture to lysis has to be
fixed, logged and treated as a covariate. An ex vivo artefact of this kind is
systematic, not random, and will track with anything that correlates with
handling, including clinic, shift and sample batch.

### 3. Purity is a measurement, not an assumption

The panel contains its own purity controls:

- **MPO** is a neutrophil azurophil granule protein. Any MPO signal means
  granulocyte contamination. It is the single most useful purity marker here.
- **CP, HP, HPX, SELENOP, GPX3 and SOD3** are secreted plasma proteins. Their
  signal is plasma carryover and reports wash quality.
- **CAT, HMOX1, NQO1 and APOBEC3A** are monocyte-biased. If you work from PBMC
  rather than purified lymphocytes, these report monocyte content as much as
  lymphocyte redox state.
- **PRDX2, SOD1, PARK7, G6PD, GPX1, GSR and GLRX** are erythrocyte-rich, so
  incomplete red cell lysis inflates them.

Report these alongside the analytes. A redox result from a preparation with
variable granulocyte and erythrocyte content is not interpretable.

## OGG1 and the oxidative glycosylases

The oxidative base excision repair enzymes are the weakest part of the panel and
you should know that before designing around them. OGG1 sits at roughly 10³ to
10⁴ copies per cell. It is genuinely low copy, and so are MUTYH, NTHL1, NEIL1 and
NEIL2. All are L3: nuclear fractionation or immunoenrichment, then PRM with a
heavy standard. NEIL3 is worse, being proliferation-restricted as well.

OGG1 does have a clean proteotypic peptide, **LDLVLPSGQSFR**, so the assay is
buildable. But a nuclear fraction and a heavy peptide standard are not optional
for it, and protein abundance is a poor proxy for repair flux in any case. Pair
it with a functional or lesion measure rather than relying on abundance alone.

## The flagship assay: γH2AX by PRM

In plasma, H2AX is a difficult target reachable only through circulating
nucleosomes. In lymphocytes it is a **core histone, present at 2 to 25 percent of
the H2A pool**, and the assay becomes straightforward.

Two peptides are specific to H2AX across the entire human proteome:

| Peptide | Residues | [M+2H]²⁺ | Use |
|---|---|---|---|
| TSATVGPK | 121-128 | 380.7136 | total H2AX, the denominator |
| ATQASQEY | 136-143 | 449.2011 | the C-terminal tail carrying Ser140 |

ATQASQEY contains **Ser140** in UniProt numbering, which is the residue the
literature calls Ser139: the ATM, ATR and DNA-PK site that defines γH2AX.
Measuring the phosphorylated and unphosphorylated forms of this peptide against
TSATVGPK gives a **stoichiometric** γH2AX occupancy, which is what foci counting
and western blotting only approximate.

Two practical points. The peptide ends in tyrosine, so a conventional C-terminal
heavy Lys or Arg standard does not exist and you need a custom labelled peptide.
And phosphopeptide recovery requires enrichment, so spike the heavy standard
before enrichment, not after, or the ratio will report recovery rather than
biology.

### Other phosphosites, and which of them trypsin can actually reach

`results/lymphocyte_phospho_markers.tsv` checks eleven canonical damage-response
phosphosites against the sequence and asks whether a usable tryptic peptide
contains them. The answer is often no, which is worth knowing before you plan the
assay rather than after.

| Site | Marker | Tryptic peptide | Verdict |
|---|---|---|---|
| H2AX Ser140 | γH2AX | ATQASQEY, 8 residues | good PRM candidate |
| CDKN1A Ser146 | p21 | QTSMTDFYHSK, 11 | good PRM candidate |
| NBN Ser343 | NBS1, ATM substrate | TTTPGPSLSQGVSVDEK, 17 | good PRM candidate |
| ATM Ser1981 | ATM autophosphorylation | SLAFEEGSQSTTISSLSEK, 19 | good PRM candidate |
| TP53 Ser15 | p53, ATM and ATR | MEEPQSDPSVEPPLSQETFSDLWK, 24 | workable, Met oxidation liability |
| TRIM28 Ser824 | KAP1, ATM substrate | 31 residues | too long; use Glu-C or Asp-N |
| CHEK1 Ser345 | CHK1, ATR substrate | 35 residues | too long; use Glu-C or Asp-N |
| RPA2 Ser4 and Ser33 | DNA-PK and ATR | 37 residues, 12 S/T/Y | too long; use Glu-C or Asp-N |
| TP53BP1 Ser25 | 53BP1, ATM substrate | 45 residues | too long; use Glu-C or Asp-N |
| CHEK2 Thr68 | CHK2, ATM substrate | none between 6 and 45 residues | trypsin cannot reach it |

Every one of these peptides carries more than one serine, threonine or tyrosine,
so site localisation needs good fragmentation. Do not assume a mass shift on the
precursor means the site you intended.

KAP1 Ser824 remains the best second readout biologically, being an abundant and
direct ATM substrate that is not a histone, but it needs a protease other than
trypsin. CHK2 Thr68 is the clearest case where an immunoassay beats targeted MS.

## Lymphocyte-specific DNA damage biology

Lymphocytes are the only cells that break their own DNA on purpose, so the panel
includes the enzymes that do it. In peripheral blood the finding is reassuring:

- **RAG1, RAG2 and DNTT** drive V(D)J recombination in developing thymocytes and
  marrow B cells. They are absent from mature peripheral lymphocytes, so they do
  not confound a peripheral blood study. They would dominate one in thymus,
  marrow or acute lymphoblastic leukaemia blasts, and DNTT is a standard marker
  there.
- **AICDA (AID)** deaminates cytosine during class-switch recombination and
  somatic hypermutation. It is absent from resting peripheral B cells and induced
  in germinal-centre ones, so it matters in tonsil or lymph node, not in blood.
- **APOBEC3A, APOBEC3B and APOBEC3G** are expressed in peripheral blood and are a
  real source of genomic uracil. APOBEC3G is the best measurable of the three.
  APOBEC3A and APOBEC3B share a catalytic domain, so check the sharing column
  before trusting a peptide to distinguish them.

**UNG** is in the panel for the same reason: it is the glycosylase that processes
AID-generated uracil, so in B-cell work it reads out programmed rather than
oxidative damage.

Mitochondrial DNA is covered by **TFAM** and **POLG**. Mitochondrial DNA carries a
higher oxidative lesion burden than nuclear DNA and is often the more sensitive
compartment, but TFAM also tracks mitochondrial content, which itself rises on
activation.

## What the peptides actually measure

Uniqueness was tested by substring search against all 20,525 human Swiss-Prot
entries, not merely against other tryptic products. That is the conservative
test: a peptide that is a tryptic product of the target but an internal stretch
of an unrelated protein is still ambiguous under a semi-tryptic search.

- 124 of 128 proteins have a **protein-specific** lead peptide.
- 4 have only a **family-level** lead: H2BC4, H3C1, H2AC11 and H2AZ1. Core
  histone variants are near-identical and no tryptic peptide separates one
  histone gene product from its family, so the panel names the family.
  H2AX is the exception that makes the γH2AX assay possible.
- 13 proteins are **sequence-identical to another accession** and cannot be
  separated by any peptide. HSPA1A and HSPA1B are one entry, so an HSPA1A result
  is an HSPA1A-plus-HSPA1B result.

## Suggested design

1. **Purified lymphocytes, not PBMC**, if the redox arm matters, because
   monocytes dominate CAT, HMOX1 and NQO1. Record the purity markers either way.
2. **Fix and log the time from draw to lysis.** Treat it as a covariate.
3. **Whole-cell lysate, DIA.** Covers the 24 L1 and much of the 36 L2 set.
4. **Nuclear fraction, PRM with heavy standards.** The only route to OGG1 and
   the other glycosylases. `results/lymphocyte_prm_transitions.tsv` gives 1,098
   transitions over 119 peptides and 60 proteins, each with light and heavy
   precursor m/z, the required label and six fragment ions above the precursor.
5. **Phospho-enrichment arm.** γH2AX Ser140 against total H2AX, plus KAP1 Ser824.
   Spike heavy standards before enrichment.
6. **Paired stimulation arm** if you want the homologous recombination proteins.
   Report resting and stimulated separately, never pooled.
7. **Orthogonal lesion measure.** The comet assay or 8-oxo-dG by LC-MS/MS on DNA
   from the same cells. Protein abundance is not repair flux, and this panel
   measures abundance.

## Reproducing

```bash
python3 src/panel_build.py lymphocyte && python3 src/validate.py lymphocyte
python3 src/panel_build.py plasma     && python3 src/validate.py plasma
python3 src/build_workbook.py         # the Excel workbook
```

The workbook carries no formulas, only values. LibreOffice could not be run in
the build environment, so a formula could not be recalculated and would have
read back as blank to pandas and most previewers; the tier counts on its
Overview sheet are therefore fixed at build time, and filtering the All markers
sheet on peptide rank 1 reproduces them.

Both must exit cleanly. `src/validate.py` recomputes every peptide mass and
fragment ion with pyteomics, re-reads every coordinate from the sequence
database, re-tests proteotypicity by brute force, and compares the digestion
against the ExPASy enzymology rule.

## Method notes

**Sequence source.** UniProtKB/Swiss-Prot as distributed in the NCBI BLAST
database `swissprot`, release dated 15 September 2026, of which 20,525 entries
are human. Rebuild with `src/fetch_swissprot.py`.

**Digestion.** Trypsin, fully specific, no cleavage before proline, zero missed
cleavages. This is the MaxQuant, Skyline and Spectronaut convention. It differs
from the ExPASy enzymology rule only at W-K-P and M-R-P, affecting 10 of the 128
proteins; the validator asserts those are the only differences.

**Masses.** Monoisotopic, cysteine fixed at carbamidomethyl, methionine
oxidation treated as a liability flag rather than a separate precursor.

**What the tiers are not.** Tier assignments are a curated reading of lymphocyte
expression, subcellular localisation and cell-cycle dependence. They are a design
prior, not a measurement. Confirm the L2 and L3 calls in your own preparation
with heavy standards before committing to a cohort. Copy-number figures are
order-of-magnitude orientation from published cell-line and primary-cell
proteomics, not values measured here.

## Files

| File | Contents |
|---|---|
| `data/panel_targets.tsv` | 128 targets: accession, gene, Swiss-Prot name, pathway, length |
| `data/lymphocyte_evidence.tsv` | Tier, expected level, basis, assay and confounder per protein |
| `data/plasma_evidence.tsv` | The same for plasma, for comparison |
| `data/human_swissprot.fasta.gz` | The 20,525 human Swiss-Prot entries used |
| `results/lymphocyte_panel_proteins.tsv` | One row per protein with digestion stats and lead peptide |
| `results/lymphocyte_panel_peptides.tsv` | 635 ranked candidate peptides |
| `results/lymphocyte_prm_transitions.tsv` | 1,098 PRM transitions for the measurable tiers |
| `results/lymphocyte_panel_qc.tsv` | The 13 internal checks |
| `results/lymphocyte_marker_list.tsv` | **The flat marker list**: every protein with its top three peptides, tier, m/z, scope and confounder |
| `results/lymphocyte_phospho_markers.tsv` | Canonical DDR phosphosites, verified against the sequence, with phospho precursor m/z and a verdict on tryptic suitability |
| `results/DNA_damage_lymphocyte_panel.xlsx` | **The whole panel as one Excel workbook**, six sheets, built by `src/build_workbook.py` |
