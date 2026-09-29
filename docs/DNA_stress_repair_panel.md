# A multiplex panel for DNA stress and DNA repair capacity

**Design document. Built outward from 8-OHdG to a panel that quantifies oxidative and
alkylation damage, the cellular response to it, and the capacity of the pathways that
protect against mutation.**

---

## 0. Executive summary

8-OHdG (properly 8-oxo-7,8-dihydro-2'-deoxyguanosine, **8-oxodG**) is the most used
oxidative DNA damage marker in the literature. It is also the single most
artefact-prone analyte in the field, and on its own it answers the wrong question.
A raised 8-oxodG cannot distinguish the two states a clinician actually cares about:

- **damage production is high, and repair is keeping up** — a flux state, largely benign
- **damage production is high, and repair is failing** — the mutagenic state

Both give the same 8-oxodG number. Separating them requires measuring repair in the
same subject at the same time. That is the design principle of this panel.

**The panel has three modules**, because a single platform cannot carry both halves of
the question. 8-oxodG is a small molecule quantified by mass spectrometry; MSH2 is a
protein epitope quantified by antibody. No cartridge multiplexes those together. Each
module is, however, genuinely multiplexed internally:

| Module | What it measures | Platform | Matrix | Plex |
|---|---|---|---|---|
| **A** | oxidative, alkylation, halogenation and deamination adducts | LC-MS/MS, one MRM method | urine, DNA hydrolysate | 18 analytes |
| **A2** | stress markers whose chemistry bars them from the MRM method | IHC, qPCR, flow, sequencing | tissue, blood, plasma | 11 markers |
| **B** | cellular damage response | flow cytometry, one stained tube | PBMC | 10 parameters |
| **C** | repair capacity | IHC + MSI/HRD + expression array | tissue, PBMC | 9 tests + 83 targets |

**Headline readout** is the Stress-to-Repair Index, defined in Section 5 as a
difference of robust z-scores rather than a ratio, and reported *pathway-resolved* so
that an abnormal result names the enzyme to follow up rather than just flagging a
number.

**Verification is built into the build.** `src/dna_stress_panel.py` recomputes every
mass-spectrometric transition from the molecular formula, checks every internal
standard for a resolvable mass shift and an intact elemental backbone, tests the flow
panel for spectral crowding, and confirms that every damage analyte has a named
enzymatic partner in the repair module. It currently reports **0 failures across
11 computable transitions and 14 of 14 damage analytes coupled**. Three findings in
Section 8 were produced by those checks and changed the design.

---

## 1. Why 8-OHdG needs company: the marker's three real problems

### 1.1 Measurement method changes the answer

This is not a minor calibration issue. A multi-laboratory comparison organised as a
consensus exercise on urinary 8-oxodG found that **immunoassay (ELISA) returns
systematically higher values than chromatographic methods, and the agreement between
the two is poor** — a difference large enough that studies using different platforms are
not comparable. The consensus recommendation was chromatographic quantification with
isotope-labelled internal standards.

**Design consequence.** This panel specifies **LC-MS/MS with a stable-isotope-labelled
internal standard for every analyte**. ELISA is excluded, not de-emphasised. Antibody
cross-reactivity against structurally related urinary components is the accepted
explanation for the positive bias, and no amount of local calibration fixes it.

### 1.2 The extraction step manufactures the analyte

Guanine oxidises during DNA isolation. Phenol-based extraction, iron contamination and
air exposure all generate 8-oxodG that was not present in the cell. Reported tissue
8-oxodG levels in the older literature span more than an order of magnitude, and the
spread is attributed largely to extraction chemistry rather than biology.

**Design consequence.** Section 6 specifies chelator-supplemented, phenol-free
extraction, and the panel carries an **independent cross-check**: the
FPG/OGG1-modified comet assay measures the same lesion class *inside intact cells*,
where no extraction has occurred. Divergence between the LC-MS/MS result and the
enzyme-modified comet result is the strongest available signal that the adduct number
is an extraction artefact. No single-analyte assay can detect its own artefact this way.

### 1.3 Urinary and genomic 8-oxodG mean different things

Urinary 8-oxodG is excreted repair product. It is a **rate** — whole-body lesion
turnover. Genomic 8-oxodG in DNA hydrolysate is a **steady-state level** — what is
present right now. High repair activity *raises* the urinary number while *lowering* the
genomic one.

**Design consequence.** The panel measures 8-oxodG in both matrices and never averages
them. Their relationship is itself informative: high urinary with low genomic is
competent repair under load; low urinary with high genomic is repair failure.

---

## 2. Module A — the adduct panel

Eighteen analytes in one liquid-chromatography tandem-mass-spectrometry method. All
transitions below were **recomputed from molecular formula** by the build script, not
transcribed. For 2'-deoxynucleosides the product ion is the protonated base, formed by
loss of the intact sugar (C5H8O3, 116.0473 Da); for the ribonucleoside the loss is
C5H8O4 (132.0423 Da).

### 2.1 Chemistry, and why each analyte is in the panel

| Analyte | Lesion chemistry | Precursor → product | Why it earns a channel |
|---|---|---|---|
| **8-oxodG** | guanine C8 oxidation, DNA | 284.10 → 168.05 | The anchor. Most characterised, most comparable to prior literature |
| **8-oxoGuo** | guanine C8 oxidation, **RNA** | 300.09 → 168.05 | Independent axis. Urinary RNA-oxidation markers were associated with long-term mortality in newly diagnosed type 2 diabetes. Rises with mitochondrial and metabolic stress even when DNA oxidation is normal |
| **8-oxoGua** | free base | 168.05 → 140.04 | Excision-flux reporter. Released both by glycosylase action and by nucleotide-pool sanitisation. Rising alongside 8-oxodG means excision is *working* |
| **FapyGua** | ring-**opened** purine | 170.07 → 152.06 | Forms preferentially at **low oxygen tension**, where 8-oxoG formation is disfavoured. Omitting it systematically under-reports purine oxidation in hypoxic tissue |
| **Tg** | pyrimidine ring saturation | 277.10 → 161.06 | Substrate of NTHL1/NEIL1, **not** OGG1. Independent read on the non-OGG1 arm of base excision repair |
| **5-hmdC** | cytosine oxidation | 258.11 → 142.06 | Pyrimidine axis — but see the trap in 2.2 |
| **5-medC** | methylation (not damage) | 242.11 → 126.07 | Not a lesion. Present solely to make 5-hmdC interpretable |
| **O6-medG** | O6 alkylation | 282.12 → 166.07 | The specific substrate of MGMT. The tightest analyte-to-enzyme pairing in the panel |
| **N7-meGua** | N7 alkylation, free base | 166.07 → 149.05 | Most abundant alkylation adduct; depurinates spontaneously so it appears in urine. Alkylating-exposure dosimeter |
| **M1dG** | malondialdehyde exocyclic adduct | 304.10 → 188.06 | Converts a lipid-peroxidation signal into a covalent DNA lesion |
| **εdA** | etheno-adenine | 276.11 → 160.06 | The only adenine-directed lesion; reports aldehyde and vinyl-monomer exposure |
| **dG** | unmodified | 268.10 → 152.06 | Denominator. Lesions reported per 10⁶ dG. Also confirms hydrolysis went to completion |
| **Creatinine** | normaliser | 114.07 → 44.05 | Mandatory for spot urine |
| **8-iso-PGF2α** | lipid peroxidation | 353.23 → 193.12 (negative) | Separates systemic oxidative stress from DNA-specific damage. Pooled reference IQR 0.18–0.40 µg/g creatinine by chemical methods in adults with BMI under 25 |
| **cdA** | 8,5'-cyclo-dA, tandem cyclopurine | 250.09 → 164.06 | Oxidative lesion repaired by **NER, not BER**. See 2.5 |
| **cdG** | 8,5'-cyclo-dG | 266.09 → 180.05 | Guanine counterpart; R and S diastereomers must be resolved |
| **5-Cl-dC** | chlorinated cytosine | 262.06 → 146.01 | Myeloperoxidase/HOCl product. Neutrophil-driven inflammation, chemically distinct from oxidation and nitration |
| **dU** | uracil in DNA | 229.08 → 113.03 | Cytosine deamination and dUTP misincorporation. A third damage mechanism: neither oxidative nor alkylative |

### 2.2 The 5-hmdC trap

5-hydroxymethyl-dC is **both** an oxidation product of cytosine **and** a programmed
epigenetic mark written deliberately by TET enzymes. Scoring it as damage without
reference to its precursor pool will misread normal epigenetic remodelling as oxidative
injury. This is why 5-medC occupies a channel: 5-hmdC is reported as a fraction of
5-medC, and the panel carries TDG and SMUG1 in the expression module so that programmed
turnover can be distinguished from oxidative chemistry.

### 2.3 Why the panel is chemically non-redundant

The analytes are deliberately spread across **different base targets and different
oxidation chemistries**, because reactive species have different preferences:

- purine oxidation, closed-ring: 8-oxodG, 8-oxoGuo, 8-oxoGua
- purine oxidation, ring-opened, oxygen-tension-dependent: FapyGua
- pyrimidine oxidation: 5-hmdC, Tg
- alkylation: O6-medG, N7-meGua
- aldehyde-derived exocyclic adducts: M1dG, εdA
- adenine-directed: εdA

A panel of ten guanine-oxidation markers would be ten measurements of one thing. This
one spans four lesion chemistries on three bases across DNA and RNA.

### 2.4 Two flagged ion overlaps

The build script identified two pairs sharing a product ion:

- **8-oxodG and 8-oxoGuo both fragment to 168.05.** This is expected and correct: both
  release the same protonated 8-oxoguanine base. Precursors differ by 16 Da, so the
  pair is fully resolved in the first quadrupole.
- **FapyGua and dG both report a product near 152.06.** Precursors are 170.07 and
  268.10 — resolved, but retention-time separation should be confirmed since FapyGua is
  the ring-opened analogue of the same nucleobase.

Neither is a collision. No analyte pair shares *both* precursor and product, which the
script verifies explicitly.

### 2.5 The cyclopurines, and why they are the most useful addition

8,5'-cyclo-2'-deoxyadenosine and its guanine counterpart form when a hydroxyl radical
attacks the sugar and the resulting radical cyclises between C5' and C8. That covalent
bond between sugar and base is the whole point:

- It **distorts the duplex**, so the lesion is recognised as bulky helix-damaging injury.
- It means these lesions are **repaired by nucleotide excision repair rather than base
  excision repair**, unlike every other oxidative lesion in the panel.
- It also means **the glycosidic bond cannot be cleaved**. The build script's neutral-loss
  check correctly refuses to compute a product ion for these two analytes for exactly
  this reason, and says so in its output.

The design value is considerable. 8-oxodG and cdA are produced by the **same oxidant**
but cleared by **different pathways**. Measuring both lets a single oxidative stress
interrogate two repair systems at once. If 8-oxodG is normal while cdA is elevated, the
oxidant load is ordinary and the NER arm is the problem — a conclusion no
guanine-oxidation marker can reach on its own. These lesions block polymerases, inhibit
gene expression, and accumulate in tissue in association with disease.

### 2.6 Four damage chemistries, not one

With these additions the adduct module spans mechanistically distinct injury types,
each with its own biological source:

| Chemistry | Analytes | Source |
|---|---|---|
| Oxidation, closed-ring purine | 8-oxodG, 8-oxoGuo, 8-oxoGua | general reactive oxygen species |
| Oxidation, ring-opened purine | FapyGua | ROS at low oxygen tension |
| Oxidation, tandem cyclisation | cdA, cdG | hydroxyl radical specifically |
| Oxidation, pyrimidine | 5-hmdC, Tg | ROS |
| **Nitration** | 8-nitroguanine (Module A2) | peroxynitrite, chronic inflammation |
| **Halogenation** | 5-Cl-dC | myeloperoxidase/HOCl, neutrophils |
| Alkylation | O6-medG, N7-meGua | nitrosamines, tobacco, chemotherapy |
| Aldehyde exocyclic | M1dG, εdA | lipid peroxidation |
| **Deamination** | dU | spontaneous, folate deficiency, antifolates |

This matters because the chemistries point at different causes. Halogenation implicates
neutrophils. Nitration implicates peroxynitrite and chronic inflammation. Deamination
implicates folate status or thymidylate-synthase inhibition. Alkylation implicates
exposure. A panel confined to guanine oxidation cannot distinguish any of these.

---

## 2A. Module A2 — markers that cannot go in the mass-spectrometry method

Eleven markers in `data/dna_stress_other_markers.csv`. Each is excluded from Module A
for a **chemical or compartmental reason**, not an arbitrary one, and several are as
clinically familiar as 8-OHdG itself.

### 2A.1 8-nitroguanine — the nitrative axis

Nitration is a damage chemistry the adduct panel otherwise misses entirely. 8-nitroguanine
is a mutagenic lesion formed during chronic inflammation via reactive nitrogen species.
It has been demonstrated at sites of carcinogenesis, accumulates in cancer-prone
inflammatory disease caused by several pathogens including human papillomavirus and
Epstein-Barr virus, and strong formation in tumour tissue has been closely associated
with poor prognosis. Given that chronic inflammation is estimated to account for
roughly a quarter of cancer cases, an inflammation-specific DNA lesion earns its place.

**Why it is not in the MRM method:** 8-nitroguanine is chemically labile and depurinates
rapidly, so immunohistochemistry with a specific antibody is the established detection
route. The marker's own instability dictates its platform.

### 2A.2 Abasic sites — measuring the bottleneck directly

An abasic site is the obligate intermediate between glycosylase excision and gap
filling, and it is **more mutagenic than the base lesion it replaced**. Section 4.2
argues that high glycosylase activity with low APE1 is worse than low glycosylase
activity; abasic site burden is how you measure that state instead of inferring it.
Quantified by aldehyde-reactive probe, with methoxyamine binding abasic sites covalently.

### 2A.3 The mitochondrial genome — a separate compartment

Three markers cover it: copy number, the 4977 bp common deletion, and polymerase-blocking
lesion frequency by long-amplicon qPCR. Mitochondrial DNA **lacks nucleotide excision
repair** and sits adjacent to the electron transport chain, so it accumulates oxidative
damage faster than nuclear DNA. The common deletion has been reported at higher levels
in the blood of breast cancer patients with concurrently lower mitochondrial DNA
content and higher oxidative damage markers.

A nuclear DNA hydrolysate does not report this compartment at all. Copy number is read
only alongside deletion burden, because copy number alone moves with both mitochondrial
biogenesis and mitochondrial loss and is ambiguous by itself.

### 2A.4 Telomere length — a different time constant

Telomeric DNA is guanine-rich, disproportionately sensitive to oxidation, and refractory
to repair. Psychological stress has been associated with higher oxidative stress, lower
telomerase activity and shorter telomeres, with the highest-stress group shorter by the
equivalent of at least a decade of additional ageing.

Its role here is **temporal**. Adducts integrate hours to days. Telomere length
integrates months to years. Including both separates an acute insult from a chronic
burden, which no single-timescale marker can do. The panel pairs it with TERT and TERF2
transcripts, since telomerase activity and telomere length answer different questions.

### 2A.5 Mutation outcome — what the whole panel is actually about

Every other marker is a lesion or an enzyme. These two are fixed mutations:

- **PIG-A mutant frequency**, by flow cytometry for loss of GPI-anchored proteins on
  erythrocytes. Reported background is roughly 2.9 to 5.56 per million mutant
  erythrocytes in healthy subjects, with increases after cigarette smoking, radiotherapy
  and occupational exposures including lead. A low, stable background is what makes
  induced mutation detectable.
- **Error-corrected (duplex) sequencing**, which returns mutation frequency *and*
  spectrum — an outcome measure and a mechanism attribution in one assay. Inter-individual
  variability in clonally expanded mutations is a known analytical challenge.

A panel that measures damage and repair but never mutation is measuring inputs and
machinery while declining to check the output.

---

---

## 3. Module B — the cellular damage response

Ten parameters on one stained tube, three lasers. The channel plan is in
`data/ddr_flow_panel.csv`; the spectral check output is in
`results/moduleB_spectral_check.csv`.

### 3.1 Two channels are not optional

**Cleaved caspase-3 (BV605).** Apoptotic cells generate pan-nuclear γH2AX from
programmed DNA fragmentation. That signal is not repair-relevant damage, and a panel
without an apoptosis gate reports cell death as genotoxic stress. Work correlating
γH2AX with radiation-induced apoptosis in lymphocytes found a **negative** correlation
between residual γH2AX and apoptosis at 48 hours, interpreted as individuals less
efficient at clearing damaged cells retaining higher basal and residual damage. Whether
a cell is dying or repairing changes the meaning of its γH2AX entirely.

**DNA content (DAPI).** γH2AX and RAD51 are strongly S/G2-dependent. Uncorrected foci
counts track proliferation rather than damage. Every foci measurement in this panel is
reported cell-cycle-resolved.

### 3.2 The damage and repair channels

| Target | Reads | Fluorochrome | Why |
|---|---|---|---|
| γH2AX (Ser139) | double-strand breaks | AF488 | Best-characterised cellular damage marker; used as a pharmacodynamic endpoint in trials |
| 53BP1 foci | DSB recognition | PE | Co-localisation with γH2AX upgrades a signal to a genuine break |
| pATM (Ser1981) | apical DSB kinase | PE-Dazzle594 | Distinguishes damage being *signalled* from damage being ignored |
| pCHK2 (Thr68) | checkpoint arm | PerCP-Cy5.5 | High γH2AX with low pCHK2 is checkpoint escape — a pro-mutagenic state |
| pKAP1 (Ser824) | DSB-specific signalling | APC-Fire750 | More DSB-specific than γH2AX; not induced by replication stress alone |
| RAD51 foci | HR execution | APC | Functional HR, as distinct from the genomic HRD scar |
| PAR polymer | PARP1 activity | AF700 | Catalytic output, the functional counterpart of PARP1 transcript |
| Ki-67 | proliferation covariate | BV711 | A damaged quiescent cell and a damaged cycling cell mean different things |

γH2AX is **sensitive but not specific** — it rises in replication stress and apoptosis
as well as at true breaks. The panel therefore never reads it alone: 53BP1
co-localisation, pKAP1 specificity, and the caspase-3 gate are what make a γH2AX number
interpretable.

### 3.3 Spectral finding

The script flagged one spillover-critical adjacency: **53BP1 (PE, 578 nm) against pATM
(PE-Dazzle594, 610 nm)**, a 32 nm gap on the same 488 nm laser. This pair requires
single-stain compensation controls and, if the instrument allows, reassignment of pATM
to a red-laser channel. Reported as a build finding rather than silently accepted.

### 3.4 Cytogenetic endpoints — the only prospectively validated damage measure

Run from the same blood draw, separate culture. These are in
`data/ddr_cytogenetic_assays.csv`.

The **cytokinesis-block micronucleus assay** carries something no other endpoint in this
panel has: **prospective validation against cancer incidence in humans.** In a
6718-subject cohort across ten countries and twenty laboratories, subjects in the
medium and high micronucleus-frequency groups showed increased all-cancer incidence
(relative risk 1.84, 95% CI 1.28–2.66 for medium; 1.53, 95% CI 1.04–2.25 for high) and
reduced cancer-free survival. The association held across national cohorts and major
cancer sites.

Two design points follow. First, the micronucleus assay is the **anchor for clinical
meaning** in the whole panel; the adduct measurements are mechanism, this is outcome.
Second, that study standardised subjects by **percentile of the distribution within each
laboratory**, precisely because inter-laboratory variability is large. Any deployment of
this panel must band micronucleus results against its own laboratory distribution, not
against a published absolute threshold.

The **alkaline comet assay** likewise has pooled international baseline frequencies and
characterised confounders available from large database-comparison efforts; use those
for banding rather than inventing local cut-offs.

The **ex-vivo challenge comet assay** is what converts the panel from a snapshot to a
capacity measurement: give PBMCs a calibrated dose, measure residual damage at 2 and
24 hours, and read repair kinetics directly. Residual damage at 24 hours is the
functional deficiency readout.

---

## 4. Module C — repair capacity, with mismatch repair as the anchor

### 4.1 The clinically established tests

Nine tests in `data/repair_clinical_tests.csv`. These are the parts of the panel that
already have regulatory and guideline standing, and mismatch repair is where that
standing is strongest.

**MMR protein immunohistochemistry — MLH1, MSH2, MSH6, PMS2.** Loss of nuclear staining
with retained internal control indicates deficient MMR. This four-antibody panel is
routine diagnostic practice, embedded both in Lynch-syndrome screening and in
checkpoint-inhibitor eligibility.

> **The heterodimer rule, which the panel enforces.** These four proteins work as two
> obligate pairs. MSH2 partners MSH6; MLH1 partners PMS2. **MSH6 loss can be secondary
> to MSH2 loss, and PMS2 loss secondary to MLH1 loss.** Scoring the four markers
> independently over-calls primary mutation events. The panel specifies pair-aware
> interpretation. Expression data reflecting MMR-protein status has been examined in
> checkpoint-inhibitor-treated cohorts outside colorectal cancer, which is a reminder
> that the four-marker readout travels beyond its original tumour setting and needs
> its interpretation rules to travel with it.

**MSI by PCR** (BAT-25, BAT-26, NR-21, NR-24, MONO-27) is the orthogonal confirmation.
It carries a specific documented trap: **minimal microsatellite shift occurs in
MSI-high tumours — notably endometrial — and is a recognised cause of false-negative
interpretation.** The panel specifies matched normal tissue and an explicit analyst
alert for small shifts, rather than relying on automated peak calling.

**MGMT promoter methylation** is the alkylation-repair readout and is predictive for
alkylating chemotherapy benefit, most established in glioblastoma where temozolomide
has been compared against radiotherapy regimens in randomised trials. The panel
requires the **quantitative** methylation value, not only the binary call, because the
methylated/unmethylated cut-off is assay-dependent and the intermediate zone is
genuinely ambiguous.

**HRD score** (loss of heterozygosity + telomeric allelic imbalance + large-scale state
transitions) reports homologous-recombination status for PARP-inhibitor decisions. The
important limitation, and the reason Module B carries RAD51 foci: **a genomic scar is a
historical record, not a statement about present-day function.** A tumour can carry a
high HRD score and have restored RAD51 loading — that is a PARP-inhibitor resistance
phenotype, and the scar score alone will misclassify it. Available HRD assays differ
and each requires validation on clinical samples before clinical use.

**POLE / POLD1 exonuclease domain sequencing** covers the mutation-protection pathway
that mismatch repair testing cannot see, and it is the single most important addition to
the repair module.

> **Why an MMR-only panel is not enough.** Polymerases ε and δ replicate the genome and
> proofread their own work through exonuclease domains. Exonuclease domain mutations
> produce an **ultramutated** tumour — sometimes exceeding a million base substitutions —
> that is nonetheless **apparently microsatellite stable**. An MMR immunohistochemistry
> and MSI workup reports such a tumour as normal. The most extreme mutator phenotype in
> human cancer is invisible to the tests that define the mismatch repair readout.
> Germline exonuclease domain mutations in these genes define polymerase
> proofreading-associated polyposis, with high penetrance and dominant inheritance.
> Proofreading and mismatch repair defects can also co-occur, and that combination
> produces the highest mutation rates observed — somatic POLE exonuclease mutations
> arising on a background of biallelic mismatch repair deficiency have been reported with
> tumour mutation rates of 237 and 123 per megabase.

Framing the problem as mismatch repair *or other mutation-protection pathways* is the
right framing, and proofreading is the pathway it points at. The panel therefore treats
replication fidelity as a **two-layer system**: the polymerase proofreads as it goes, and
mismatch repair corrects what proofreading missed. Testing only the second layer leaves
the first unexamined, and the first is where the largest mutation burdens come from.
POLE, POLD1 and PCNA are in the expression panel under pathway code `PROOF`.

**Mutational signature deconvolution** is the only element that attributes observed
mutations to a *specific failed pathway* — MMR deficiency, HR deficiency, oxidative
damage and alkylation each leave distinguishable single-base-substitution signatures.
Tumour mutational burden, by contrast, is a consequence of accumulated repair failure
and is mechanistically uninformative on its own.

### 4.2 The multiplexed expression panel

83 targets across eleven repair pathways plus six normaliser candidates, in
`data/repair_expression_panel.csv`. Two tiers, sized by the build script against real
plate formats:

| Tier | Targets | Format | Spare wells |
|---|---|---|---|
| Core clinical | 46 | 48-plex qPCR array | 2 |
| Extended research | 83 | 96-plex array | 13 |

Pathway coverage:

| Code | Pathway | Genes |
|---|---|---|
| BER | base excision repair | 20 |
| HR | homologous recombination | 11 |
| SIG | damage signalling | 9 |
| NER | nucleotide excision repair | 8 |
| MMR | **mismatch repair** | 7 |
| NHEJ | non-homologous end joining | 7 |
| ICL | interstrand crosslink repair | 4 |
| DR | direct damage reversal | 3 |
| TLS | translesion synthesis | 3 |
| PROOF | **polymerase proofreading** | 3 |
| TEL | telomere maintenance | 2 |

Base excision repair carries the most genes deliberately: it is the pathway that
handles the panel's anchor lesion, and it has enough sequential steps that knowing
*which* step is limiting matters. A glycosylase excises the base; APE1 incises the
resulting abasic site; POLB fills the gap; XRCC1 scaffolds; LIG3 seals. **High
glycosylase activity with low APE1 is worse than low glycosylase activity**, because it
converts a modified base into an abasic site, which is more mutagenic than the lesion
it replaced. A pathway-level score would miss that; a step-resolved panel sees it.

Translesion synthesis is included for the same reason with opposite sign. REV3L and
REV1 provide damage **tolerance**, not repair — they let replication proceed past a
lesion at the cost of fidelity. **High TLS with low repair is the actively mutagenic
configuration**, and it looks like coping if you only measure whether replication
completed.

Normalisation uses the geometric mean of stability-selected housekeepers from the
six-gene candidate set, not a single reference gene. B2M is included as a candidate but
flagged: it is interferon-responsive, so inflammation moves it, and inflammation is
correlated with oxidative stress — exactly the confounding structure that would bias
this panel.

### 4.3 The damage-to-repair coupling, which is the point

Every damage analyte is mapped to the enzymes that handle it. The build script verifies
this mapping and reports **14 of 14 damage analytes coupled**:

| Damage analyte | Pathway | Enzymatic partners in panel |
|---|---|---|
| 8-oxodG | BER | OGG1, MUTYH |
| 8-oxoGuo | BER | NUDT1 (MTH1) |
| 8-oxoGua | BER | NUDT1 (MTH1) |
| FapyGua | BER | NEIL1, NEIL2 |
| 5-hmdC | BER | SMUG1, TDG |
| Tg | BER | NTHL1 |
| O6-medG | DR | **MGMT** |
| N7-meGua | BER | MPG |
| M1dG | NER | XPC, ERCC1, ERCC4 |
| εdA | NER | XPC, ERCC1, ERCC4 |
| **cdA** | **NER** | XPA, XPC, ERCC1, ERCC2, ERCC4, ERCC5, DDB2 |
| **cdG** | **NER** | XPA, XPC, ERCC1, ERCC2, ERCC4, ERCC5, DDB2 |
| 5-Cl-dC | BER | SMUG1, TDG |
| dU | BER | UNG, SMUG1, DUT, TYMS |

The RNA analyte is the interesting case. **Oxidised RNA is not repaired in place** —
there is no RNA equivalent of base excision repair. Its only enzymatic defence is
MTH1 hydrolysing 8-oxo-GTP out of the nucleotide pool before incorporation. So 8-oxoGuo
couples to NUDT1 and to nothing else, and that is a statement about biology rather than
a gap in the panel.

---

## 5. Quantification

### 5.1 Per-analyte normalisation

- **Spot urine**: divide by creatinine. Non-negotiable — unnormalised spot-urine
  concentrations are not comparable between subjects. Report as nmol/mol creatinine.
- **24-hour urine**: report as absolute excretion per 24 h, and as per kg body weight.
- **DNA hydrolysate**: report as lesions per 10⁶ dG, using the measured dG channel as
  denominator.
- **Flow**: report as median fluorescence intensity ratio to unstimulated control, and
  as percentage of cells above a control-defined threshold, both cell-cycle-gated.

### 5.2 Robust z-scoring on the log scale

Adduct distributions are right-skewed and outlier-prone. The panel scores on the log
scale against a laboratory-specific reference **median and interquartile range**, not
mean and standard deviation:

```
z_i = ( ln(x_i) - ln(median_ref,i) ) / ( IQR_ref,i / 1.349 )
```

The 1.349 divisor converts an interquartile range into a standard-deviation equivalent
for a normal distribution, so the score remains interpretable on the familiar z scale
while resisting the outliers that wreck a parametric z.

### 5.3 Precision weighting

Analytes are not measured equally well. An analyte at 25% coefficient of variation
should not carry the same weight in a composite as one at 5%. Weights are inverse to
log-scale measurement variance:

```
w_i = ( 1 / ln(1+CV_i)^2 ) / sum_j ( 1 / ln(1+CV_j)^2 )
```

In the worked example this gives 8-oxodG (6% CV) a weight of 0.392 and M1dG (15% CV) a
weight of 0.068 — a factor of six, reflecting how much each measurement can be trusted.

### 5.4 The three composite indices

- **ODI, Oxidative Damage Index** — precision-weighted mean of Module A lesion z-scores.
- **RCI, Repair Capacity Index** — reported *per pathway* (RCI_MMR, RCI_BER, RCI_NER,
  RCI_HR, RCI_NHEJ, RCI_DR) and as a mean across pathways.
- **SRI, Stress-to-Repair Index** — the headline number.

### 5.5 Why the headline index is a difference and not a ratio

The intuitive construction is damage divided by repair. It is wrong, and the reason is
worth stating because it is a trap that has appeared in published composite scores.

A z-scored denominator **crosses zero by construction** — that is what z-scoring means.
So ODI/RCI is undefined at the population median and explodes in a neighbourhood around
it, giving enormous scores to subjects who are close to average. The panel instead uses

```
SRI = z(ODI) - z(RCI)
```

On the log scale a difference is exactly what a ratio is reaching for
(`ln(a/b) = ln a - ln b`), and it is finite everywhere. Interpretation:

| SRI | State |
|---|---|
| > +0.5 | **decompensated** — damage exceeds repair capacity; the mutagenic risk state |
| −0.5 to +0.5 | **compensated** — damage matched by repair |
| < −0.5 | **repair reserve** |

### 5.6 What makes the result actionable

Run `python3 src/dna_stress_panel.py` for the worked arithmetic on illustrative values.
In that example SRI = +1.51, decompensated. The composite is driven by RCI_BER (−1.45)
and RCI_DR (−0.90), and the coupling table names OGG1, MUTYH, NUDT1 and MGMT as the
specific enzymes to follow up. **A single scalar damage score would have flagged the
subject without localising the defect.** Pathway resolution is the difference between a
result and a next step.

*The numbers in that example are invented to demonstrate the arithmetic. They are not
measurements and carry no reference-interval authority.*

---

## 6. Pre-analytics, where this panel is most likely to fail

Pre-analytical error dominates every other error source for these analytes. This
section is not boilerplate.

### 6.1 Urine

- First morning void, or timed 24-hour collection. Do not mix designs within a study.
- Aliquot and freeze at −80 °C within 2 hours.
- **Avoid repeated freeze-thaw.** Aliquot at first freeze into single-use volumes.
- Measure creatinine on the same aliquot as the adducts.
- Record BMI, smoking status and time of collection. All three move these analytes;
  urinary 8-isoprostane in particular shows a significant positive association with BMI.

### 6.2 Blood and DNA extraction — the critical path

- EDTA tubes; isolate PBMCs within 2 hours.
- **Extraction chemistry determines the 8-oxodG result.** Use a chelator-supplemented,
  phenol-free method (salting-out or anion-exchange). Add a metal chelator such as
  deferoxamine to suppress Fenton chemistry during lysis. Minimise air exposure and
  work cold.
- Enzymatic hydrolysis to nucleosides, not acid hydrolysis: acid conditions oxidise
  guanine and destroy the ring-opened lesions.
- Run the FPG/OGG1-modified comet assay on cells from the same draw as the artefact
  cross-check described in 1.2.

### 6.3 Tissue

FFPE is acceptable for MMR immunohistochemistry, MSI, MGMT methylation and HRD.

**FFPE is not acceptable for oxidative adduct quantification.** Formalin fixation and
prolonged block storage oxidise DNA. Adduct measurements require fresh or
snap-frozen tissue. This restriction is why Module A and Module C can run on different
matrices from the same subject and must not be forced onto one specimen type.

---

## 7. Validation plan before any clinical use

| Element | Target |
|---|---|
| Linearity | ≥ 5 calibration levels, r² ≥ 0.99, across the expected physiological range |
| Lower limit of quantification | CV ≤ 20% and accuracy 80–120% at the LLOQ |
| Intra-assay precision | CV ≤ 10% for Module A nucleosides; ≤ 15% for free bases |
| Inter-assay precision | CV ≤ 15% across ≥ 20 independent runs |
| Matrix effect | post-column infusion assessment for every analyte in every matrix |
| Internal standard M+3 bleed | quantify the analyte's natural M+3 contribution into the IS channel for the six analytes flagged at +3 nominal shift |
| Carryover | ≤ 20% of LLOQ in a blank following the highest calibrator |
| Reference interval | ≥ 120 healthy subjects per stratum, sex- and age-stratified, per CLSI-style design |
| Stability | freeze-thaw ×3, bench-top 24 h, long-term −80 °C at 6 and 12 months |
| Orthogonal confirmation | Module A 8-oxodG against FPG-modified comet on paired samples |
| Laboratory banding | establish local percentile distributions for micronucleus and comet endpoints rather than adopting published absolute cut-offs |

---

## 8. What the build checks found

Two findings from the automated verification changed the design rather than being
written up after the fact.

**The internal-standard threshold was initially wrong.** The first implementation
required an exact mass shift of ≥ 3.0 Da and failed four internal standards whose
shifts were 2.991–2.997 Da. Those are `[13C1,15N2]` and `[15N3]` labels — **+3 nominal**
mass shifts, fully resolvable on a unit-resolution quadrupole. The check was testing
exact mass where the instrument resolves nominal mass. It now computes the integer
nominal shift, requires ≥ 3, and issues an advisory (`PASS_VERIFY_M3`) at exactly +3
so that M+3 isotopic bleed-through from the analyte into the standard's channel is
verified during method development. Six analytes carry that advisory; it appears in the
validation plan in Section 7.

**Two internal standards were upgraded.** FapyGua and M1dG are both guanine-derived, so
`[15N5]`-labelled precursors are commercially available. Both were moved from +3 to
+5 nominal shift, removing the bleed-through concern entirely for those two channels.

**One spectral adjacency was flagged and left visible.** 53BP1 (PE) and pATM
(PE-Dazzle594) sit 32 nm apart on the same laser. Rather than quietly reassigning it,
the build reports it as spillover-critical so that compensation controls are treated as
mandatory for that pair.

**A validation rule refused to apply itself, correctly.** The glycosidic neutral-loss
check computes a product ion by subtracting the intact sugar from the precursor. For the
two cyclopurines it declines to do so and reports the reason: the C5'-C8 covalent bond
that defines the lesion is precisely what prevents the sugar leaving as a neutral
fragment. The validator encoding that distinction is what stops a plausible-looking but
chemically impossible transition from being accepted for cdA and cdG.

**One coupling gap turned out to be biology.** The first coupling run reported 8-oxoGuo
as uncoupled. That was correct as written but incomplete: MTH1 sanitises 8-oxo-GTP as
well as 8-oxo-dGTP, so the RNA analyte does have an enzymatic partner. The annotation
was corrected, and the underlying point — oxidised RNA is not repaired in place — is now
stated explicitly in Section 4.3 instead of appearing as a missing entry.

---

## 9. Honest limitations

- **No composite index in this panel is clinically validated.** ODI, RCI and SRI are
  defined here as a coherent scoring framework, not as validated clinical instruments.
  They require prospective outcome validation before any clinical decision rests on
  them. The individual components differ enormously in maturity: MMR immunohistochemistry
  and MSI are routine diagnostics, the micronucleus assay has prospective cancer-incidence
  data, and most Module A adducts beyond 8-oxodG are research-grade.
- **Reference intervals are laboratory-specific.** Published absolute values for these
  analytes are not transferable between laboratories, which is the central lesson of
  the standardisation literature cited throughout.
- **Transcript abundance is not enzyme activity.** Module C's expression panel measures
  message, not function. That is why the design pairs transcripts with functional
  readouts wherever one exists: PARP1 with PAR polymer, RAD51 transcript with RAD51
  foci, CHEK2 with pCHK2, H2AX with γH2AX, MGMT with promoter methylation. Where no
  functional counterpart exists, the transcript result is weaker evidence and should be
  reported as such.
- **The ex-vivo challenge assay measures lymphocyte repair capacity**, which is a proxy
  for the tissue of interest, not a measurement of it.
- **Mutational signature attribution needs mutation count.** Exome-scale data give
  unstable signature fits; the attribution layer in Section 4.1 needs genome-scale
  sequencing to be reliable.
- **Module A2 spans wildly different maturity levels.** Telomere length and
  mitochondrial copy number have large literatures and poor standardisation; PIG-A has a
  well-characterised background frequency; error-corrected sequencing is still being
  established for human biomonitoring. The nitrative and halogenative markers are
  research-grade.
- **Cell-free DNA is context, not a damage measure.** It reports that cells died, not
  that their genomes were damaged, and leukocyte lysis during sample handling contaminates
  the signal. It is included with that caveat attached rather than scored as damage.
- **Cost and throughput are not addressed here.** A three-module panel requiring
  LC-MS/MS, flow cytometry, immunohistochemistry, fragment analysis and an expression
  array is a research instrument. Deciding which modules survive into a deployable
  assay is a separate exercise from designing the complete one.

---

## 10. Files

```
data/dna_damage_adducts.csv        18 Module A analytes: formula, transition, internal standard
data/dna_stress_other_markers.csv  11 Module A2 markers on other platforms
data/ddr_flow_panel.csv            10 Module B flow parameters with laser/emission assignment
data/ddr_cytogenetic_assays.csv     5 cytogenetic and comet endpoints
data/repair_clinical_tests.csv      9 clinically established repair tests with pitfalls
data/repair_expression_panel.csv   83 repair genes across 11 pathways, 2 tiers
src/dna_stress_panel.py            builder, validator and scoring framework
src/build_report_pdf.py            renders the report to PDF
results/                           generated verification output
```

## Sources

Literature located through PubMed. Key references, with DOI links:

- Urinary 8-oxodG measurement consensus and the ELISA-versus-chromatography
  discrepancy — Cooke *et al.*, *FASEB J*
  ([DOI](https://doi.org/10.1096/fj.09-147124))
- Micronucleus frequency as a prospectively validated cancer-risk biomarker, 6718
  subjects — Bonassi *et al.*, *Carcinogenesis* 2006
  ([DOI](https://doi.org/10.1093/carcin/bgl177))
- Genetic damage biomarkers for cancer epidemiology — Fenech, *Toxicology* 2002
  ([DOI](https://doi.org/10.1016/s0300-483x(02)00480-8))
- Comet assay in human biomonitoring, international database comparison and baseline
  frequencies — *Mutat Res Rev Mutat Res*
  ([DOI](https://doi.org/10.1016/j.mrrev.2021.108371)) and
  ([DOI](https://doi.org/10.1016/j.mrrev.2019.108288))
- Urinary nucleic-acid oxidation markers and long-term mortality in type 2 diabetes —
  *Diabetes Care* ([DOI](https://doi.org/10.2337/dc11-1620)); RNA oxidation follow-up
  ([DOI](https://doi.org/10.1016/j.jdiacomp.2018.12.004))
- γH2AX as a quantitative DSB biomarker — *Aging*
  ([DOI](https://doi.org/10.18632/aging.100284)); γH2AX and other histone modifications
  in the clinic — *Biochim Biophys Acta*
  ([DOI](https://doi.org/10.1016/j.bbagrm.2012.02.021))
- γH2AX, apoptosis and radiosensitivity in patient lymphocytes — *Int J Radiat Biol*
  ([DOI](https://doi.org/10.1080/09553002.2025.2548465)); DSB repair capacity and
  normal-tissue toxicity — *J Radiat Res*
  ([DOI](https://doi.org/10.1093/jrr/rrae081))
- OGG1-driven base excision repair in disease — *Curr Opin Allergy Clin Immunol*
  ([DOI](https://doi.org/10.1097/ACI.0000000000000135))
- MMR protein expression and checkpoint-inhibitor treated disease — *J Cancer Res Clin
  Oncol* ([DOI](https://doi.org/10.1007/s00432-022-04002-4))
- Minimal microsatellite shift as an MSI interpretation pitfall — *Mod Pathol*
  ([DOI](https://doi.org/10.1038/s41379-018-0179-3))
- Homologous recombination deficiency biomarkers and assay limitations — *J Pers Med*
  ([DOI](https://doi.org/10.3390/jpm11070612))
- Temozolomide versus radiotherapy in elderly glioblastoma, the MGMT-predictive setting
  — *Lancet Oncol* ([DOI](https://doi.org/10.1016/S1470-2045(12)70265-6)); glioblastoma
  overview ([DOI](https://doi.org/10.1016/B978-0-12-802997-8.00023-2))
- Urinary 8-isoprostane reference values and BMI association, systematic review and
  meta-analysis — *Toxicol Lett*
  ([DOI](https://doi.org/10.1016/j.toxlet.2020.04.006))

Markers added in the expansion:

- 8-nitroguanine in inflammation-related carcinogenesis — Hiraku, *Environ Health Prev
  Med* ([DOI](https://doi.org/10.1007/s12199-009-0118-5))
- 8,5'-cyclopurine-2'-deoxynucleosides: formation, measurement, NER-dependent repair and
  biological effects — Jaruga & Dizdaroglu, *DNA Repair*
  ([DOI](https://doi.org/10.1016/j.dnarep.2008.06.005))
- Polymerase ε and δ proofreading mutations defining a hypermutated, microsatellite-stable
  cancer class — Briggs & Tomlinson, *J Pathol*
  ([DOI](https://doi.org/10.1002/path.4185)); germline PMS2 with somatic POLE exonuclease
  mutation and the resulting mutation rates — *J Pathol*
  ([DOI](https://doi.org/10.1002/path.4957))
- Mitochondrial common deletion elevated in blood with oxidative stress — *Mitochondrion*
  ([DOI](https://doi.org/10.1016/j.mito.2015.12.001))
- Telomere shortening with life stress, oxidative stress and telomerase activity — Epel
  *et al.*, *PNAS* ([DOI](https://doi.org/10.1073/pnas.0407162101))
- PIG-A gene mutation assay in human biomonitoring, including background mutant
  frequencies — *Environ Mol Mutagen* ([DOI](https://doi.org/10.1002/em.22577))
- Duplex (error-corrected) sequencing for in-vivo mutation measurement — *Biol Reprod*
  ([DOI](https://doi.org/10.1093/biolre/ioaf029))
- Direct detection and quantification of abasic sites, including methoxyamine binding —
  *Nucl Med Biol* ([DOI](https://doi.org/10.1016/j.nucmedbio.2009.07.007))

Mass-spectrometric transitions in Section 2 are computed from molecular formula by
`src/dna_stress_panel.py`. Nominal transitions for free bases and for 8-isoprostane are
literature-typical values marked `lit.` in the build output and must be confirmed on the
local instrument during method development.
