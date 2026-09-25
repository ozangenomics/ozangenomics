# PROTAC architecture and the SAR of *Primula* root constituents

**A computational-chemistry assessment: can anything in primrose root serve as a degrader module?**

---

## 0. Executive answer

**No constituent of *Primula* root is a PROTAC or a PROTAC-like molecule.** Across all 19
constituents and class models evaluated, the maximum ECFP4 Tanimoto similarity to any
PROTAC reference compound is **0.133**. The conventional threshold at which ECFP4
similarity begins to imply shared bioactivity is ~0.4. There is essentially zero
structural overlap.

That is the uninteresting half of the answer. The interesting half is that when the
extract is decomposed into its structural classes and each is scored against the
*property* requirements of the three PROTAC modules, the classes separate cleanly, and
one of them lands somewhere genuinely useful:

| *Primula* class | Best-fit PROTAC role | Verdict |
|---|---|---|
| Triterpenoid sapogenin (protoprimulagenin A) | rigid linker / scaffold | **Plausible** — conformationally locked, Fsp3 = 1.00, 0 rotatable bonds, 2 orthogonal exit vectors |
| Primeverose disaccharide | cleavable hydrophilic linker | **Plausible** — PEG-like polarity at half the rotatable-bond cost, glycosidase-cleavable |
| Phenolic glycosides (primulaverin, primverin, gaultherin, androsin) | none | Prodrug architecture, not a degrader module — but a *masking* lesson |
| Intact saponins (primulasaponins, priverosaponin B 22-acetate) | none | Far beyond any permeability envelope; **but** a credible delivery-adjuvant class |
| Flavonoids | none | Aerial-part comparators; no degrader relevance |

The single most defensible translational proposal from this analysis is **not** that a
*Primula* constituent becomes a warhead, linker, or E3 ligand. It is that
*Primula*-type saponins belong to the structural class being developed as
**endosomal-escape enhancers**, which is the exact liability of the cyclic-peptide
degrader in the attached schematic. Section 6 argues, from the aglycone, why *Primula*
is nonetheless a **poorer** starting point for that purpose than *Saponaria* or
*Quillaja*.

---

## 1. PROTAC chemistry, in the terms the schematic uses

A PROTAC (PROteolysis-TArgeting Chimera) is an event-driven, catalytic bifunctional
molecule. It does not inhibit its target; it re-routes it.

**The three modules.**

1. **Warhead** — binds the protein of interest (POI). It need not be an inhibitor; it
   only needs affinity and a solvent-exposed vector. In the attached schematic the
   warhead is a **cyclic peptide**, which is the design choice you make when the target
   is a flat protein-protein interface with no small-molecule pocket. Cyclic peptides
   buy you interface-sized contact area at the cost of mass and polarity.
2. **Linker** — sets the geometry of the ternary complex. This is the module that
   actually decides whether the molecule works. Its length, rigidity, and attachment
   vectors determine whether the E3's catalytic lysine-presenting surface can reach a
   surface lysine on the POI.
3. **E3 ligand** — recruits an E3 ubiquitin ligase. In practice this is a small set:
   cereblon (CRBN, via the glutarimide of thalidomide/lenalidomide/pomalidomide), VHL
   (via a hydroxyproline peptidomimetic such as VH032), and a longer tail of MDM2, IAP,
   DCAF15/16, KEAP1.

**The mechanism that constrains the chemistry.**

- Degradation requires a productive **ternary complex** POI·PROTAC·E3, not merely two
  binary binding events. Positive cooperativity (α > 1) from new POI–E3 contacts is what
  separates good degraders from inert bifunctionals.
- The **hook effect**: at high concentration, binary complexes outcompete the ternary
  complex and degradation falls off. Dose-response is bell-shaped, not sigmoidal.
- The molecule is **catalytic** — one PROTAC can cycle through many POI copies — which
  is why sub-stoichiometric exposure can suffice, and why potency is measured as
  DC50/Dmax rather than IC50.

**The property space this forces you into.**

PROTACs are **beyond-Rule-of-5 (bRo5)**. The measured reference points in this analysis:

| Compound | MW | cLogP | TPSA | HBD | HBA | RotB | Fsp3 | QED |
|---|---|---|---|---|---|---|---|---|
| Pomalidomide (CRBN ligand) | 273.2 | −0.33 | 109.6 | 2 | 5 | 1 | 0.23 | 0.54 |
| VH032 (VHL ligand) | 472.6 | 2.25 | 111.6 | 3 | 6 | 6 | 0.50 | 0.60 |
| Vepdegestrant / ARV-471 (clinical ER degrader) | 723.9 | 6.05 | 96.4 | 2 | 7 | 7 | 0.40 | 0.22 |
| Bavdegalutamide / ARV-110 (clinical AR degrader) | 812.3 | 3.70 | 181.2 | 2 | 12 | 9 | 0.46 | 0.30 |
| Cyclosporin A (oral macrocycle benchmark) | 1202.6 | 3.27 | 278.8 | 5 | 12 | 15 | 0.79 | 0.15 |

Note what the two clinical degraders have in common and what they do *not*. MW 724–812
and 7–9 rotatable bonds, yes. But **HBD = 2 for both**, and TPSA 96–181. Hydrogen-bond
donor count is the variable that oral bRo5 chemistry defends most aggressively, because
desolvation of donors is what kills passive permeability. Cyclosporin A is orally
bioavailable at MW 1203 precisely because N-methylation and intramolecular hydrogen
bonding hold it to **HBD = 5**. Keep that number in mind for Section 5.

---

## 2. What is actually in the root, and a pharmacognostic warning

According to PubMed, the two most directly relevant analyses are:

- Trendafilova *et al.* (2026), *Plants*, on *Primula veris* subsp. *columnae*, which by
  UHPLC-MS/MS and NMR reports from the **underground parts**: a new triterpene saponin
  primulasaponin VI, together with primulasaponins I, III and V; **priverosaponin B
  22-acetate**; **primulaverin**; **primverin**; **gaultherin**; and
  **3-methoxy-4-primeverosylacetophenone** (androsin). Surface flavonoids and flavonoid
  glycosides were isolated from the **aerial** parts.
  [DOI](https://doi.org/10.3390/plants15152259)
- Müller, Ganzera and Stuppner (2005), *J. Chromatogr. A*, the validated LC-ELSD/MS
  method for *Primula elatior* and *Primula veris* root, reporting **total saponin content
  up to 14.9% in *P. veris* roots** and **primeverin as the dominant phenolic glycoside at
  0.64–1.42%**, with the two species distinguishable by saponin pattern.
  [DOI](https://doi.org/10.1016/j.chroma.2005.10.067)

**The warning.** The question asks about ***Primula vulgaris***. The European
Pharmacopoeia drug *Primulae radix* is defined as *Primula veris* L. and/or *Primula
elatior* (L.) Hill — **not** *P. vulgaris*. Essentially all quantitative constituent data,
including both references above, derive from *P. veris* / *P. elatior*. *P. vulgaris*
shares the chemotype family (13,28-epoxy-oleanane saponins plus primeveroside-linked
phenolics), but congener distribution and absolute saponin content are **not**
transferable between these species, and the 2005 study explicitly shows saponin pattern
is species-discriminating. Any *P. vulgaris* programme must re-run the constituent
analysis on authenticated *P. vulgaris* material rather than inheriting *P. veris* data.
This document therefore treats the constituent list as **genus-level** and flags it
throughout.

---

## 3. Method

Structures were assembled from the literature constituent list and validated by an
independent check: **every drawn structure's RDKit-computed molecular formula was
required to match the literature-reported formula.** All 13 pass
(`results/formula_validation.csv`). Reference PROTAC structures were taken from ChEMBL
(thalidomide CHEMBL468, lenalidomide CHEMBL848, pomalidomide CHEMBL43452, vepdegestrant
CHEMBL5095210, bavdegalutamide CHEMBL4862963, cyclosporine CHEMBL160).

Three deliberate epistemic choices, because they affect how much weight the numbers bear:

1. **Stereochemistry is asserted only where verifiable.** The phenolic glycosides and
   flavonoids carry full stereochemistry. The triterpenoid sapogenins are written
   **constitution-only**. Every descriptor used here (MW, TPSA, HBD, HBA, rotatable
   bonds, Fsp3, Crippen logP) is a constitution-level property, and default Morgan
   fingerprints do not use chirality, so no conclusion in this document depends on
   stereochemistry that was not independently confirmed.
2. **Intact saponins are modelled as a series, not as named congeners.** Machine-readable
   structures for individual primulasaponins are not reliably available. Instead
   `src/build_models.py` programmatically glycosylates protoprimulagenin A with 1–6 sugar
   residues, bracketing the property envelope of the whole class (MW 621 → 1340). This is
   legitimate because for isomeric glycans the descriptors above depend on glycan
   *composition*, not on linkage regiochemistry.
3. **The 4-OMe / 5-OMe assignment between primulaverin and primverin** is inconsistent
   across secondary sources. They are treated as a positional isomer pair; since they are
   isomers, every computed descriptor is identical and no conclusion turns on which is
   which.

Reproduce with `python3 src/run_analysis.py`.

---

## 4. Result 1 — structural similarity is nil, and that is informative

Maximum ECFP4 Tanimoto to any PROTAC reference, per constituent:

| Constituent | max Tc | nearest reference |
|---|---|---|
| Methyl salicylate | 0.133 | Thalidomide |
| Primulaverin | 0.133 | Bavdegalutamide |
| Gaultherin | 0.130 | Thalidomide |
| Primverin | 0.124 | Bavdegalutamide |
| Androsin | 0.120 | VH032 |
| Methyl 4-methoxysalicylate | 0.103 | Bavdegalutamide |
| Kaempferol | 0.102 | Vepdegestrant |
| Rutin | 0.099 | VH032 |
| Quercetin | 0.089 | Vepdegestrant |
| Priverogenin B 22-acetate (model) | 0.086 | VH032 |
| Primeverose | 0.079 | PEG3 linker |
| Protoprimulagenin A | 0.066 | piperazine-piperidine linker |
| Saponin models (1–6 sugars) | 0.059–0.066 | linker / cyclosporin A |

Nothing approaches 0.4. The practical reading: **there is no "natural CRBN binder" or
"natural VHL binder" hiding in primrose root**, and no virtual-screening shortcut where
this extract is concerned. The glutarimide of the CRBN ligands and the hydroxyproline
thiazole of VH032 are synthetic pharmacophores with no natural counterpart here.

Note also that the *nearest neighbour* assignments are chemically sensible even at low
absolute similarity — the saponins' nearest neighbours are the linkers and cyclosporin A
(aliphatic, sp3-rich, polar), the flavonoids' nearest neighbour is vepdegestrant (the
only phenol-bearing reference). The fingerprint is behaving; the compounds are simply far
apart.

---

## 5. Result 2 — the property-space analysis, module by module

### 5.1 The sapogenin as a rigid linker scaffold — the strongest structural fit

Protoprimulagenin A (13,28-epoxy-oleanane-3,16-diol), the shared aglycone:

| | protoprimulagenin A | PEG3 linker | piperazine-piperidine linker |
|---|---|---|---|
| MW | 458.7 | 150.2 | 183.3 |
| cLogP | 6.35 | −1.00 | −0.11 |
| TPSA | 49.7 | 58.9 | 27.3 |
| HBD | 2 | 2 | 2 |
| **Rotatable bonds** | **0** | 7 | 2 |
| **Fsp3** | **1.00** | 1.00 | 1.00 |

**Zero rotatable bonds, Fsp3 = 1.00, and two chemically differentiated hydroxyl exit
vectors** (C-3 secondary equatorial, C-16 secondary) on a pentacyclic cage. This is a
conformationally locked molecular ruler.

Why that matters for degraders: the field's linker evolution has run steadily *away* from
flexible PEG and alkyl chains and *toward* rigid, sp3-rich spacers — piperazines,
piperidines, spirocycles, bicyclo[1.1.1]pentanes. The reasons are (i) the entropic
penalty of freezing out a floppy linker on ternary complex formation, (ii) the need to
reproducibly present the warhead and E3 ligand at a fixed distance and angle, and (iii)
that rotatable-bond count and Fsp3 track with both permeability and developability in
bRo5 space. An oleanane cage delivers all three, plus differentiated exit vectors so the
two ends are synthetically orthogonal.

**The disqualifying liability is lipophilicity.** cLogP = 6.35 for the bare aglycone. A
degrader is already carrying 700–800 Da of warhead and E3 ligand; ARV-471 is at cLogP
6.05 for the *entire molecule*. Grafting a cLogP 6.35 spacer into that is not viable —
you would be building toward cLogP > 10, with the attendant metabolic instability,
plasma-protein binding, and aggregation. Any serious use of this scaffold requires
polarity to be engineered back in (retaining one or two sugars, or introducing the
C-23/C-16/C-28 oxidations that the *Saponaria* biosynthesis literature characterises),
and that trades away the very rigidity that made it attractive.

**Verdict: structurally elegant, chemically expensive.** It is a real idea, not a good
first idea.

### 5.2 Primeverose as a cleavable hydrophilic linker — the most practical fit

Primeverose (6-*O*-β-D-xylopyranosyl-D-glucose) against the standard hydrophilic linker:

| | primeverose | PEG3 linker |
|---|---|---|
| MW | 312.3 | 150.2 |
| cLogP | −4.76 | −1.00 |
| TPSA | 169.3 | 58.9 |
| HBD | 7 | 2 |
| **Rotatable bonds** | **3** | **7** |
| Fsp3 | 1.00 | 1.00 |

Two quantitative observations.

**First, sugars buy polarity at a lower conformational-entropy cost than PEG.** Primeverose
delivers 2.9× the TPSA of PEG3 with **fewer than half the rotatable bonds** (3 vs 7).
Pyranose rings carry their oxygens on a locked chair; PEG carries them on a freely
rotating chain. If the design goal is "add hydrophilicity without adding floppiness" —
which is precisely the bRo5 degrader problem — a pyranose is the more efficient unit.

**Second, the glycosidic bond is a cleavable trigger.** β-glycosidase-cleavable linkers
are established antibody-drug-conjugate chemistry (β-glucuronide linkers in particular),
and *Primula* root itself demonstrates the mechanism: endogenous β-glycosidase hydrolysis
of primulaverin and gaultherin releases the volatile aglycones. A sugar linker is
therefore not just a spacer but a **conditional release element** — one could mask a
degrader's polarity or target its release to tissue with the appropriate glycosidase
activity.

**Verdict: the most directly transferable idea in the extract**, and notably it is a
*motif* transfer, not a compound transfer.

### 5.3 The phenolic glycosides — a prodrug architecture, not a module

| | aglycone | glycoside | Δ |
|---|---|---|---|
| Methyl 4-methoxysalicylate → primulaverin | MW 182.2, cLogP +1.19, TPSA 55.8, HBD 1 | MW 476.4, cLogP −2.88, TPSA 193.8, HBD 6 | +294 Da, −4.07 logP |
| Methyl salicylate → gaultherin | MW 152.2, cLogP +1.18, TPSA 46.5, HBD 1 | MW 446.4, cLogP −2.88, TPSA 184.6, HBD 6 | +294 Da, −4.06 logP |
| Apocynin → androsin | MW 166.2, cLogP +1.60, TPSA 46.5, HBD 1 | MW 460.4, cLogP −2.46, TPSA 184.6, HBD 6 | +294 Da, −4.06 logP |

These are not degrader modules. At MW 446–476 they are in the right mass range for an
"E3 ligand plus linker" fragment, but at TPSA ~185–194 and HBD 6 they are roughly twice
as polar as either clinical degrader, and the aglycones are small phenolic fragments with
no reported E3 affinity.

What they *do* illustrate, with unusual clarity because the transformation is constant
across all three pairs (+294 Da, −4.07 logP, HBD 1 → 6), is a **uniform, enzymatically
reversible polarity switch**. That is the structural logic a bRo5 degrader prodrug would
need. The plant is running the masking strategy that medicinal chemistry applies to
overly lipophilic or overly permeable compounds — in reverse, and with a clean traceless
release.

### 5.4 The intact saponins — outside every permeability envelope

The modelled glycosylation series, against the oral macrocycle benchmark:

| | MW | cLogP | TPSA | **HBD** | HBA | RotB |
|---|---|---|---|---|---|---|
| Cyclosporin A (orally bioavailable) | 1202.6 | 3.27 | 278.8 | **5** | 12 | 15 |
| Saponin model, 3 sugars | 899.1 | 1.49 | 246.7 | **9** | 16 | 7 |
| Saponin model, 4 sugars | 1061.3 | −0.68 | 325.8 | **12** | 21 | 10 |
| Saponin model, 5 sugars | 1193.4 | −2.22 | 384.8 | **14** | 25 | 12 |
| Saponin model, 6 sugars | 1339.5 | −3.37 | 443.7 | **16** | 29 | 14 |

This is the decisive comparison in the whole analysis. At **essentially identical mass**,
cyclosporin A carries **5** hydrogen-bond donors and a *Primula*-type pentaglycoside
carries **14**, with TPSA 385 against 279 and cLogP −2.2 against +3.3.

Cyclosporin is the existence proof that a 1200 Da molecule can cross a membrane. It does
so by N-methylating its amides and folding to bury the rest, holding HBD to 5 and
maintaining positive logP. A saponin does the opposite of every one of those things. Its
polarity is displayed on a solvent-exposed glycan that cannot be conformationally hidden.

**Passive permeability is therefore not merely poor — it is structurally precluded.** An
intact *Primula* saponin cannot be a cyclic-peptide warhead substitute, cannot be carried
as a linker, and cannot be dragged across a membrane as part of any conjugate. This also
means the amphiphilicity that gives these compounds their pharmacological activity
(membrane interaction, hemolysis, the gastric vagal reflex behind the traditional
expectorant use) is a *surface* activity, not an intracellular one.

---

## 6. The one genuinely promising angle, and the structural reason it is second-best

The cyclic-peptide degrader in the attached schematic has a specific, well-known
liability: a macrocyclic peptide warhead pushes the conjugate to MW 1200–1800, where
passive diffusion fails and material that does enter the cell is **trapped in
endosomes**. Endosomal entrapment, not binding affinity, is the limiting step.

There is an established class of compounds for exactly this problem, and it is this
structural class. According to PubMed, Lacchini *et al.* (2025), *Plant Biotechnology
Journal*, describe plant-derived triterpenoid saponins with **endosomal escape-enhancing
(EEE)** properties that improve cytosolic delivery of protein- and DNA-based
therapeutics by disrupting endosomal membranes — naming **QS-21** (the *Quillaja* vaccine
adjuvant) and **SO1861** from *Saponaria officinalis*, and characterising the cytochrome
P450s performing the **C23, C28 and C16 oxidations** on the oleanane skeleton required to
build SO1861. [DOI](https://doi.org/10.1111/pbi.70122)

*Primula* saponins are the same architectural class: an oleanane-derived aglycone bearing
a branched oligosaccharide, amphiphilic, membrane-active. So the honest proposal is that
these compounds belong on the **delivery** side of a large-degrader programme — as
co-formulated permeation or endosomal-escape enhancers — and not in the degrader molecule
at all.

**But the SAR argues *Primula* is the wrong genus for it.** The EEE saponins are built on
aglycones bearing a **C-28 carboxylic acid**, which carries the second (ester-linked)
glycan chain that makes them *bisdesmosidic*. In protoprimulagenin A, **C-28 is consumed
by the 13,28-ether bridge that defines the entire *Primula* sapogenin class.** It is not
a carboxylic acid and cannot be acylglycosylated. The structural consequence:

- *Primula*-type saponins are constitutionally **barred from forming the C-28 ester glycan
  chain** that characterises QS-21 and SO1861.
- They are therefore **monodesmosidic** (glycan at C-3 only). Monodesmosidic saponins are
  the strongly hemolytic, membrane-lytic ones; converting them to bisdesmosides is the
  classical way to attenuate that. For *Primula*, that route is closed at the aglycone.

This predicts a **narrower therapeutic window** for *Primula* saponins as delivery
enhancers than for *Saponaria*/*Quillaja* saponins, for a reason visible directly in the
aglycone connectivity. It is a falsifiable structural prediction, and it is the most
useful thing this analysis produces.

*Confidence:* the C-28 ether assignment in the *Primula* sapogenin class and the
mono-/bisdesmosidic hemolysis SAR are both well established; the specific inference
joining them is mine and has not been experimentally tested.

---

## 7. What I would actually do next

1. **Authenticate the species.** Re-run the Müller LC-ELSD/MS method on voucher-verified
   *P. vulgaris* root. Do not carry *P. veris* composition into a *P. vulgaris* programme.
   Everything downstream depends on this.
2. **Test the falsifiable prediction directly.** Hemolysis assay plus an endosomal-escape
   readout (e.g. saporin potentiation) on purified *Primula* monodesmosidic saponin versus
   SO1861 as the bisdesmosidic control. If the therapeutic-window gap predicted in
   Section 6 does not appear, the C-28 argument is wrong and should be discarded.
3. **Prototype the primeverose linker.** Synthesise a pomalidomide–primeveroside–warhead
   conjugate and compare it head-to-head against the matched PEG3 analogue on ternary
   complex formation, DC50/Dmax, and permeability. The hypothesis under test is the
   Section 5.2 claim: equal or better solubility at lower conformational entropy cost.
4. **Evaluate the sapogenin as a rigid spacer** *in silico* first. Enumerate C-3 and C-16
   exit-vector geometries, measure the accessible distance and angle range against the
   linker geometries in solved CRBN and VHL ternary-complex structures, and only commit
   to synthesis if the oleanane cage spans a useful region of that space. Budget for
   polarity rescue from the outset.
5. **Do not pursue** a *Primula* constituent as a warhead or E3 ligand. Section 4 closes
   that question.

---

## 8. Limitations

- Constituent identities are genus-level (*P. veris* / *P. elatior* / *P. veris* subsp.
  *columnae*); *P. vulgaris*-specific composition was not available and is **not** assumed
  to be identical. See Section 2.
- Intact saponins are class models spanning 1–6 sugars, not named congeners. They bracket
  the property envelope; they are not structures of specific isolates.
- Sapogenin stereochemistry is unassigned by design. No conclusion depends on it.
- cLogP is Crippen (atom-contribution). For amphiphiles of this size, calculated logP is
  unreliable in absolute terms; the comparisons here are used ordinally.
- Descriptor-space arguments bound what is *plausible*. They do not establish binding,
  degradation, or permeability, all of which require measurement.
- PubChem was unreachable under this environment's egress policy, so structures were
  drawn from literature descriptions and validated by molecular-formula agreement rather
  than retrieved from a structure database.

---

## References

Retrieved via PubMed:

1. Trendafilova A, Ivanova V, Novakovic M, Jadranin M, Staleva P, Dimitrova P, Veleva R,
   Topouzova-Hristova T, Paunova-Krasteva T. *Comparative Phytochemical Profiling and In
   Vitro Antibiofilm Activity of the Aerial and Underground Parts of Primula veris subsp.
   columnae (Ten.) Lüdi.* Plants (Basel). 2026;15(15):2259.
   [DOI](https://doi.org/10.3390/plants15152259)
2. Müller A, Ganzera M, Stuppner H. *Analysis of phenolic glycosides and saponins in
   Primula elatior and Primula veris (primula root) by liquid chromatography, evaporative
   light scattering detection and mass spectrometry.* J Chromatogr A. 2005;1112(1-2):218-23.
   [DOI](https://doi.org/10.1016/j.chroma.2005.10.067)
3. Lacchini E, Qu T, Moses T, Volkov AN, Goossens A. *Engineering Gypsophila elegans hairy
   root cultures to produce endosomal escape-enhancing saponins.* Plant Biotechnol J.
   2025;23(8):3068-3082. [DOI](https://doi.org/10.1111/pbi.70122)

Reference structures: ChEMBL v34 (CHEMBL468, CHEMBL848, CHEMBL43452, CHEMBL5095210,
CHEMBL4862963, CHEMBL160).
