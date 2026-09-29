"""
Multiplex panel builder and validator for DNA stress (damage) and DNA repair capacity.

Design intent
-------------
8-oxo-dG ("8-OHdG") cannot be multiplexed with a repair protein such as MSH2: one is a
small molecule measured by mass spectrometry, the other an epitope measured by antibody.
The panel is therefore modular by analyte class and matrix, and each module is internally
multiplexable in a single acquisition:

    Module A  oxidative / alkylation adducts   LC-MS/MS, one MRM method     urine + DNA hydrolysate
    Module B  cellular damage response        one stained tube, flow       PBMC
    Module C  repair capacity                 IHC + MSI/HRD + expression   tissue / PBMC

Verification built into the run (no external dependencies):
  1. every Module A precursor m/z is recomputed from the molecular formula
  2. every 2'-deoxynucleoside / ribonucleoside product ion is recomputed from the
     glycosidic neutral loss and checked against the declared transition
  3. no two analytes share a precursor->product pair that chromatography must resolve
  4. each analyte's stable-isotope internal standard is checked for a resolvable mass
     shift and for an intact elemental backbone
  5. flow panel emission peaks are checked for spectral crowding
  6. each core damage analyte is checked for an enzymatic partner in the repair module
  7. panel tiers are checked against real plate formats

Usage:  python3 src/dna_stress_panel.py
"""

import csv
import math
import os
import re
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")

# ----------------------------------------------------------------------------------
# Monoisotopic masses (CODATA / AME2020 rounded to 7 dp)
# ----------------------------------------------------------------------------------
MASS = {
    "H": 1.0078250319, "D": 2.0141017780,
    "C": 12.0000000000, "13C": 13.0033548378,
    "N": 14.0030740052, "15N": 15.0001088984,
    "O": 15.9949146221, "18O": 17.9991604,
    "P": 30.97376151,   "S": 31.97207069,
}
PROTON = 1.00727646688

# glycosidic neutral losses giving the protonated base as product ion
LOSS_DEOXYRIBOSE = 5 * MASS["C"] + 8 * MASS["H"] + 3 * MASS["O"]   # C5H8O3  116.0473
LOSS_RIBOSE      = LOSS_DEOXYRIBOSE + MASS["O"]                     # C5H8O4  132.0423

TOKEN = re.compile(r"(?:\[(\d+[A-Z][a-z]?)\]|([A-Z][a-z]?))(\d*)")


def parse_formula(formula):
    """Parse a formula with optional bracketed isotopes, e.g. '[13C]1C9H16[15N]2O7'."""
    counts = Counter()
    pos = 0
    for m in TOKEN.finditer(formula):
        if m.start() != pos:
            raise ValueError("unparsed text in %r at %d" % (formula, pos))
        pos = m.end()
        sym = m.group(1) or m.group(2)
        n = int(m.group(3)) if m.group(3) else 1
        counts[sym] += n
    if pos != len(formula):
        raise ValueError("trailing text in %r" % formula)
    return counts


def mono_mass(formula):
    counts = parse_formula(formula)
    total = 0.0
    for sym, n in counts.items():
        if sym not in MASS:
            raise ValueError("unknown isotope %r in %r" % (sym, formula))
        total += MASS[sym] * n
    return total


def heavy_atom_skeleton(formula):
    """Element composition ignoring isotopic labelling and ignoring hydrogen.

    Two compounds sharing this skeleton are the same molecule up to labelling and
    exchangeable hydrogen, which is the requirement for a valid internal standard.
    """
    skel = Counter()
    for sym, n in parse_formula(formula).items():
        el = re.sub(r"^\d+", "", sym)
        if el in ("H", "D"):
            continue
        skel[el] += n
    return skel


def mz(formula, mode, charge=1):
    m = mono_mass(formula)
    return (m + PROTON) if mode == "positive" else (m - PROTON)


# ----------------------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------------------
def load(name):
    with open(os.path.join(DATA, name), newline="") as fh:
        return list(csv.DictReader(fh))


# ----------------------------------------------------------------------------------
# Check 1-4 : Module A mass-spectrometric integrity
# ----------------------------------------------------------------------------------
def validate_module_a(adducts, tol_ppm=5.0, tol_nominal=0.15):
    """Recompute every transition from first principles and report per-analyte status."""
    out = []
    for a in adducts:
        row = {
            "abbrev": a["abbrev"],
            "analyte": a["analyte"],
            "formula": a["formula"],
            "species_type": a["species_type"],
            "mode": a["precursor_mode"],
        }
        prec = mz(a["formula"], a["precursor_mode"])
        row["precursor_mz_calc"] = round(prec, 4)

        # product ion: computable only for nucleosides, where the sugar is lost intact
        st = a["species_type"]
        if st == "2'-deoxynucleoside":
            calc = prec - LOSS_DEOXYRIBOSE
            row["product_basis"] = "loss of C5H8O3 (116.0473)"
        elif st == "ribonucleoside":
            calc = prec - LOSS_RIBOSE
            row["product_basis"] = "loss of C5H8O4 (132.0423)"
        else:
            calc = None
            row["product_basis"] = "not glycosidic - literature transition, confirm on instrument"

        declared = float(a["product_nominal"])
        row["product_mz_declared"] = declared
        if calc is None:
            row["product_mz_calc"] = ""
            row["product_check"] = "NOT_APPLICABLE"
            row["product_delta_mDa"] = ""
        else:
            row["product_mz_calc"] = round(calc, 4)
            delta = abs(calc - declared)
            row["product_delta_mDa"] = round(delta * 1000, 1)
            row["product_check"] = "PASS" if delta <= tol_nominal else "FAIL"

        # internal standard: resolvable shift and intact skeleton.
        #
        # The threshold is applied to the NOMINAL (integer) mass shift, not the exact
        # shift. A unit-resolution quadrupole separates nominal masses, so a [13C1,15N2]
        # label whose exact shift is 2.9974 Da is a +3 channel and is perfectly
        # resolvable; testing the exact shift against 3.0 would reject it spuriously.
        # A +3 label is the practical minimum because the analyte's own natural M+3
        # isotopologue bleeds into the standard's channel, so +3 is passed with an
        # advisory to verify that bleed-through during method development.
        is_f = a["is_formula"]
        is_mass = mono_mass(is_f)
        an_mass = mono_mass(a["formula"])
        shift = is_mass - an_mass
        nominal_shift = int(round(is_mass) - round(an_mass))
        row["internal_standard"] = a["internal_standard"]
        row["is_precursor_mz_calc"] = round(mz(is_f, a["precursor_mode"]), 4)
        row["is_mass_shift_Da"] = round(shift, 4)
        row["is_nominal_shift"] = nominal_shift
        same_skeleton = heavy_atom_skeleton(is_f) == heavy_atom_skeleton(a["formula"])
        if not same_skeleton:
            row["is_check"] = "FAIL_SKELETON"
        elif nominal_shift < 3:
            row["is_check"] = "FAIL_SHIFT_TOO_SMALL"
        elif nominal_shift == 3:
            row["is_check"] = "PASS_VERIFY_M3"
        else:
            row["is_check"] = "PASS"
        row["matrix"] = a["matrix"]
        row["role_in_panel"] = a["role_in_panel"]
        out.append(row)
    return out


def check_transition_collisions(validated, tol=0.2):
    """Two analytes sharing a precursor->product pair must be resolved by retention time."""
    collisions = []
    for i in range(len(validated)):
        for j in range(i + 1, len(validated)):
            a, b = validated[i], validated[j]
            if a["mode"] != b["mode"]:
                continue
            dp = abs(a["precursor_mz_calc"] - b["precursor_mz_calc"])
            dq = abs(a["product_mz_declared"] - b["product_mz_declared"])
            if dp <= tol and dq <= tol:
                collisions.append((a["abbrev"], b["abbrev"], round(dp, 3), round(dq, 3)))
    return collisions


def check_isobaric_precursors(validated, tol=0.2):
    """Same precursor but different product ion: acceptable, but flag it as comp-critical."""
    flags = []
    for i in range(len(validated)):
        for j in range(i + 1, len(validated)):
            a, b = validated[i], validated[j]
            if a["mode"] != b["mode"]:
                continue
            dp = abs(a["precursor_mz_calc"] - b["precursor_mz_calc"])
            dq = abs(a["product_mz_declared"] - b["product_mz_declared"])
            if dp <= tol and dq > tol:
                flags.append((a["abbrev"], b["abbrev"], "shared precursor, distinct product"))
            if dp > tol and dq <= tol:
                flags.append((a["abbrev"], b["abbrev"], "shared product, distinct precursor"))
    return flags


# ----------------------------------------------------------------------------------
# Check 5 : flow cytometry spectral crowding
# ----------------------------------------------------------------------------------
def validate_flow_panel(flow, min_gap_nm=40):
    per_laser = defaultdict(list)
    for f in flow:
        per_laser[f["laser_nm"]].append((int(f["emission_peak_nm"]), f["target"], f["fluorochrome"]))
    report = []
    for laser in sorted(per_laser, key=int):
        chans = sorted(per_laser[laser])
        for k in range(len(chans) - 1):
            gap = chans[k + 1][0] - chans[k][0]
            report.append({
                "laser_nm": laser,
                "channel_low": "%s (%s, %dnm)" % (chans[k][1], chans[k][2], chans[k][0]),
                "channel_high": "%s (%s, %dnm)" % (chans[k + 1][1], chans[k + 1][2], chans[k + 1][0]),
                "gap_nm": gap,
                "status": "OK" if gap >= min_gap_nm else "SPILLOVER_CRITICAL",
            })
    return report, {l: len(v) for l, v in per_laser.items()}


# ----------------------------------------------------------------------------------
# Check 6 : damage analyte -> repair enzyme coupling
# ----------------------------------------------------------------------------------
def build_damage_repair_pairs(adducts, genes):
    partner_index = defaultdict(list)
    for g in genes:
        for token in [t.strip() for t in g["damage_analyte_partner"].split("/")]:
            if token and token.lower() != "none":
                partner_index[token].append(g)
    pairs = []
    for a in adducts:
        if a["lesion_class"] in ("normaliser", "methylated base (not damage)",
                                 "unmodified nucleoside", "lipid peroxidation (not DNA)"):
            continue
        hits = partner_index.get(a["abbrev"], [])
        pairs.append({
            "damage_analyte": a["abbrev"],
            "lesion_class": a["lesion_class"],
            "n_repair_partners": len(hits),
            "repair_genes": "; ".join(g["gene"] for g in hits) or "(none in panel)",
            "pathways": "; ".join(sorted({g["pathway_code"] for g in hits})) or "-",
            "status": "COUPLED" if hits else "UNCOUPLED",
        })
    return pairs


# ----------------------------------------------------------------------------------
# Check 7 : tier sizing against plate formats
# ----------------------------------------------------------------------------------
PLATE_FORMATS = [(48, "48-plex qPCR array"), (96, "96-plex qPCR / nCounter panel"),
                 (384, "384-well high-throughput array")]


def size_tiers(genes):
    core = [g for g in genes if g["tier"] == "core"]
    ext = [g for g in genes if g["tier"] == "extended"]
    ref = [g for g in genes if g["tier"] == "normaliser"]
    tiers = [
        ("Core clinical", len(core) + len(ref), core + ref),
        ("Extended research", len(core) + len(ext) + len(ref), core + ext + ref),
    ]
    rows = []
    for name, n, _ in tiers:
        fit = next(((cap, lab) for cap, lab in PLATE_FORMATS if n <= cap), (None, "exceeds 384"))
        rows.append({
            "tier": name, "n_targets": n, "fits_format": fit[1],
            "format_capacity": fit[0] if fit[0] else "",
            "spare_wells": (fit[0] - n) if fit[0] else "",
        })
    return rows, {"core": len(core), "extended": len(ext), "normaliser": len(ref)}


# ----------------------------------------------------------------------------------
# Composite scoring framework
# ----------------------------------------------------------------------------------
def robust_z(values, ref_median, ref_iqr):
    """Robust z on the log scale. Adduct distributions are right-skewed and
    outlier-prone, so median/IQR is used rather than mean/SD."""
    if ref_iqr <= 0:
        raise ValueError("reference IQR must be positive")
    scale = ref_iqr / 1.349  # IQR -> SD equivalent for a normal distribution
    return [(math.log(v) - math.log(ref_median)) / scale for v in values]


def precision_weights(cvs):
    """Inverse-variance weights on the log scale: an analyte measured with a 25%
    CV should not carry the same weight in a composite as one measured at 5%."""
    raw = [1.0 / (math.log1p(cv) ** 2) for cv in cvs]
    tot = sum(raw)
    return [r / tot for r in raw]


def composite(zs, weights):
    return sum(z * w for z, w in zip(zs, weights))


def stress_repair_index(odi_z, rci_z):
    """Damage-to-repair balance.

    Defined as a DIFFERENCE of z-scores, not a ratio. A ratio ODI/RCI is unstable
    because a z-scored denominator crosses zero by construction, so the quotient is
    undefined at the population median and explodes near it. On the log scale a
    difference of z-scores is exactly the quantity a ratio is reaching for, and it
    is finite everywhere.

        > 0  damage exceeds repair capacity   (decompensated; mutagenic risk state)
        ~ 0  damage matched by repair         (compensated)
        < 0  repair capacity exceeds damage   (reserve)
    """
    return odi_z - rci_z


# ----------------------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------------------
def write_csv(path, rows, fieldnames=None):
    if not rows:
        return
    fieldnames = fieldnames or list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def main():
    os.makedirs(RESULTS, exist_ok=True)
    adducts = load("dna_damage_adducts.csv")
    flow = load("ddr_flow_panel.csv")
    cyto = load("ddr_cytogenetic_assays.csv")
    clin = load("repair_clinical_tests.csv")
    genes = load("repair_expression_panel.csv")

    lines = []
    def say(s=""):
        lines.append(s)
        print(s)

    say("=" * 78)
    say("MULTIPLEX PANEL FOR DNA STRESS AND DNA REPAIR CAPACITY")
    say("build and verification report")
    say("=" * 78)

    # ---- Module A -------------------------------------------------------------
    say()
    say("MODULE A - oxidative and alkylation adducts (LC-MS/MS, single MRM method)")
    say("-" * 78)
    va = validate_module_a(adducts)
    write_csv(os.path.join(RESULTS, "moduleA_transition_validation.csv"), va)

    nfail = sum(1 for r in va if r["product_check"] == "FAIL") + \
            sum(1 for r in va if r["is_check"].startswith("FAIL"))
    n_lit = sum(1 for r in va if r["product_check"] == "NOT_APPLICABLE")
    n_m3 = sum(1 for r in va if r["is_check"] == "PASS_VERIFY_M3")
    say("%-11s %9s %9s %9s %-7s  %-20s %7s %-14s" %
        ("analyte", "precursor", "prod_calc", "prod_decl", "product", "internal standard",
         "shift", "IS check"))
    ABBREV = {"PASS": "PASS", "FAIL": "FAIL", "NOT_APPLICABLE": "lit."}
    for r in va:
        say("%-11s %9.4f %9s %9.2f %-7s  %-20s %+7.3f %-14s" %
            (r["abbrev"], r["precursor_mz_calc"],
             ("%.4f" % r["product_mz_calc"]) if r["product_mz_calc"] != "" else "-",
             r["product_mz_declared"], ABBREV[r["product_check"]],
             r["internal_standard"], r["is_mass_shift_Da"], r["is_check"]))
    say()
    say("  glycosidic losses used: deoxyribose %.4f Da, ribose %.4f Da"
        % (LOSS_DEOXYRIBOSE, LOSS_RIBOSE))
    say("  transitions recomputed from formula: %d of %d analytes verifiable this way"
        % (len(va) - n_lit, len(va)))
    say("  %d free bases / non-nucleosides carry literature transitions marked 'lit.'" % n_lit)
    say("  failures: %d    internal standards at +3 needing M+3 bleed check: %d"
        % (nfail, n_m3))

    col = check_transition_collisions(va)
    say()
    if col:
        say("  UNRESOLVABLE TRANSITION COLLISIONS (would require chromatographic separation):")
        for a, b, dp, dq in col:
            say("    %s vs %s  dPrecursor=%.3f dProduct=%.3f" % (a, b, dp, dq))
    else:
        say("  transition collisions: none - every analyte has a unique precursor/product pair")

    iso = check_isobaric_precursors(va)
    if iso:
        say("  shared-ion pairs requiring chromatographic resolution or distinct dwell:")
        for a, b, why in iso:
            say("    %-10s %-10s %s" % (a, b, why))

    matrices = Counter(a["matrix"] for a in adducts)
    say()
    say("  matrix partitioning:")
    for m, n in matrices.most_common():
        say("    %-26s %d analytes" % (m, n))

    # ---- Module B -------------------------------------------------------------
    say()
    say("MODULE B - cellular damage response (one stained tube, flow cytometry)")
    say("-" * 78)
    spec, per_laser = validate_flow_panel(flow)
    write_csv(os.path.join(RESULTS, "moduleB_spectral_check.csv"), spec)
    say("  %d parameters across %d lasers: %s"
        % (len(flow), len(per_laser),
           ", ".join("%snm:%d" % (k, v) for k, v in sorted(per_laser.items(), key=lambda x: int(x[0])))))
    crit = [s for s in spec if s["status"] != "OK"]
    for s in spec:
        mark = " <-- spillover critical" if s["status"] != "OK" else ""
        say("    %4snm  %3dnm gap  %-42s -> %s%s"
            % (s["laser_nm"], s["gap_nm"], s["channel_low"].split(" (")[0],
               s["channel_high"].split(" (")[0], mark))
    say("  channels needing single-stain compensation controls: %d" % len(crit))
    say("  mandatory non-damage channels present: %s"
        % ", ".join(f["target"] for f in flow if f["detector_role"] in ("exclusion gate", "cell-cycle gate")))

    say()
    say("  cytogenetic endpoints (separate culture, same blood draw): %d" % len(cyto))
    for c in cyto:
        say("    %-42s %s" % (c["assay"][:42], c["role_in_panel"]))

    # ---- Module C -------------------------------------------------------------
    say()
    say("MODULE C - repair capacity")
    say("-" * 78)
    say("  established clinical tests carried into the panel: %d" % len(clin))
    for c in clin:
        say("    %-38s %-22s %s" % (c["test"][:38], c["pathway"][:22], c["role_in_panel"]))

    tiers, tier_counts = size_tiers(genes)
    write_csv(os.path.join(RESULTS, "moduleC_tier_sizing.csv"), tiers)
    say()
    say("  expression panel: %d targets (core %d, extended %d, normalisers %d)"
        % (len(genes), tier_counts["core"], tier_counts["extended"], tier_counts["normaliser"]))
    for t in tiers:
        say("    %-20s %3d targets -> %-32s spare wells: %s"
            % (t["tier"], t["n_targets"], t["fits_format"], t["spare_wells"]))

    bypath = Counter(g["pathway_code"] for g in genes)
    say()
    say("  pathway coverage:")
    PATHNAME = {"MMR": "mismatch repair", "BER": "base excision repair",
                "NER": "nucleotide excision repair", "HR": "homologous recombination",
                "NHEJ": "non-homologous end joining", "DR": "direct damage reversal",
                "ICL": "interstrand crosslink repair", "TLS": "translesion synthesis",
                "SIG": "damage signalling", "REF": "reference / normaliser"}
    for code, n in sorted(bypath.items(), key=lambda x: -x[1]):
        say("    %-5s %-30s %2d genes" % (code, PATHNAME.get(code, code), n))

    # ---- damage-repair coupling ----------------------------------------------
    say()
    say("DAMAGE-TO-REPAIR COUPLING")
    say("-" * 78)
    pairs = build_damage_repair_pairs(adducts, genes)
    write_csv(os.path.join(RESULTS, "damage_repair_coupling.csv"), pairs)
    for p in pairs:
        say("  %-10s %-32s %-6s %s"
            % (p["damage_analyte"], p["lesion_class"][:32], p["pathways"], p["repair_genes"][:44]))
    unc = [p for p in pairs if p["status"] == "UNCOUPLED"]
    say()
    say("  damage analytes with a named enzymatic partner in Module C: %d of %d"
        % (len(pairs) - len(unc), len(pairs)))
    if unc:
        say("  uncoupled (interpreted at pathway level only): %s"
            % ", ".join(p["damage_analyte"] for p in unc))

    # ---- scoring demonstration ----------------------------------------------
    say()
    say("COMPOSITE SCORING FRAMEWORK (worked example on SYNTHETIC values)")
    say("-" * 78)
    say("  The numbers below are invented to demonstrate the arithmetic. They are")
    say("  not measurements and carry no reference-interval authority.")
    say()
    demo_analytes = ["8-oxodG", "8-oxoGuo", "Tg", "O6-medG", "M1dG"]
    demo_obs      = [5.8, 2.4, 0.9, 0.22, 0.14]      # nmol/mol creatinine (illustrative)
    demo_ref_med  = [3.2, 1.6, 0.8, 0.18, 0.11]
    demo_ref_iqr  = [1.8, 0.9, 0.5, 0.10, 0.07]
    demo_cv       = [0.06, 0.07, 0.12, 0.10, 0.15]

    zs = [robust_z([o], m, i)[0] for o, m, i in zip(demo_obs, demo_ref_med, demo_ref_iqr)]
    ws = precision_weights(demo_cv)
    odi = composite(zs, ws)
    say("  %-10s %9s %9s %9s %7s %7s" % ("analyte", "observed", "ref med", "robust z", "CV", "weight"))
    for name, o, m, z, cv, w in zip(demo_analytes, demo_obs, demo_ref_med, zs, demo_cv, ws):
        say("  %-10s %9.3f %9.3f %9.3f %6.0f%% %7.3f" % (name, o, m, z, cv * 100, w))
    say("  -> Oxidative Damage Index (ODI)      = %+.3f" % odi)

    demo_rci_pathways = {"MMR": +0.10, "BER": -1.45, "NER": -0.20, "HR": +0.05,
                         "NHEJ": -0.10, "DR": -0.90}
    rci = sum(demo_rci_pathways.values()) / len(demo_rci_pathways)
    for k, v in demo_rci_pathways.items():
        say("  RCI_%-5s = %+.3f" % (k, v))
    say("  -> Repair Capacity Index (RCI)       = %+.3f" % rci)
    dsr = stress_repair_index(odi, rci)
    say("  -> Stress-to-Repair Index (ODI - RCI)= %+.3f" % dsr)
    state = ("decompensated - damage exceeds repair" if dsr > 0.5 else
             "compensated" if dsr > -0.5 else "repair reserve")
    say("  -> interpretation: %s" % state)
    say()
    say("  Pathway resolution is what makes this actionable: the composite is driven")
    say("  by RCI_BER (%+.2f) and RCI_DR (%+.2f), and the coupling table above names"
        % (demo_rci_pathways["BER"], demo_rci_pathways["DR"]))
    say("  OGG1/MUTYH/NUDT1 and MGMT as the specific enzymes to follow up. A single")
    say("  scalar damage score would not have localised the defect.")

    write_csv(os.path.join(RESULTS, "scoring_worked_example.csv"), [{
        "analyte": n, "observed_illustrative": o, "ref_median_illustrative": m,
        "robust_z": round(z, 4), "assay_cv": cv, "precision_weight": round(w, 4),
    } for n, o, m, z, cv, w in zip(demo_analytes, demo_obs, demo_ref_med, zs, demo_cv, ws)])

    # ---- summary -------------------------------------------------------------
    say()
    say("=" * 78)
    say("BUILD SUMMARY")
    say("=" * 78)
    total_checks = len(va) * 2 + 1 + 1
    say("  Module A analytes                     %d (incl. %d normalisers)"
        % (len(adducts), sum(1 for a in adducts if a["lesion_class"] == "normaliser")))
    say("  Module A transition failures          %d" % nfail)
    say("  Module A transition collisions        %d" % len(col))
    say("  Module B flow parameters              %d" % len(flow))
    say("  Module B spillover-critical pairs     %d" % len(crit))
    say("  Module B cytogenetic endpoints        %d" % len(cyto))
    say("  Module C clinical tests               %d" % len(clin))
    say("  Module C expression targets           %d" % len(genes))
    say("  Repair pathways covered               %d" % len([k for k in bypath if k != "REF"]))
    say("  Damage analytes enzymatically coupled %d of %d" % (len(pairs) - len(unc), len(pairs)))
    verdict = "ALL CHECKS PASSED" if (nfail == 0 and not col) else "CHECKS FAILED - see above"
    say("  Verdict                               %s" % verdict)

    with open(os.path.join(RESULTS, "panel_build_report.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0 if (nfail == 0 and not col) else 1


if __name__ == "__main__":
    raise SystemExit(main())
