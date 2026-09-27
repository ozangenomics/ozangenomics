"""Process plate-format ELISA optical densities: 4PL fit, back-calculation, QC flags.

Developed and tested against the SIMULATED dataset in data/simulated_elisa/. When run
on that input every output file is prefixed SIMULATED_ and carries '# SIMULATED DATA'
header lines, because results derived from simulated input are themselves simulated.

Per plate:
  1. fit a four-parameter logistic to the standards (Levenberg-Marquardt, pure Python)
  2. back-calculate every well:  x = C * ((A - D) / (y - D) - 1) ** (1 / B)
  3. per-sample mean, duplicate CV, flags (below LLOQ, above ULOQ, CV, saturation,
     outside curve asymptotes), standard recovery, QC recovery
  4. if a ground-truth file is present, score the back-calculated unknowns against it

Usage:
    python3 src/process_elisa.py [--input data/simulated_elisa/SIMULATED_elisa_long.csv]
                                 [--out results] [--cv-limit 15] [--recovery 80 120]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics as st
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_INPUT = os.path.join(ROOT, "data", "simulated_elisa", "SIMULATED_elisa_long.csv")
DEFAULT_OUT = os.path.join(ROOT, "results")


# ---------------------------------------------------------------------------
# 4PL model and fitter
# ---------------------------------------------------------------------------
def four_pl(x: float, A: float, B: float, C: float, D: float) -> float:
    if x <= 0:
        return A
    return D + (A - D) / (1.0 + (x / C) ** B)


def inverse_four_pl(y: float, A: float, B: float, C: float, D: float):
    """Return concentration, or None if y lies outside the open interval (A, D)."""
    lo, hi = min(A, D), max(A, D)
    if not (lo < y < hi):
        return None
    return C * ((A - D) / (y - D) - 1.0) ** (1.0 / B)


def _solve(mat, vec):
    """Gaussian elimination with partial pivoting for small dense systems."""
    n = len(vec)
    M = [row[:] + [vec[i]] for i, row in enumerate(mat)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(M[r][i]))
        M[i], M[piv] = M[piv], M[i]
        if abs(M[i][i]) < 1e-14:
            raise ZeroDivisionError("singular normal matrix")
        for r in range(n):
            if r != i:
                f = M[r][i] / M[i][i]
                for c in range(i, n + 1):
                    M[r][c] -= f * M[i][c]
    return [M[i][n] / M[i][i] for i in range(n)]


WEIGHTS = {"none": lambda y: 1.0, "1/y": lambda y: 1.0 / max(y, 1e-3),
           "1/y2": lambda y: 1.0 / max(y, 1e-3) ** 2}


def fit_four_pl(xs, ys, weighting="1/y2", max_iter=500, tol=1e-10):
    """Weighted least-squares 4PL fit by Levenberg-Marquardt.
    Weights are 1, 1/y or 1/y^2 on the observed response; 1/y^2 is the usual
    choice for immunoassays because response variance scales with signal.
    B and C are fitted in log space so they stay positive."""
    wfun = WEIGHTS[weighting]
    sw = [math.sqrt(wfun(y)) for y in ys]

    def unpack(q):
        return q[0], math.exp(q[1]), math.exp(q[2]), q[3]

    def residuals(q):
        A, B, C, D = unpack(q)
        return [w * (four_pl(x, A, B, C, D) - y) for x, y, w in zip(xs, ys, sw)]

    pos = [x for x in xs if x > 0]
    q = [min(ys), 0.0, math.log(math.exp(st.mean(math.log(x) for x in pos))), max(ys)]
    r = residuals(q)
    sse = sum(v * v for v in r)
    lam = 1e-3
    for _ in range(max_iter):
        # numerical Jacobian
        J = []
        h = 1e-6
        for k in range(4):
            qk = q[:]
            qk[k] += h
            rk = residuals(qk)
            J.append([(a - b) / h for a, b in zip(rk, r)])  # column k
        JTJ = [[sum(J[i][m] * J[j][m] for m in range(len(r))) for j in range(4)] for i in range(4)]
        JTr = [sum(J[i][m] * r[m] for m in range(len(r))) for i in range(4)]
        improved = False
        for _ in range(30):
            Amat = [[JTJ[i][j] + (lam * JTJ[i][i] if i == j else 0.0) for j in range(4)] for i in range(4)]
            try:
                step = _solve(Amat, [-v for v in JTr])
            except ZeroDivisionError:
                lam *= 10
                continue
            q_new = [a + b for a, b in zip(q, step)]
            r_new = residuals(q_new)
            sse_new = sum(v * v for v in r_new)
            if sse_new < sse:
                q, r, sse = q_new, r_new, sse_new
                lam = max(lam / 10, 1e-12)
                improved = True
                break
            lam *= 10
        if not improved or sse < tol:
            break
    A, B, C, D = unpack(q)
    # report R^2 on the unweighted scale so it is comparable across weightings
    raw_sse = sum((four_pl(x, A, B, C, D) - y) ** 2 for x, y in zip(xs, ys))
    ybar = st.mean(ys)
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    r2 = 1.0 - raw_sse / ss_tot if ss_tot > 0 else float("nan")
    return {"A": A, "B": B, "C": C, "D": D, "sse": raw_sse, "r2": r2, "n": len(xs),
            "weighting": weighting}


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------
def read_commented_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        header = [l for l in fh if l.startswith("#")]
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    return header, rows


def write_csv(path, fieldnames, rows, header_lines):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        for line in header_lines:
            fh.write(line.rstrip("\n") + "\n")
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def fmt(v, nd=3):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    return f"{v:.{nd}f}"


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--input", default=DEFAULT_INPUT)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--cv-limit", type=float, default=15.0, help="duplicate CV%% flag limit")
    ap.add_argument("--recovery", type=float, nargs=2, default=(80.0, 120.0),
                    help="acceptable recovery window %% for standards and QC")
    ap.add_argument("--saturation", type=float, default=4.0)
    ap.add_argument("--weighting", choices=sorted(WEIGHTS), default="1/y2",
                    help="least-squares weighting on the response (default 1/y2)")
    ap.add_argument("--qc-targets", default=None,
                    help="JSON file with qc_levels_pg_mL; defaults to the parameters file "
                         "beside the input if present")
    ap.add_argument("--ground-truth", default=None,
                    help="ground-truth CSV for scoring; defaults to the file beside the input")
    args = ap.parse_args()

    header_lines, wells = read_commented_csv(args.input)
    simulated = any("SIMULATED" in l for l in header_lines) or \
        any(r.get("data_status") == "SIMULATED" for r in wells)
    prefix = "SIMULATED_" if simulated else ""
    out_header = ([l for l in header_lines if "SIMULATED" in l or "No instrument" in l] +
                  ["# Derived by src/process_elisa.py from " + os.path.relpath(args.input, ROOT)])
    os.makedirs(args.out, exist_ok=True)

    in_dir = os.path.dirname(args.input)
    qc_path = args.qc_targets or os.path.join(in_dir, prefix + "elisa_parameters.json")
    qc_targets = {}
    if os.path.exists(qc_path):
        with open(qc_path, encoding="utf-8") as fh:
            qc_targets = json.load(fh).get("qc_levels_pg_mL", {})
    gt_path = args.ground_truth or os.path.join(in_dir, prefix + "elisa_ground_truth.csv")

    rec_lo, rec_hi = args.recovery
    by_plate = defaultdict(list)
    for r in wells:
        r["od_f"] = float(r["od"])
        by_plate[int(r["plate"])].append(r)

    fit_rows, std_rows, sample_rows = [], [], []
    for plate in sorted(by_plate):
        pw = by_plate[plate]
        stds = [r for r in pw if r["role"] == "standard"]
        xs = [float(r["nominal_pg_mL"]) for r in stds]
        ys = [r["od_f"] for r in stds]
        fit = fit_four_pl(xs, ys, args.weighting)
        A, B, C, D = fit["A"], fit["B"], fit["C"], fit["D"]
        lloq, uloq = min(xs), max(xs)
        blanks = [r["od_f"] for r in pw if r["role"] == "blank"]
        blank_mean = st.mean(blanks) if blanks else float("nan")
        fit_rows.append({
            "plate": plate, "A": fmt(A, 4), "B": fmt(B, 4), "C": fmt(C, 2), "D": fmt(D, 4),
            "weighting": fit["weighting"], "r2": fmt(fit["r2"], 5), "sse": fmt(fit["sse"], 5),
            "n_standards": fit["n"],
            "lloq_pg_mL": fmt(lloq, 2), "uloq_pg_mL": fmt(uloq, 1),
            "blank_mean_od": fmt(blank_mean), "blank_vs_A_od": fmt(blank_mean - A),
        })

        # standard recovery (per well)
        for r in stds:
            nominal = float(r["nominal_pg_mL"])
            bc = inverse_four_pl(r["od_f"], A, B, C, D)
            rec = None if bc is None else 100.0 * bc / nominal
            std_rows.append({
                "plate": plate, "well": r["well"], "standard": r["sample_id"],
                "nominal_pg_mL": fmt(nominal, 2), "od": fmt(r["od_f"]),
                "back_calc_pg_mL": fmt(bc, 2), "recovery_pct": fmt(rec, 1),
                "flag": "" if rec is not None and rec_lo <= rec <= rec_hi else "RECOVERY",
            })

        # samples: group by sample_id
        groups = defaultdict(list)
        for r in pw:
            if r["role"] in ("unknown", "qc", "blank"):
                groups[r["sample_id"]].append(r)
        for sid, rs in groups.items():
            ods = [r["od_f"] for r in rs]
            mean_od = st.mean(ods)
            cv = 100.0 * st.stdev(ods) / mean_od if len(ods) > 1 and mean_od > 0 else float("nan")
            flags = []
            if any(o >= args.saturation for o in ods):
                flags.append("SATURATED")
            role = rs[0]["role"]
            # CV is not meaningful on the blank's near-zero signal
            if role != "blank" and len(ods) > 1 and not math.isnan(cv) and cv > args.cv_limit:
                flags.append("CV")
            conc = inverse_four_pl(mean_od, A, B, C, D)
            if role != "blank":
                if conc is None:
                    flags.append("BELOW_A" if mean_od <= min(A, D) else "ABOVE_D")
                elif conc < lloq:
                    flags.append("BELOW_LLOQ")
                elif conc > uloq:
                    flags.append("ABOVE_ULOQ")
            rec = None
            if role == "qc" and sid in qc_targets and conc is not None:
                rec = 100.0 * conc / float(qc_targets[sid])
                if not rec_lo <= rec <= rec_hi:
                    flags.append("QC_RECOVERY")
            reportable = conc is not None and not any(
                f in flags for f in ("SATURATED", "BELOW_LLOQ", "ABOVE_ULOQ", "BELOW_A", "ABOVE_D"))
            sample_rows.append({
                "plate": plate, "sample_id": sid, "role": role,
                "wells": ";".join(r["well"] for r in rs), "n": len(rs),
                "mean_od": fmt(mean_od), "cv_pct": fmt(cv, 1),
                "conc_pg_mL": fmt(conc, 2), "reportable": "yes" if reportable else "no",
                "qc_target_pg_mL": fmt(float(qc_targets[sid]), 1) if sid in qc_targets else "",
                "qc_recovery_pct": fmt(rec, 1), "flags": ";".join(flags),
            })

    write_csv(os.path.join(args.out, prefix + "elisa_fit_parameters.csv"),
              list(fit_rows[0].keys()), fit_rows, out_header)
    write_csv(os.path.join(args.out, prefix + "elisa_standards_recovery.csv"),
              list(std_rows[0].keys()), std_rows, out_header)
    write_csv(os.path.join(args.out, prefix + "elisa_results.csv"),
              list(sample_rows[0].keys()), sample_rows, out_header)

    # ------------------------------------------------------------------ console
    print(f"{'SIMULATED input - ' if simulated else ''}processed {len(wells)} wells on "
          f"{len(by_plate)} plate(s) from {os.path.relpath(args.input, ROOT)}")
    print(f"\n4PL fits, weighting {args.weighting} (OD = D + (A-D)/(1+(x/C)^B)):")
    for f in fit_rows:
        print(f"  plate {f['plate']}: A={f['A']} B={f['B']} C={f['C']} D={f['D']} "
              f"R2={f['r2']}  blank-A={f['blank_vs_A_od']}")
    n_std_flag = sum(1 for s in std_rows if s["flag"])
    print(f"\nstandards outside {rec_lo:.0f}-{rec_hi:.0f}% recovery: {n_std_flag}/{len(std_rows)}")
    for s in std_rows:
        if s["flag"]:
            print(f"  plate {s['plate']} {s['well']} {s['standard']} recovery {s['recovery_pct']}%")
    qcs = [s for s in sample_rows if s["role"] == "qc"]
    print("\nQC recovery:")
    for q in qcs:
        print(f"  plate {q['plate']} {q['sample_id']:8s} {q['conc_pg_mL']:>9s} pg/mL "
              f"({q['qc_recovery_pct']}%) {q['flags']}")
    unk = [s for s in sample_rows if s["role"] == "unknown"]
    flag_counts = defaultdict(int)
    for s in unk:
        for f in s["flags"].split(";"):
            if f:
                flag_counts[f] += 1
    print(f"\nunknowns: {len(unk)}, reportable: {sum(1 for s in unk if s['reportable'] == 'yes')}")
    for k, v in sorted(flag_counts.items()):
        print(f"  {k}: {v}")

    # ------------------------------------------------------------------ scoring
    if os.path.exists(gt_path):
        _, truth = read_commented_csv(gt_path)
        tmap = {(int(t["plate"]), t["sample_id"]): float(t["true_pg_mL"]) for t in truth}
        score_rows, biases = [], []
        range_correct = range_total = 0
        for s in unk:
            key = (s["plate"], s["sample_id"])
            if key not in tmap:
                continue
            true = tmap[key]
            est = float(s["conc_pg_mL"]) if s["conc_pg_mL"] else None
            plate_fit = next(f for f in fit_rows if f["plate"] == s["plate"])
            in_range_true = float(plate_fit["lloq_pg_mL"]) <= true <= float(plate_fit["uloq_pg_mL"])
            bias = None if est is None else 100.0 * (est - true) / true
            if s["reportable"] == "yes" and bias is not None:
                biases.append(bias)
            range_total += 1
            range_correct += int(in_range_true == (s["reportable"] == "yes"))
            score_rows.append({
                "plate": s["plate"], "sample_id": s["sample_id"],
                "true_pg_mL": fmt(true, 2), "est_pg_mL": s["conc_pg_mL"],
                "bias_pct": fmt(bias, 1), "reportable": s["reportable"],
                "true_in_range": "yes" if in_range_true else "no", "flags": s["flags"],
            })
        write_csv(os.path.join(args.out, prefix + "elisa_scoring.csv"),
                  list(score_rows[0].keys()), score_rows, out_header)
        abs_b = [abs(b) for b in biases]
        print(f"\nscoring against {os.path.relpath(gt_path, ROOT)}:")
        print(f"  reportable unknowns scored: {len(biases)}")
        print(f"  median |bias|: {st.median(abs_b):.1f}%   90th pct |bias|: "
              f"{sorted(abs_b)[int(0.9 * (len(abs_b) - 1))]:.1f}%   max |bias|: {max(abs_b):.1f}%")
        print(f"  within ±20%: {sum(1 for b in abs_b if b <= 20)}/{len(abs_b)}")
        print(f"  mean signed bias: {st.mean(biases):+.1f}%")
        print(f"  range classification (reportable == truly within LLOQ-ULOQ): "
              f"{range_correct}/{range_total}")
    else:
        print("\nno ground-truth file found; scoring skipped")


if __name__ == "__main__":
    main()
