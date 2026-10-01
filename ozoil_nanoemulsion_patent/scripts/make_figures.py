"""Regenerate FIG. 1–9 of the PCT draft from the source workbook.

Inputs : the formulation workbook (sheets 'PSA', 'Distributions-LD', and the six
         DLS metafiles in 'Distribution-DLS'), unzipped so that the EMF files are
         available as emf/image1.emf ... image6.emf.
Outputs: figures/FIG1 ... FIG9 PNG files and data/DLS_curves_extracted.csv.

The DLS curves are recovered from the Malvern EMF metafiles by parsing the
POLYLINE16 records and mapping pixel coordinates to the chart axes
(x: log10 scale, grid lines at 10/100/1000/10000 nm; y: linear, baseline = 0 %).
"""
import struct, json, csv, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import openpyxl

XLSX = sys.argv[1] if len(sys.argv) > 1 else 'formulation.xlsx'
EMF_DIR = sys.argv[2] if len(sys.argv) > 2 else 'emf'
OUT = sys.argv[3] if len(sys.argv) > 3 else 'figures'

def emf_polylines(path):
    data = open(path, 'rb').read(); off = 0; polys = []
    while off < len(data):
        t, sz = struct.unpack_from('<II', data, off); body = data[off+8:off+sz]
        if t == 87:  # EMR_POLYLINE16
            n = struct.unpack_from('<I', body, 16)[0]
            pts = struct.unpack_from('<%dh' % (2*n), body, 20)
            polys.append([(pts[k], pts[k+1]) for k in range(0, 2*n, 2)])
        if sz == 0: break
        off += sz
    return [p for p in polys if len(p) == 49]  # 49-point size-distribution curves

def xmap(px): return 10 ** (1 + (px - 82) / ((614 - 82) / 3))
ymax = {1: 30, 2: 20, 3: 20, 4: 25, 5: 25, 6: 25}; ybase = {1: 175, 2: 175, 3: 175, 4: 189, 5: 189, 6: 189}
labels5 = ['1 pass', '3 passes', '5 passes', '7 passes', '10 passes']; labels3 = labels5[:3]
tests = {1: ('20260729A', 'Ex. 1 (F1: 1.2% lecithin / 0.5% PS80), F12Y–H30Z, 30,000 psi', labels5),
         2: ('20260729B', 'Ex. 2 (F1), F12Y–H30Z, 20,000 psi', labels5),
         3: ('20260729C', 'Ex. 3 (F1), F20Y–H30Z, 27,500 psi', labels5),
         4: ('20260730A', 'Ex. 4 (F2: 1.2% lecithin / 1.0% PS80), F12Y–H30Z, 30,000 psi', labels3),
         5: ('20260803A', 'Ex. 5 (F3: 2.4% lecithin / 0.5% PS80), F12Y–H30Z, 30,000 psi', labels3),
         6: ('20260803B', 'Ex. 6 (F4: 2.4% lecithin / 1.0% PS80), F12Y–H30Z, 30,000 psi', labels3)}
styles = [('k', '-'), ('k', '--'), ('k', '-.'), ('k', ':'), ('0.5', '-'), ('0.5', '--')]

rows_out = []
for i in range(1, 7):
    series = [([xmap(x) for x, y in p], [(ybase[i]-y)/(ybase[i]-49)*ymax[i] for x, y in p]) for p in emf_polylines(f'{EMF_DIR}/image{i}.emf')]
    fig, ax = plt.subplots(figsize=(7.5, 4))
    for (xs, ys), lab, (c, ls) in zip(series, tests[i][2], styles):
        ax.plot(xs, ys, color=c, ls=ls, lw=1.6, label=lab)
        rows_out += [(tests[i][0], lab, round(x, 1), round(y, 2)) for x, y in zip(xs, ys)]
    ax.set_xscale('log'); ax.set_xlim(10, 10000); ax.set_ylim(0, ymax[i]+2)
    ax.set_xlabel('Hydrodynamic diameter (d, nm)'); ax.set_ylabel('Intensity (%)')
    ax.set_title(f'FIG. {i+1}  DLS intensity size distribution — Test {tests[i][0]}\n{tests[i][1]}', fontsize=10)
    ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(f'{OUT}/FIG{i+1}_DLS_{tests[i][0]}.png', dpi=220); plt.close(fig)
with open('data/DLS_curves_extracted.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['test', 'passes', 'diameter_nm', 'intensity_pct']); w.writerows(rows_out)

wb = openpyxl.load_workbook(XLSX, data_only=True); ws = wb['Distributions-LD']
rows = list(ws.iter_rows(min_row=2, values_only=True)); d = np.array([r[0] for r in rows])
cols = {ws.cell(1, c).value: np.array([r[c-1] for r in rows]) for c in range(2, 6)}
ldlab = {'20260729A-Unp': 'F1 (1.2% lecithin / 0.5% PS80) — coarse pre-emulsion', '20260730A-Unp': 'F2 (1.2% lecithin / 1.0% PS80) — coarse pre-emulsion',
         '20260803A-Unp': 'F3 (2.4% lecithin / 0.5% PS80) — coarse pre-emulsion', '20260803B-Unp': 'F4 (2.4% lecithin / 1.0% PS80) — coarse pre-emulsion'}
fig, ax = plt.subplots(figsize=(7.5, 4))
for (k, v), (c, ls) in zip(cols.items(), styles): ax.plot(d, v, color=c, ls=ls, lw=1.6, label=ldlab[k])
ax.set_xscale('log'); ax.set_xlim(0.01, 100); ax.set_xlabel('Particle diameter (µm)'); ax.set_ylabel('Volume frequency (%)')
ax.set_title('FIG. 1  Laser-diffraction volume size distribution of unprocessed coarse pre-emulsions\n(Horiba LA-950, RI = 1.47 + 0.01i)', fontsize=10)
ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig(f'{OUT}/FIG1_LD_unprocessed.png', dpi=220); plt.close(fig)

ws = wb['PSA']; data = {}; cur = None
for r in ws.iter_rows(min_row=4, values_only=True):
    if r[1]: cur = r[1]; data[cur] = []
    if isinstance(r[6], (int, float)) and isinstance(r[11], (int, float)): data[cur].append((r[6], r[11], r[12]))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.8)); mk = ['o', 's', '^', 'D', 'v', 'P']
for (k, v), m, (c, ls) in zip(data.items(), mk, styles):
    a1.plot([x[0] for x in v], [x[1] for x in v], marker=m, color=c, ls=ls, lw=1.3, label=k)
    a2.plot([x[0] for x in v], [x[2] for x in v], marker=m, color=c, ls=ls, lw=1.3, label=k)
a1.set_xlabel('Number of passes'); a1.set_ylabel('Z-average (nm)'); a1.set_title('(a) Z-average vs. passes', fontsize=10); a1.grid(alpha=0.3)
a2.set_xlabel('Number of passes'); a2.set_ylabel('PdI'); a2.set_title('(b) Polydispersity index vs. passes', fontsize=10); a2.grid(alpha=0.3); a2.axhline(0.1, color='0.6', ls=':', lw=1)
a2.legend(frameon=False, fontsize=8)
fig.suptitle('FIG. 8  Evolution of Z-average and PdI with number of Microfluidizer passes (Tests 20260729A–20260803B)', fontsize=10)
fig.tight_layout(); fig.savefig(f'{OUT}/FIG8_Zavg_PdI_vs_passes.png', dpi=220); plt.close(fig)
print('done')
