"""Assemble the lymphocyte marker panel into a single Excel workbook.

Reads the pipeline outputs in results/ and writes
results/DNA_damage_lymphocyte_panel.xlsx with one sheet per view:

    Overview          what the tiers mean, live counts, method notes
    Measurable panel  the L1 and L2 proteins with their lead peptide
    All markers       every protein with its top three peptides
    Phospho markers   canonical DDR sites, with tryptic suitability
    PRM transitions   precursor and product ions for the measurable set
    Controls          proteins that are not analytes here, and why

Run src/panel_build.py lymphocyte first.
"""
from __future__ import annotations

import csv
import os

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FONT = 'Arial'

# Tier colours: blue for the two measurable tiers, amber for enrichment,
# grey for what cannot be measured in resting lymphocytes.
TIER_FILL = {
    'L1': 'C6DBEF', 'L2': 'E3EEF9', 'L3': 'FDEBD0',
    'L4-stim': 'EAEAEA', 'L4-absent': 'D9D9D9',
}
HEAD_FILL = PatternFill('solid', fgColor='1F3864')
HEAD_FONT = Font(name=FONT, size=10, bold=True, color='FFFFFF')
THIN = Side(style='thin', color='BFBFBF')
BORDER = Border(bottom=THIN)


def path(*p):
    return os.path.join(HERE, *p)


def read(name):
    with open(path('results', name)) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_sheet(ws, headers, rows, widths, numfmt=None, tier_col=None,
                wrap_cols=()):
    """Write a header row and body, with filter, freeze pane and styling."""
    numfmt = numfmt or {}
    ws.append(headers)
    for c, _h in enumerate(headers, 1):
        cell = ws.cell(row=1, column=c)
        cell.fill, cell.font = HEAD_FILL, HEAD_FONT
        cell.alignment = Alignment(vertical='center', wrap_text=True)
    ws.row_dimensions[1].height = 28

    for r in rows:
        ws.append(r)
    for r in range(2, ws.max_row + 1):
        tier = ws.cell(row=r, column=tier_col).value if tier_col else None
        fill = (PatternFill('solid', fgColor=TIER_FILL[tier])
                if tier in TIER_FILL else None)
        for c in range(1, len(headers) + 1):
            cell = ws.cell(row=r, column=c)
            cell.font = Font(name=FONT, size=10)
            cell.border = BORDER
            if fill:
                cell.fill = fill
            if c in numfmt:
                cell.number_format = numfmt[c]
            if c in wrap_cols:
                cell.alignment = Alignment(wrap_text=True, vertical='top')
            else:
                cell.alignment = Alignment(vertical='top')
    for c, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(c)].width = w
    ws.freeze_panes = 'A2'
    if ws.max_row > 1:
        ws.auto_filter.ref = (f'A1:{get_column_letter(len(headers))}'
                              f'{ws.max_row}')


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


def main():
    markers = read('lymphocyte_marker_list.tsv')
    phos = read('lymphocyte_phospho_markers.tsv')
    trans = read('lymphocyte_prm_transitions.tsv')
    prot = read('lymphocyte_panel_proteins.tsv')

    PATHWAY = {
        'DDR_signalling': 'DNA damage response signalling',
        'DSB_repair': 'Double-strand break repair',
        'BER': 'Base excision repair',
        'BER_oxidative': 'Oxidative base excision repair',
        'NER': 'Nucleotide excision repair',
        'MMR': 'Mismatch repair',
        'nuclear_DAMP': 'Chromatin and nuclear DAMPs',
        'oxidative_stress': 'Oxidative stress and redox',
        'stress_chaperone': 'Stress chaperones',
        'mtDNA_maintenance': 'Mitochondrial DNA maintenance',
        'lymphocyte_deaminase': 'Lymphocyte deaminases',
        'lymphocyte_programmed_breaks': 'Programmed breaks, development',
    }
    wb = Workbook()

    # ---- Overview ---------------------------------------------------------
    ws = wb.active
    ws.title = 'Overview'
    ws.sheet_view.showGridLines = False
    ws['A1'] = 'DNA damage, oxidation and repair marker panel'
    ws['A1'].font = Font(name=FONT, size=16, bold=True, color='1F3864')
    ws['A2'] = 'Isolated peripheral blood lymphocytes, LC-MS/MS'
    ws['A2'].font = Font(name=FONT, size=11, italic=True, color='595959')

    rowi = 4
    for line in [
        'Tier meanings',
    ]:
        ws.cell(row=rowi, column=1, value=line).font = Font(
            name=FONT, size=11, bold=True)
        rowi += 1
    tiers = [
        ('L1', 'Abundant. Whole-cell lysate, DIA or PRM, no enrichment.'),
        ('L2', 'Deep or fractionated proteome, or targeted PRM.'),
        ('L3', 'Needs nuclear fractionation or immunoenrichment.'),
        ('L4-stim', 'Absent from resting cells. Needs mitogen stimulation.'),
        ('L4-absent', 'Not a lymphocyte protein. Use as a control marker.'),
    ]
    tier_counts = {}
    for m in markers:
        if m['peptide_rank'] == '1':
            tier_counts[m['tier']] = tier_counts.get(m['tier'], 0) + 1
    for label, desc in tiers:
        ws.cell(row=rowi, column=1, value=label).fill = PatternFill(
            'solid', fgColor=TIER_FILL[label])
        ws.cell(row=rowi, column=1).font = Font(name=FONT, size=10, bold=True)
        ws.cell(row=rowi, column=2, value=desc).font = Font(name=FONT, size=10)
        # Static counts, not formulas: LibreOffice is unavailable in the
        # build environment, so a formula could not be recalculated and would
        # read back as blank to pandas and most previewers. The All markers
        # sheet is the source of truth; its autofilter reproduces these.
        ws.cell(row=rowi, column=3,
                value=tier_counts.get(label, 0)).font = Font(
            name=FONT, size=10, bold=True)
        ws.cell(row=rowi, column=4, value='proteins').font = Font(
            name=FONT, size=10, color='595959')
        rowi += 1
    ws.cell(row=rowi, column=2, value='Measurable total, L1 plus L2').font = (
        Font(name=FONT, size=10, bold=True))
    ws.cell(row=rowi, column=3,
            value=tier_counts.get('L1', 0) + tier_counts.get('L2', 0)).font = (
        Font(name=FONT, size=10, bold=True))
    rowi += 2

    notes = [
        ('Sheets', ''),
        ('Measurable panel',
         'The L1 and L2 proteins with their lead peptide. Start here.'),
        ('All markers',
         'Every protein with its top three peptides, ranked.'),
        ('Phospho markers',
         'Canonical damage-response sites, checked against the sequence.'),
        ('PRM transitions',
         'Precursor and product ions, light and heavy, for the measurable '
         'tiers.'),
        ('Controls',
         'Proteins that are not analytes in this matrix, and what they '
         'report instead.'),
        ('', ''),
        ('Method', ''),
        ('Sequences',
         'UniProtKB/Swiss-Prot via the NCBI BLAST distribution, release '
         '15 September 2026; 20,525 human entries.'),
        ('Digestion',
         'Trypsin, fully specific, no cleavage before proline, zero missed '
         'cleavages. The MaxQuant, Skyline and Spectronaut convention.'),
        ('Masses',
         'Monoisotopic. Cysteine fixed carbamidomethyl, +57.021464. '
         'Methionine oxidation flagged as a liability, not a separate '
         'precursor.'),
        ('Proteotypicity',
         'Tested by substring search against all 20,525 human entries, not '
         'only against other tryptic products.'),
        ('Scope',
         'protein_specific means unique to this protein. family_level means '
         'the peptide measures a protein family, named in shared_with. '
         'Core histone variants are family_level by necessity.'),
        ('', ''),
        ('Counts',
         'The tier counts above are fixed values computed when this file was '
         'built, not live formulas. Filter the All markers sheet on Peptide '
         'rank = 1 to reproduce them.'),
        ('Caveat',
         'Tier assignments are a curated reading of lymphocyte expression, '
         'localisation and cell-cycle dependence. They are a design prior, '
         'not a measurement. Confirm L2 and L3 calls in your own '
         'preparation with heavy standards before committing to a cohort.'),
    ]
    for label, text in notes:
        ws.cell(row=rowi, column=1, value=label).font = Font(
            name=FONT, size=10, bold=True)
        c = ws.cell(row=rowi, column=2, value=text)
        c.font = Font(name=FONT, size=10)
        c.alignment = Alignment(wrap_text=True, vertical='top')
        ws.row_dimensions[rowi].height = 30 if len(text) > 95 else 15
        rowi += 1
    for col, w in zip('ABCD', (22, 95, 8, 10)):
        ws.column_dimensions[col].width = w

    # ---- Measurable panel -------------------------------------------------
    lead = [m for m in markers if m['peptide_rank'] == '1']
    meas = [m for m in lead if m['tier'] in ('L1', 'L2')]
    rows = [[PATHWAY[m['pathway']], m['gene'], m['accession'],
             m['protein_name'], m['tier'], m['peptide'], m['residues'],
             num(m['mz_2plus']), num(m['mz_3plus']), m['scope'],
             m['liabilities'], m['expected_level'], m['assay'],
             m['confounder']] for m in meas]
    write_sheet(
        wb.create_sheet('Measurable panel'),
        ['Pathway', 'Gene', 'UniProt', 'Protein', 'Tier', 'Lead peptide',
         'Residues', 'm/z 2+', 'm/z 3+', 'Scope', 'Liabilities',
         'Expected level', 'Assay', 'Confounder'],
        rows, [30, 11, 10, 42, 7, 26, 11, 11, 11, 16, 24, 32, 34, 40],
        numfmt={8: '0.0000', 9: '0.0000'}, tier_col=5,
        wrap_cols=(12, 13, 14))

    # ---- All markers ------------------------------------------------------
    rows = [[PATHWAY[m['pathway']], m['gene'], m['accession'],
             m['protein_name'], m['tier'], m['expected_level'], m['assay'],
             m['confounder'], int(m['peptide_rank']), m['peptide'],
             m['residues'], num(m['mz_2plus']), num(m['mz_3plus']),
             m['scope'], m['shared_with'], m['liabilities']]
            for m in markers]
    write_sheet(
        wb.create_sheet('All markers'),
        ['Pathway', 'Gene', 'UniProt', 'Protein', 'Tier', 'Expected level',
         'Assay', 'Confounder', 'Peptide rank', 'Peptide', 'Residues',
         'm/z 2+', 'm/z 3+', 'Scope', 'Shared with', 'Liabilities'],
        rows,
        [30, 11, 10, 42, 10, 32, 34, 40, 8, 26, 11, 11, 11, 16, 34, 26],
        numfmt={12: '0.0000', 13: '0.0000'}, tier_col=5,
        wrap_cols=(6, 7, 8))

    # ---- Phospho markers --------------------------------------------------
    rows = [[p['gene'], p['accession'], p['site'], p['marker'], p['kinase'],
             p['peptide'], p['residues'],
             int(p['length']) if p['length'] else '',
             int(p['missed_cleavages']) if p['missed_cleavages'] != '' else '',
             num(p['unmod_mz_2plus']), num(p['phospho_mz_2plus']),
             num(p['phospho_mz_3plus']),
             int(p['STY_count']) if p['STY_count'] else '', p['verdict']]
            for p in phos]
    ws = wb.create_sheet('Phospho markers')
    write_sheet(
        ws,
        ['Gene', 'UniProt', 'Site', 'Marker', 'Kinase', 'Tryptic peptide',
         'Residues', 'Length', 'Missed cleav.', 'Unmod m/z 2+',
         'Phospho m/z 2+', 'Phospho m/z 3+', 'S/T/Y in peptide', 'Verdict'],
        rows, [11, 10, 8, 32, 16, 34, 12, 9, 12, 13, 14, 14, 14, 54],
        numfmt={10: '0.0000', 11: '0.0000', 12: '0.0000'},
        wrap_cols=(14,))
    good = PatternFill('solid', fgColor='D5E8D4')
    warn = PatternFill('solid', fgColor='FDEBD0')
    for r in range(2, ws.max_row + 1):
        v = ws.cell(row=r, column=14).value or ''
        f = good if v.startswith('good') else warn
        for c in range(1, 15):
            ws.cell(row=r, column=c).fill = f

    # ---- PRM transitions --------------------------------------------------
    rows = [[t['gene'], t['accession'], t['peptide'],
             t['quantification_scope'], int(t['precursor_charge']),
             num(t['precursor_mz_light']),
             num(t['precursor_mz_heavy']) if t['precursor_mz_heavy'] else '',
             t['heavy_label'], t['product_ion'], int(t['product_charge']),
             num(t['product_mz'])] for t in trans]
    write_sheet(
        wb.create_sheet('PRM transitions'),
        ['Gene', 'UniProt', 'Peptide', 'Scope', 'Prec. charge',
         'Prec. m/z light', 'Prec. m/z heavy', 'Heavy label', 'Product ion',
         'Prod. charge', 'Product m/z'],
        rows, [11, 10, 26, 17, 12, 15, 15, 26, 12, 12, 13],
        numfmt={6: '0.0000', 7: '0.0000', 11: '0.0000'})

    # ---- Controls ---------------------------------------------------------
    ctrl = [m for m in lead if m['tier'].startswith('L4')]
    rows = [[m['gene'], m['accession'], m['protein_name'], m['tier'],
             m['expected_level'], m['confounder'], m['assay'], m['peptide'],
             num(m['mz_2plus'])] for m in ctrl]
    write_sheet(
        wb.create_sheet('Controls'),
        ['Gene', 'UniProt', 'Protein', 'Tier', 'Expected level',
         'What it reports instead', 'Use', 'Peptide', 'm/z 2+'],
        rows, [11, 10, 44, 12, 34, 54, 46, 26, 11],
        numfmt={9: '0.0000'}, tier_col=4, wrap_cols=(5, 6, 7))

    out = path('results', 'DNA_damage_lymphocyte_panel.xlsx')
    wb.save(out)
    print(f'wrote {out}')
    print(f'  Measurable panel {len(meas)} | All markers {len(markers)} | '
          f'Phospho {len(phos)} | Transitions {len(trans)} | '
          f'Controls {len(ctrl)}')
    return out


if __name__ == '__main__':
    main()
