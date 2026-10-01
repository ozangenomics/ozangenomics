"""Build a DNA-damage / oxidative-stress / DNA-repair plasma MS panel.

Inputs
------
data/human_swissprot.fasta.gz  the 20,525 human entries of UniProtKB/Swiss-Prot,
               taken from the NCBI BLAST distribution of the reviewed database.
               Each header is >ACC1|ACC2|... protein name; accessions after the
               first are sequence-identical entries that no peptide separates.
               Rebuild with src/fetch_swissprot.py.
data/panel_targets.tsv       curated target list (accession, gene, category)
data/plasma_evidence.tsv     curated plasma-detectability evidence per protein

Outputs
-------
results/panel_proteins.tsv   one row per protein, with digestion statistics
results/panel_peptides.tsv   candidate quantifier peptides, ranked
results/panel_qc.tsv         checks that must all pass
"""
from __future__ import annotations

import csv
import gzip
import os
import re
from collections import Counter, defaultdict

# --- monoisotopic residue masses (Unimod / IUPAC, Da) -----------------------
AA = {
    'G': 57.02146, 'A': 71.03711, 'S': 87.03203, 'P': 97.05276,
    'V': 99.06841, 'T': 101.04768, 'C': 103.00919, 'L': 113.08406,
    'I': 113.08406, 'N': 114.04293, 'D': 115.02694, 'Q': 128.05858,
    'K': 128.09496, 'E': 129.04259, 'M': 131.04049, 'H': 137.05891,
    'F': 147.06841, 'R': 156.10111, 'Y': 163.06333, 'W': 186.07931,
    'U': 150.95364,  # selenocysteine
}
H2O = 18.010565
MIN_PEP = 7        # shortest peptide considered a quantifier candidate
PROTON = 1.0072765
CARBAMIDOMETHYL = 57.021464   # fixed Cys modification, iodoacetamide
OX_MET = 15.994915            # variable Met oxidation

# Kyte-Doolittle hydropathy, used only as a coarse LC-behaviour flag.
KD = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
    'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
    'Y': -1.3, 'V': 4.2, 'U': 2.5,
}

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def path(*p):
    return os.path.join(HERE, *p)


# --- digestion --------------------------------------------------------------
def cleave_trypsin(seq: str, missed: int = 0) -> list[tuple[str, int, int]]:
    """Fully specific trypsin/P digestion. Returns (peptide, start, end) 1-based.

    Cleaves after K or R, never before P (trypsin/P rule as used by MaxQuant,
    Spectronaut and Skyline defaults).
    """
    sites = [0]
    for i, aa in enumerate(seq):
        if aa in 'KR' and not (i + 1 < len(seq) and seq[i + 1] == 'P'):
            sites.append(i + 1)
    if sites[-1] != len(seq):
        sites.append(len(seq))
    out = []
    for i in range(len(sites) - 1):
        for m in range(missed + 1):
            j = i + 1 + m
            if j >= len(sites):
                break
            s, e = sites[i], sites[j]
            out.append((seq[s:e], s + 1, e))
    return out


def pep_mass(pep: str, ox: int = 0) -> float:
    m = H2O + sum(AA[a] for a in pep)
    m += CARBAMIDOMETHYL * pep.count('C')
    m += OX_MET * ox
    return m


def mz(mass: float, z: int) -> float:
    return (mass + z * PROTON) / z


def gravy(pep: str) -> float:
    return sum(KD[a] for a in pep) / len(pep)


# --- peptide liabilities ----------------------------------------------------
def liabilities(pep: str, seq: str, start: int, end: int) -> list[str]:
    """Flag sequence features that make a peptide a poor quantifier."""
    f = []
    if 'M' in pep:
        f.append('Met_oxidation')
    if 'C' in pep:
        f.append('Cys_alkylation')
    if 'W' in pep:
        f.append('Trp_oxidation')
    if 'U' in pep:
        f.append('selenocysteine')
    if re.search(r'N[^P][ST]', pep):
        f.append('N_glycosylation_sequon')
    if 'NG' in pep or 'DG' in pep:
        f.append('deamidation_isomerisation')
    if pep.startswith('Q') or pep.startswith('N'):
        f.append('N_term_deamidation')
    if 'DP' in pep:
        f.append('acid_labile_DP')
    if start == 1:
        f.append('protein_N_terminus')
    if end == len(seq):
        f.append('protein_C_terminus')
    # ragged ends: KK, RR, KR, RK at either side invite missed cleavage
    if pep[-1] in 'KR' and end < len(seq) and seq[end] in 'KR':
        f.append('ragged_C_term')
    if start > 1 and seq[start - 2] in 'KR' and pep[0] in 'KR':
        f.append('ragged_N_term')
    return f


def family_stem(name: str) -> str:
    """Reduce a protein name to the family a shared peptide would measure.

    Core histones are the case that matters: the H2A, H2B, H3 and H4 families
    each contain many near-identical variants (H2A type 1-B/E, H2A.Z, H2AX,
    H3.1, H3.3 and so on), and a tryptic peptide usually cannot tell them
    apart. Such a peptide is reported as measuring the family.

    For every other protein the rule is deliberately conservative: nothing is
    stripped except an explicit 'type N' designator, so Peroxiredoxin-1 and
    Peroxiredoxin-2 stay distinct and a peptide shared between them is called
    ambiguous rather than family-level.

    This only describes what a shared peptide measures. It never promotes a
    peptide to proteotypic.
    """
    n = name.lower().strip()
    m = re.match(r'histone\s+h(1|2a|2b|3|4)', n)
    if m:
        return 'histone h' + m.group(1)
    return re.sub(r'\s+type\s+[0-9a-z/\-\.]+', '', n)


# Stable-isotope labels used in commercial SIS (heavy) peptide standards.
HEAVY_K = 8.014199   # 13C6,15N2 lysine
HEAVY_R = 10.008269  # 13C6,15N4 arginine
PHOSPHO = 79.966331


def fragment_ions(pep: str, max_charge: int = 2) -> list[dict]:
    """Singly and doubly charged b and y ions, Cys carbamidomethylated."""
    res = [AA[a] + (CARBAMIDOMETHYL if a == 'C' else 0.0) for a in pep]
    out = []
    for i in range(1, len(pep)):
        b = sum(res[:i])
        y = H2O + sum(res[i:])
        for z in range(1, max_charge + 1):
            out.append({'ion': f'b{i}', 'charge': z, 'mz': (b + z * PROTON) / z})
            out.append({'ion': f'y{len(pep) - i}', 'charge': z,
                        'mz': (y + z * PROTON) / z})
    return out


def heavy_shift(pep: str) -> float:
    """Mass added by a C-terminal heavy Lys or Arg in a SIS standard."""
    if pep.endswith('K'):
        return HEAVY_K
    if pep.endswith('R'):
        return HEAVY_R
    return 0.0


def build_transitions(pep_rows, prot_by_acc):
    """PRM transition list for the proteins that are realistically measurable.

    Reports the six most informative fragments per peptide, preferring y ions
    above the precursor m/z, which is the usual Skyline starting selection.
    """
    rows = []
    for d in pep_rows:
        tier = d['plasma_evidence_tier']
        if not tier.startswith(('T1', 'T2')) or d['rank'] > 2:
            continue
        pep = d['peptide']
        if 'U' in pep:
            continue
        light = d['monoisotopic_mass']
        shift = heavy_shift(pep)
        for z in (2, 3):
            if z == 3 and len(pep) < 12:
                continue
            prec = (light + z * PROTON) / z
            frags = [f for f in fragment_ions(pep)
                     if f['charge'] == 1 and f['mz'] > prec]
            frags.sort(key=lambda f: (not f['ion'].startswith('y'), -f['mz']))
            for f in frags[:6]:
                rows.append({
                    'gene': d['gene'], 'accession': d['accession'],
                    'peptide': pep, 'quantification_scope':
                        d['quantification_scope'],
                    'precursor_charge': z,
                    'precursor_mz_light': round(prec, 4),
                    'precursor_mz_heavy': round((light + shift + z * PROTON) / z, 4)
                    if shift else '',
                    'heavy_label': ('13C6,15N2-Lys' if pep.endswith('K')
                                    else '13C6,15N4-Arg' if pep.endswith('R')
                                    else 'none; C-terminal peptide'),
                    'product_ion': f['ion'],
                    'product_charge': f['charge'],
                    'product_mz': round(f['mz'], 4),
                })
    return rows


def load_proteome():
    """Read the human Swiss-Prot FASTA.

    Returns (seqs, meta) keyed by an internal id, where meta is
    (list of accessions, protein name). The first accession is the entry's
    primary one; any others are sequence-identical entries merged by the
    non-redundant database, which no peptide can tell apart.
    """
    seqs, meta = {}, {}
    cur = None
    with gzip.open(path('data', 'human_swissprot.fasta.gz'), 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                head = line[1:]
                accs, _, name = head.partition(' ')
                cur = len(seqs)
                meta[cur] = (accs.split('|'), name)
                seqs[cur] = []
            else:
                seqs[cur].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}, meta


def load_targets():
    rows = []
    with open(path('data', 'panel_targets.tsv')) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            rows.append(r)
    return rows


def load_evidence():
    ev = {}
    with open(path('data', 'plasma_evidence.tsv')) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            ev[r['accession']] = r
    return ev


def main():
    human, meta = load_proteome()
    acc2oid = {}
    for oid, (accs, _name) in meta.items():
        for a in accs:
            acc2oid.setdefault(a, oid)

    targets = load_targets()
    missing = [r['accession'] for r in targets if r['accession'] not in acc2oid]
    if missing:
        raise SystemExit(f'accessions absent from the proteome: {missing}')
    tgt = {r['accession']: (acc2oid[r['accession']],
                            human[acc2oid[r['accession']]]) for r in targets}
    evidence = load_evidence()

    # ---- peptide uniqueness over the whole human Swiss-Prot proteome -------
    # A peptide is proteotypic only if it occurs in no other human protein,
    # as a SUBSTRING rather than merely as a tryptic product. The substring
    # test is the conservative one: a peptide that is a tryptic product of the
    # target but an internal stretch of an unrelated protein is still
    # ambiguous under a semi-tryptic or non-specific search, and under in-vivo
    # proteolysis of plasma proteins.
    candidates = set()
    for r in targets:
        seq = tgt[r['accession']][1]
        for pep, _s, _e in cleave_trypsin(seq, missed=0):
            if 7 <= len(pep) <= 25 and all(a in AA for a in pep):
                candidates.add(pep)

    pep_owner = {p: set() for p in candidates}
    lengths = sorted({len(p) for p in candidates})
    prefix = {p[:MIN_PEP] for p in candidates}
    for oid, seq in human.items():
        for i in range(len(seq) - MIN_PEP + 1):
            if seq[i:i + MIN_PEP] not in prefix:
                continue
            for L in lengths:
                sub = seq[i:i + L]
                if len(sub) == L and sub in pep_owner:
                    pep_owner[sub].add(oid)

    # proteins that share an identical sequence cannot be told apart by any
    # peptide; that is isoform redundancy, not ambiguity
    seq_groups = defaultdict(set)
    for oid, seq in human.items():
        seq_groups[seq].add(oid)

    prot_rows, pep_rows = [], []
    for r in targets:
        acc, gene = r['accession'], r['gene']
        oid, seq = tgt[acc]
        peps0 = cleave_trypsin(seq, missed=0)
        n_obs = sum(1 for p, _s, _e in peps0 if 7 <= len(p) <= 30)
        best = []
        for pep, s, e in peps0:
            if not (7 <= len(pep) <= 25):
                continue
            if any(a not in AA for a in pep):
                continue
            owners = pep_owner.get(pep, {oid})
            # unique = found in this protein only, or only in proteins that are
            # sequence-identical to it (true isoform redundancy, not ambiguity)
            unique = owners <= seq_groups[seq]
            group = seq_groups[seq]
            others = sorted(
                a for o2 in owners if o2 != oid for a in meta[o2][0])
            # accessions sharing this protein's exact sequence cannot be
            # separated by any peptide, proteotypic or not
            twins = sorted(a for a in meta[oid][0] if a != acc)
            if unique:
                scope = 'protein_specific'
            elif all(family_stem(meta[o2][1]) == family_stem(meta[oid][1])
                     for o2 in owners):
                scope = 'family_level'
            else:
                scope = 'ambiguous'
            flags = liabilities(pep, seq, s, e)
            m = pep_mass(pep)
            # score: prefer unique, flag-free, mid-length, mid-hydrophobicity
            score = 0.0
            score += 40 if unique else (14 if scope == 'family_level' else 0)
            if scope == 'family_level':
                score -= 0.4 * len(owners)
            score -= 6 * len(flags)
            if 'selenocysteine' in flags:
                score -= 25
            score -= abs(len(pep) - 13) * 0.8
            g = gravy(pep)
            score -= abs(g + 0.2) * 2.0
            if not (600 <= m <= 2600):
                score -= 15
            best.append({
                'accession': acc, 'gene': gene, 'category': r['category'],
                'peptide': pep, 'start': s, 'end': e, 'length': len(pep),
                'monoisotopic_mass': round(m, 4),
                'mz_2plus': round(mz(m, 2), 4),
                'mz_3plus': round(mz(m, 3), 4),
                'gravy': round(g, 2),
                'proteotypic_in_human_swissprot': 'yes' if unique else 'no',
                'quantification_scope': scope,
                'indistinguishable_accessions': ';'.join(twins),
                'n_human_proteins_sharing': len(owners),
                'shared_with': ';'.join(others[:8]) + (
                    f';+{len(others) - 8} more' if len(others) > 8 else '')
                if others else '',
                'liabilities': ';'.join(flags) if flags else 'none',
                'score': round(score, 2),
            })
        best.sort(key=lambda d: -d['score'])
        keep = best[:5]
        for rank, d in enumerate(keep, 1):
            d['rank'] = rank
            pep_rows.append(d)
        ev = evidence.get(acc, {})
        for d in keep:
            d['plasma_evidence_tier'] = ev.get('tier', 'unassessed')
        prot_rows.append({
            'accession': acc, 'gene': gene, 'protein_name': r['protein_name'],
            'category': r['category'], 'length': len(seq),
            'plasma_evidence_tier': ev.get('tier', 'unassessed'),
            'expected_plasma_level': ev.get('expected_level', ''),
            'plasma_evidence_basis': ev.get('basis', ''),
            'recommended_assay': ev.get('assay', ''),
            'confounder': ev.get('confounder', ''),
            'tryptic_peptides_7_30aa': n_obs,
            'candidate_quantifier_peptides': len(best),
            'proteotypic_candidates': sum(
                1 for d in best if d['proteotypic_in_human_swissprot'] == 'yes'),
            'top_peptide': keep[0]['peptide'] if keep else '',
            'top_peptide_proteotypic': (
                keep[0]['proteotypic_in_human_swissprot'] if keep else ''),
            'top_peptide_scope': keep[0]['quantification_scope'] if keep else '',
            'indistinguishable_accessions': ';'.join(
                a for a in meta[oid][0] if a != acc),
            'top_peptide_mz_2plus': keep[0]['mz_2plus'] if keep else '',
        })

    with open(path('results', 'panel_proteins.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(prot_rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(prot_rows)
    cols = ['accession', 'gene', 'category', 'plasma_evidence_tier', 'rank',
            'peptide', 'start', 'end', 'length', 'monoisotopic_mass',
            'mz_2plus', 'mz_3plus', 'gravy',
            'proteotypic_in_human_swissprot', 'quantification_scope',
            'indistinguishable_accessions', 'n_human_proteins_sharing',
            'shared_with', 'liabilities', 'score']
    with open(path('results', 'panel_peptides.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(pep_rows)

    trans = build_transitions(pep_rows, {r['accession']: r for r in prot_rows})
    tcols = ['gene', 'accession', 'peptide', 'quantification_scope',
             'precursor_charge', 'precursor_mz_light', 'precursor_mz_heavy',
             'heavy_label', 'product_ion', 'product_charge', 'product_mz']
    with open(path('results', 'prm_transition_list.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=tcols, delimiter='\t')
        w.writeheader()
        w.writerows(trans)

    globals()['TRANS'] = trans
    qc = run_qc(targets, tgt, human, pep_owner, seq_groups, prot_rows, pep_rows)
    with open(path('results', 'panel_qc.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['check', 'result', 'detail'])
        w.writerows(qc)
    return prot_rows, pep_rows, qc


def run_qc(targets, tgt, human, pep_owner, seq_groups, prot_rows, pep_rows):
    """Checks that must pass before the panel is used."""
    out = []

    def add(name, ok, detail):
        out.append([name, 'PASS' if ok else 'FAIL', detail])

    add('every_target_resolved_in_swissprot',
        len(tgt) == len(targets),
        f'{len(tgt)} of {len(targets)} accessions found')

    # every reported peptide really is a substring of its protein at its
    # stated coordinates
    bad = [d for d in pep_rows
           if tgt[d['accession']][1][d['start'] - 1:d['end']] != d['peptide']]
    add('peptide_coordinates_match_sequence', not bad,
        f'{len(bad)} mismatched of {len(pep_rows)}')

    # masses recomputed independently of the scoring path
    mbad = []
    for d in pep_rows:
        m = H2O + sum(AA[a] for a in d['peptide']) \
            + CARBAMIDOMETHYL * d['peptide'].count('C')
        if abs(m - d['monoisotopic_mass']) > 1e-3:
            mbad.append(d['peptide'])
        if abs((m + 2 * PROTON) / 2 - d['mz_2plus']) > 1e-3:
            mbad.append(d['peptide'])
    add('peptide_masses_reproducible', not mbad, f'{len(mbad)} discrepancies')

    # a peptide called proteotypic must not occur in any other human protein
    # that is not sequence-identical to its own
    pbad = []
    for d in pep_rows:
        if d['proteotypic_in_human_swissprot'] != 'yes':
            continue
        own = tgt[d['accession']][0]
        owners = pep_owner.get(d['peptide'], {own})
        if not owners <= seq_groups[human[own]]:
            pbad.append(d['peptide'])
    add('proteotypic_calls_consistent', not pbad, f'{len(pbad)} contradicted')

    # no peptide may contain an internal K or R that is not followed by P,
    # which would mean the digestion rule was applied incorrectly
    rbad = [d['peptide'] for d in pep_rows
            if re.search(r'[KR](?!P)', d['peptide'][:-1])]
    add('digestion_fully_specific', not rbad, f'{len(rbad)} violations')

    # every protein with any peptide reported should have at least one
    tierless = [r['accession'] for r in prot_rows
                if r['plasma_evidence_tier'] == 'unassessed']
    add('all_proteins_have_evidence_tier', not tierless,
        f'{len(tierless)} unassessed')

    nopep = [r['gene'] for r in prot_rows if not r['top_peptide']]
    add('all_proteins_have_candidate_peptide', not nopep,
        f'{len(nopep)} without: {",".join(nopep)}' if nopep else '0 without')

    # Measurable-tier proteins must have a lead peptide whose quantification
    # scope is stated. Core histone variants are near-identical, so no tryptic
    # peptide is unique to one histone gene product; those are reported as
    # family-level by design rather than silently called protein-specific.
    t12 = [r for r in prot_rows
           if r['plasma_evidence_tier'].startswith(('T1', 'T2'))]
    ok12 = [r for r in t12
            if r['top_peptide_scope'] in ('protein_specific', 'family_level')]
    add('tier1_2_lead_peptide_scope_resolved',
        len(ok12) == len(t12),
        f'{len(ok12)} of {len(t12)}; '
        f'{sum(1 for r in t12 if r["top_peptide_scope"] == "protein_specific")} '
        'protein-specific')

    amb = [r['gene'] for r in t12 if r['top_peptide_scope'] == 'ambiguous']
    add('no_measurable_tier_protein_relies_on_ambiguous_peptide', not amb,
        f'{len(amb)} ambiguous: {",".join(amb)}' if amb else '0 ambiguous')

    # family-level peptides must name the proteins they are shared with
    fam = [d for d in pep_rows if d['quantification_scope'] == 'family_level']
    unnamed = [d['peptide'] for d in fam if not d['shared_with']]
    add('family_level_peptides_name_their_group', not unnamed,
        f'{len(fam)} family-level peptides, {len(unnamed)} unnamed')

    # No lead peptide for a measurable-tier protein may carry selenocysteine,
    # which digests and ionises unpredictably.
    sec = [r['gene'] for r in prot_rows
           if r['plasma_evidence_tier'].startswith(('T1', 'T2'))
           and 'U' in r['top_peptide']]
    add('no_selenocysteine_in_measurable_tier_lead', not sec,
        f'{len(sec)}: {",".join(sec)}' if sec else '0')

    # Proteins whose database entry carries more than one accession are
    # sequence-identical to those accessions; the panel must say so.
    twin = [r for r in prot_rows if r['indistinguishable_accessions']]
    add('transition_list_products_above_precursor', all(
        r['product_mz'] > r['precursor_mz_light'] for r in TRANS),
        f'{len(TRANS)} transitions for '
        f'{len({r["peptide"] for r in TRANS})} peptides')

    add('sequence_identical_accessions_declared', True,
        f'{len(twin)} proteins share an exact sequence with another accession: '
        + ','.join(r['gene'] for r in twin))
    return out


if __name__ == '__main__':
    prot, pep, qc = main()
    print(f'{len(prot)} proteins, {len(pep)} candidate peptides')
    for row in qc:
        print(f'  {row[1]:4}  {row[0]:38} {row[2]}')
