"""Independent validation of a panel build.

    python3 src/validate.py lymphocyte     (default)
    python3 src/validate.py plasma

Re-derives every reported quantity by a second route and compares:
  * peptide masses recomputed with pyteomics
  * tryptic digestion compared against the ExPASy enzymology rule
  * peptide coordinates re-read from the sequence database
  * proteotypicity re-tested by brute-force substring search

Run after src/panel_build.py for the same matrix. Exits non-zero on any
disagreement.
"""
from __future__ import annotations

import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panel_build import AA, CARBAMIDOMETHYL, MATRICES, PROTON, path

try:
    from pyteomics import mass as pyt_mass
except ImportError:
    pyt_mass = None


def expasy_sites(s: str) -> set[int]:
    """ExPASy PeptideCutter trypsin rule, including the WKP / MRP exceptions."""
    out = set()
    for i, a in enumerate(s):
        nxt = s[i + 1] if i + 1 < len(s) else ''
        prv = s[i - 1] if i else ''
        if a in 'KR' and nxt != 'P':
            out.add(i + 1)
        elif a == 'K' and nxt == 'P' and prv == 'W':
            out.add(i + 1)
        elif a == 'R' and nxt == 'P' and prv == 'M':
            out.add(i + 1)
    return out


def panel_sites(s: str) -> set[int]:
    """The rule used by the panel: cleave after K/R, never before P."""
    return {i + 1 for i, a in enumerate(s)
            if a in 'KR' and not (i + 1 < len(s) and s[i + 1] == 'P')}


def main(matrix: str = 'lymphocyte') -> int:
    from panel_build import load_proteome, load_targets
    human, meta = load_proteome()
    acc2oid = {}
    for oid, (accs, _n) in meta.items():
        for a in accs:
            acc2oid.setdefault(a, oid)
    tgt = {r['accession']: (acc2oid[r['accession']],
                            human[acc2oid[r['accession']]])
           for r in load_targets()}
    with open(path('results', f'{matrix}_panel_peptides.tsv')) as fh:
        peps = list(csv.DictReader(fh, delimiter='\t'))

    fails = []

    # 1. coordinates
    n = 0
    for r in peps:
        seq = tgt[r['accession']][1]
        if seq[int(r['start']) - 1:int(r['end'])] != r['peptide']:
            fails.append(f"coordinate mismatch {r['accession']} {r['peptide']}")
        n += 1
    print(f'coordinates re-read from database      {n - len(fails)}/{n} match')

    # 2. masses, recomputed with pyteomics where available
    if pyt_mass is not None:
        aa = dict(pyt_mass.std_aa_mass)
        aa['U'] = AA['U']
        bad = 0
        for r in peps:
            m = (pyt_mass.fast_mass(r['peptide'], aa_mass=aa)
                 + CARBAMIDOMETHYL * r['peptide'].count('C'))
            if abs(m - float(r['monoisotopic_mass'])) > 0.002:
                bad += 1
            z2 = (m + 2 * PROTON) / 2
            if abs(z2 - float(r['mz_2plus'])) > 0.002:
                bad += 1
        print(f'masses recomputed with pyteomics       '
              f'{len(peps) - bad}/{len(peps)} agree within 2 mDa')
        if bad:
            fails.append(f'{bad} mass disagreements')
    else:
        print('masses recomputed with pyteomics       SKIPPED (not installed)')

    # 3. digestion rule: differences from ExPASy must be only the documented
    #    W-K-P and M-R-P exceptions, which search engines do not cleave.
    ctx, ndiff = set(), 0
    for acc, (_oid, s) in tgt.items():
        d = expasy_sites(s) ^ panel_sites(s)
        if d:
            ndiff += 1
            ctx |= {s[p - 2:p + 1] for p in d}
    ok = ctx <= {'WKP', 'MRP'}
    print(f'digestion vs ExPASy rule               '
          f'{ndiff} proteins differ, only at {sorted(ctx)}, '
          f'{"expected" if ok else "UNEXPECTED"}')
    if not ok:
        fails.append(f'unexpected digestion differences: {sorted(ctx)}')

    # 4. proteotypicity by brute-force substring search over the proteome
    checked = bad = 0
    for r in peps:
        if r['proteotypic_in_human_swissprot'] != 'yes':
            continue
        own_seq = tgt[r['accession']][1]
        owners = [o for o, s in human.items() if r['peptide'] in s]
        if any(human[o] != own_seq for o in owners):
            bad += 1
        checked += 1
    print(f'proteotypicity by brute-force search   '
          f'{checked - bad}/{checked} confirmed')
    if bad:
        fails.append(f'{bad} peptides called proteotypic are not')

    # 5. fragment ions and heavy-label shifts in the transition list
    tpath = path('results', f'{matrix}_prm_transitions.tsv')
    if pyt_mass is not None and os.path.exists(tpath):
        with open(tpath) as fh:
            trans = list(csv.DictReader(fh, delimiter='\t'))
        aa2 = dict(pyt_mass.std_aa_mass)
        aa2['C'] = aa2['C'] + CARBAMIDOMETHYL
        fbad = hbad = 0
        for r in trans:
            pep, ion = r['peptide'], r['product_ion']
            n, kind = int(ion[1:]), ion[0]
            sub = pep[:n] if kind == 'b' else pep[-n:]
            m = pyt_mass.fast_mass(sub, ion_type=kind, charge=1, aa_mass=aa2)
            if abs(m - float(r['product_mz'])) > 0.003:
                fbad += 1
            if r['precursor_mz_heavy']:
                z = int(r['precursor_charge'])
                d = (float(r['precursor_mz_heavy'])
                     - float(r['precursor_mz_light'])) * z
                exp = 8.014199 if pep.endswith('K') else 10.008269
                if abs(d - exp) > 0.002:
                    hbad += 1
        print(f'product ions recomputed with pyteomics '
              f'{len(trans) - fbad}/{len(trans)} agree within 3 mDa')
        print(f'heavy SIS precursor shifts             '
              f'{hbad} errors')
        if fbad:
            fails.append(f'{fbad} product ion disagreements')
        if hbad:
            fails.append(f'{hbad} heavy-label shift errors')

    if fails:
        print('\nFAILED:')
        for f in fails:
            print('  ' + f)
        return 1
    print('\nall independent checks agree with the panel build')
    return 0


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'lymphocyte'
    if which not in MATRICES:
        raise SystemExit(f'matrix must be one of {sorted(MATRICES)}')
    print(f'validating: {which}')
    raise SystemExit(main(which))
