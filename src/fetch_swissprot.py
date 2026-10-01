"""Rebuild data/human_swissprot.fasta.gz from the NCBI BLAST distribution.

UniProtKB/Swiss-Prot is distributed by NCBI as a BLAST database in the public
bucket s3://ncbi-blast-databases. This script downloads that database, decodes
the index and header files, and writes out the human entries.

The UniProt REST API was not reachable from the environment this panel was
built in, which is why the sequences come from this mirror rather than from
uniprot.org directly. The content is the same reviewed database.

Usage:  python3 src/fetch_swissprot.py [workdir]
"""
from __future__ import annotations

import bisect
import gzip
import os
import re
import struct
import sys
import urllib.request

BUCKET = 'https://ncbi-blast-databases.s3.amazonaws.com'
EXTS = ('pin', 'phr', 'psq')
# BLAST version 5 packs residues as indices into this alphabet
ALPHABET = '-ABCDEFGHIKLMNPQRSTVWXYZU*OJ'
ACC = re.compile(rb'\x1a\x06([A-Z][0-9][A-Z0-9]{3}[0-9])'
                 rb'|\x1a\x0a([A-Z][0-9][A-Z0-9]{3}[0-9][A-Z0-9]{4})')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def download(workdir: str) -> str:
    os.makedirs(workdir, exist_ok=True)
    with urllib.request.urlopen(f'{BUCKET}/latest-dir', timeout=60) as r:
        release = r.read().decode().strip()
    print(f'release {release}')
    for ext in EXTS:
        dest = os.path.join(workdir, f'swissprot.{ext}')
        if os.path.exists(dest):
            continue
        url = f'{BUCKET}/{release}/swissprot.{ext}'
        print(f'  downloading {url}')
        urllib.request.urlretrieve(url, dest)
    return release


def read_index(pin: bytes):
    p = 0

    def i32():
        nonlocal p
        v = struct.unpack('>i', pin[p:p + 4])[0]
        p += 4
        return v

    def pstr():
        nonlocal p
        n = i32()
        v = pin[p:p + n]
        p += n
        return v.decode()

    ver, dbtype = i32(), i32()
    if (ver, dbtype) != (5, 1):
        raise SystemExit(f'unexpected BLAST db version/type {ver}/{dbtype}')
    i32()                      # volume
    title, _lmdb, date = pstr(), pstr(), pstr()
    n = i32()
    p += 8                     # total residues
    i32()                      # longest sequence
    hoff = struct.unpack(f'>{n + 1}i', pin[p:p + 4 * (n + 1)])
    p += 4 * (n + 1)
    soff = struct.unpack(f'>{n + 1}i', pin[p:p + 4 * (n + 1)])
    return title, date, n, hoff, soff


def main(workdir: str) -> int:
    release = download(workdir)
    pin = open(os.path.join(workdir, 'swissprot.pin'), 'rb').read()
    phr = open(os.path.join(workdir, 'swissprot.phr'), 'rb').read()
    psq = open(os.path.join(workdir, 'swissprot.psq'), 'rb').read()
    title, date, n, hoff, soff = read_index(pin)
    print(f'{title}, {date}, {n} entries')

    human = sorted({bisect.bisect_right(hoff, m.start()) - 1
                    for m in re.finditer(rb'\[Homo sapiens\]', phr)})
    print(f'{len(human)} human entries')

    out = os.path.join(HERE, 'data', 'human_swissprot.fasta.gz')
    with gzip.open(out, 'wt', compresslevel=9) as fh:
        for oid in human:
            head = phr[hoff[oid]:hoff[oid + 1]]
            accs = []
            for m in ACC.finditer(head):
                a = (m.group(1) or m.group(2)).decode()
                if a not in accs:
                    accs.append(a)
            texts = [t.decode('latin1')
                     for t in re.findall(rb'[\x20-\x7e]{12,}', head)]
            hs = [t for t in texts if 'Homo sapiens' in t]
            t = hs[0] if hs else (texts[0] if texts else '')
            m = re.search(r'RecName: Full=([^;\[]+)', t)
            name = m.group(1).strip() if m else t[:60]
            seq = ''.join(ALPHABET[b]
                          for b in psq[soff[oid]:soff[oid + 1] - 1])
            fh.write(f">{'|'.join(accs)} {name}\n{seq}\n")
    print(f'wrote {out} from release {release}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else '/tmp/swissprot'))
