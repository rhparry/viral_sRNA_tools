#!/usr/bin/env python3
"""
v0.1  Rhys Parry <r.parry@uq.edu.au>

Takes a FASTQ (gzipped or not), computes:
  • histogram of read lengths
  • counts of first-base bias and last-base bias per length
  • percentage of total reads for each of those counts

Outputs a TSV:
  Length  Count  A_first  C_first  G_first  T_first  A_last  C_last  G_last  T_last
         A_first_% ... T_first_%  A_last_% ... T_last_%
"""

import sys
import gzip
import argparse
from collections import Counter, defaultdict

def open_fastq(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    else:
        return open(path, 'r')

def main():
    p = argparse.ArgumentParser(description="FASTQ first+last-base bias histogram")
    p.add_argument('fastq', help="Input FASTQ file (can be .gz)")
    args = p.parse_args()

    # Counters
    hist = Counter()
    first = {b: Counter() for b in 'ACGT'}
    last  = {b: Counter() for b in 'ACGT'}
    total_reads = 0

    with open_fastq(args.fastq) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip('\n')
            # skip '+' line and quals
            fh.readline()
            fh.readline()

            L = len(seq)
            if L == 0:
                continue

            total_reads += 1
            hist[L] += 1

            f = seq[0].upper()
            l = seq[-1].upper()
            if f in first:
                first[f][L] += 1
            if l in last:
                last[l][L] += 1

    # prepare header
    bases = ['A','C','G','T']
    hdr = ['Length','Count'] \
        + [f"{b}_first" for b in bases] \
        + [f"{b}_last"  for b in bases] \
        + [f"{b}_first_pct" for b in bases] \
        + [f"{b}_last_pct"  for b in bases]
    print('\t'.join(hdr))

    # output rows sorted by length
    for L in sorted(hist):
        cnt = hist[L]
        row = [str(L), str(cnt)]
        # raw counts
        for b in bases:
            row.append(str(first[b].get(L, 0)))
        for b in bases:
            row.append(str(last[b].get(L, 0)))
        # percentages (of total reads)
        for b in bases:
            pct = first[b].get(L, 0) / total_reads * 100
            row.append(f"{pct:.3f}")
        for b in bases:
            pct = last[b].get(L, 0) / total_reads * 100
            row.append(f"{pct:.3f}")

        print('\t'.join(row))


if __name__ == '__main__':
    main()
