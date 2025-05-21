#!/usr/bin/env python3
"""
filter_3end_clipping.py v1.5 Rhys Parry <r.parry@uq.edu.au>

Processes an end-to-end or local-aligned BAM to separate all mapped reads
(bypass any length filtering) based on presence of non-template 3'-end
additions. Detects soft-clips (local mode) or terminal mismatches (E2E via MD tag),
then tallies single-nucleotide and dinucleotide motifs of the tail.

Usage:
    filter_3end_clipping.py <input.bam> <noTail.bam> <withTail.bam>

  - <noTail.bam>   : mapped reads without 3'-end additions
  - <withTail.bam> : mapped reads with 3'-end soft-clips or mismatches

Outputs summary (to stderr):
  • Total processed reads
  • No-tail reads count
  • Tail reads count
  • Single-base counts of tail nucleotides
  • Dinucleotide counts of last two tail nucleotides

Notes:
  • Unmapped reads are skipped; every mapped read is examined.
  • Soft-clips (CIGAR=4) yield tail sequence directly from SEQ.
  • Hard-clips (CIGAR=5) remove bases — fallback to MD mismatches.
  • End-to-end mismatches: MD tag encodes mismatches, trailing matches;
    single-base mismatch at end appears as "<num><base>0";
    two-base mismatches at end appear as "0<base1>0<base2>0".

Dependencies:
  pysam
"""
import sys
from collections import Counter
import pysam
import re

def usage():
    sys.stderr.write("Usage: filter_3end_clipping.py <input.bam> <noTail.bam> <withTail.bam>\n")
    sys.exit(1)


def process(in_bam, out_no, out_yes):
    sam = pysam.AlignmentFile(in_bam, "rb")
    no_tail = pysam.AlignmentFile(out_no,  "wb", template=sam)
    with_tail = pysam.AlignmentFile(out_yes, "wb", template=sam)

    total = no_count = yes_count = 0
    single = Counter()
    di = Counter()
    bases = set('ACGT')
    # regex for two mismatches at end: 0B10C0 (trailing 0)
    re_two = re.compile(r"0([ACGT])0([ACGT])0$")
    # regex for single mismatch at end: <num><base>0
    re_one = re.compile(r"\d+([ACGT])0$")

    for read in sam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        total += 1
        seq = read.query_sequence or ''

        # extract tail: check CIGAR first
        cig = read.cigartuples or []
        tail = ''
        if cig and cig[-1][0] == 4:
            # soft-clip: tail bases available
            clip_len = cig[-1][1]
            tail = seq[-clip_len:]
        elif cig and cig[-1][0] == 5:
            # hard-clip: no SEQ for tail
            tail = ''
        else:
            # no clip: check MD for mismatches at end
            try:
                md = read.get_tag('MD')
                # try two-base mismatch first
                m2 = re_two.search(md)
                m1 = None
                if not m2:
                    # then single-base
                    m1 = re_one.search(md)
                if m2:
                    tail = m2.group(1) + m2.group(2)
                elif m1:
                    tail = m1.group(1)
            except KeyError:
                tail = ''

        if tail:
            with_tail.write(read)
            yes_count += 1
            t = tail.upper()
            # single-base tally
            if t[-1] in bases:
                single[t[-1]] += 1
            # dinucleotide tally
            if len(t) >= 2 and all(b in bases for b in t[-2:]):
                di[t[-2:]] += 1
        else:
            no_tail.write(read)
            no_count += 1

    sam.close()
    no_tail.close()
    with_tail.close()

    # summary
    sys.stderr.write(f"Total reads processed: {total}\n")
    sys.stderr.write(f"No-tail reads        : {no_count}\n")
    sys.stderr.write(f"Tail reads           : {yes_count}\n")
    sys.stderr.write("\nSingle-base counts:\n")
    for b in 'ACGT':
        sys.stderr.write(f"{b}: {single.get(b,0)}\n")
    sys.stderr.write("\nDinucleotide counts:\n")
    for b1 in 'ACGT':
        for b2 in 'ACGT':
            key = b1 + b2
            sys.stderr.write(f"{key}: {di.get(key,0)}\n")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        usage()
    process(*sys.argv[1:4])
