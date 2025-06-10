#!/usr/bin/env python3
"""
tailclip_filter.py v0.9 Rhys Parry r.parry@uq.edu.au

Separates reads based on presence of any non-template end additions
indicated by 3'-end soft or hard clipping in the CIGAR (use local
alignment). In addition, tallies raw counts of all single-nucleotide
and dinucleotide combinations present in the trailing soft-clipped
sequence at the 3' end.

Usage:
    tailclip_filter.py <input.bam> <noClip.bam> <withClip.bam>

  - <noClip.bam>   : mapped reads **without** trailing soft/hard clips
  - <withClip.bam> : mapped reads **with** 3'-end soft or hard clipping

Outputs summary (to stderr):
  • Total mapped reads processed
  • No clip reads
  • Clipped reads
  • Single-nucleotide counts (A, C, G, T) in the last clipped base
  • Dinucleotide counts (AA, AC, ..., TT) from the last two clipped bases

Notes:
  • Unmapped reads are skipped for speed and not written to any BAM.
  • Must map with Bowtie2 `--local` so soft-clipped bases appear in SEQ.
  • Hard-clipped bases (CIGAR 5) are removed from SEQ and cannot be counted.
  • If soft-clip length <1, that read is not included in clip-base counts.
  • If soft-clip length ==1, only single-base counts are incremented.
  • If soft-clip length >=2, both single- and di-nucleotide counts are incremented.

Dependencies:
  pysam
"""
import sys
from collections import Counter
import pysam

def usage():
    sys.stderr.write("Usage: filter_3end_clipping.py <input.bam> <noClip.bam> <withClip.bam>\n")
    sys.exit(1)


def process(input_bam, no_clip_bam, with_clip_bam):
    sam = pysam.AlignmentFile(input_bam, "rb")
    no_clip = pysam.AlignmentFile(no_clip_bam,   "wb", template=sam)
    with_clip = pysam.AlignmentFile(with_clip_bam, "wb", template=sam)

    total = 0
    no_clip_count = 0
    with_clip_count = 0
    single_counts = Counter()
    dinuc_counts = Counter()
    bases = ['A','C','G','T']

    for read in sam.fetch(until_eof=True):
        # skip unmapped reads entirely
        if read.is_unmapped:
            continue

        total += 1
        cig = read.cigartuples or []
        # soft-clip = 4, hard-clip = 5
        if cig and cig[-1][0] in (4,5):
            with_clip.write(read)
            with_clip_count += 1
            # count only soft-clips for base content
            if cig[-1][0] == 4:
                clip_len = cig[-1][1]
                seq = read.query_sequence or ''
                if clip_len >= 1 and len(seq) >= clip_len:
                    tail = seq[-clip_len:].upper()
                    # single-base count
                    last_base = tail[-1]
                    if last_base in bases:
                        single_counts[last_base] += 1
                    # dinucleotide count
                    if clip_len >= 2:
                        dinuc = tail[-2:]
                        if all(b in bases for b in dinuc):
                            dinuc_counts[dinuc] += 1
        else:
            no_clip.write(read)
            no_clip_count += 1

    sam.close()
    no_clip.close()
    with_clip.close()

    # Summary
    sys.stderr.write(f"Total mapped reads: {total}\n")
    sys.stderr.write(f"No clip reads     : {no_clip_count}\n")
    sys.stderr.write(f"Clipped reads     : {with_clip_count}\n")
    sys.stderr.write("\nSingle-base counts (soft-clip tails):\n")
    for b in bases:
        sys.stderr.write(f"{b}: {single_counts.get(b,0)}\n")
    sys.stderr.write("\nDinucleotide counts (soft-clip tails):\n")
    for b1 in bases:
        for b2 in bases:
            key = b1 + b2
            sys.stderr.write(f"{key}: {dinuc_counts.get(key,0)}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage()
    input_bam, no_clip_bam, with_clip_bam = sys.argv[1:]
    process(input_bam, no_clip_bam, with_clip_bam)
