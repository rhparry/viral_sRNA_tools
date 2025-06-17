#!/usr/bin/env python3
"""
tailclip_filter.py v2.5 Rhys Parry r.parry@uq.edu.au

Separates reads based on 3'-end soft- or hard-clipping and classifies non-template additions (NTAs)
as mono-, di-, tri-, tetra-, or longer nucleotide tails. Uses reference genome to determine if the clipped bases
match the genome or are likely non-template additions.

USAGE:
    tailclip_filter.py <input.bam> <noClip.bam> <withClip.bam> <reference.fa> [output_prefix]

INPUTS:
    - <input.bam>:        BAM file mapped using local alignment (e.g., Bowtie2 --local)
    - <noClip.bam>:       Output BAM for reads without 3'-end soft/hard clipping
    - <withClip.bam>:     Output BAM for reads with 3'-end soft/hard clipping
    - <reference.fa>:     Reference genome in FASTA format (must be indexed with `samtools faidx`)
    - [output_prefix]:    Optional output prefix for the TSV file and plots (default: withClip.bam)

REQUIREMENTS:
    - BAM file must be aligned to the same reference FASTA file provided as input
    - Reference FASTA **must be indexed** (i.e., `<reference.fa>.fai` must exist)
    - For multi-sequence FASTA files, lack of an index will cause failures or incorrect results

OUTPUT:
    - Summary to stderr
    - TSV file with counts of mono-, di-, tri-, tetra-nucleotide NTAs and totals
    - Plots in SVG and PDF format:
        • Combined proportional + bit score sequence logo (top: proportion, bottom: bit score)
        • Proportional bar plot by type (for QC)

DEPENDENCIES:
    - pysam
    - pandas
    - matplotlib
    - numpy
    - logomaker >= 0.8

LOGO PLOT NOTES:
    - Sequence logo plots normalize bit scores per position using information content
    - Positions start from 1, not 0
    - Mono, Di, Tri, and Tetra logos shown side-by-side with relative widths (1:2:3:4)

TOY EXAMPLES:
    Reference: ...TCG

    1. Read ends with soft-clip: G
         - If ref = G → Not counted (could be templated)
         - If ref != G → Counted as Mono-NTA

    2. Read ends with GG
         - Ref = CG → GG != CG → Counted as Di-NTA
         - Ref = GG → GG == GG → Not counted (matches template)

    3. Read ends with GTG, Ref = ATG:
         - G (last) != A (ref) → mismatch
         - T == T (ref)
         - G == G (ref)
         - → Counted as Tri-NTA (≥1 mismatch)

    4. Read ends with ATG, Ref = ATG:
         - All match → Not counted (potential templated match)
         - Important note: Cannot distinguish NTA from perfect templated match
"""


import sys
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker as lm
from collections import Counter
from pysam import FastaFile

bases = ['A', 'C', 'G', 'T']

COLOR_SCHEME = {
    'A': '#FDB863',
    'C': '#5E3C99',
    'G': '#E66101',
    'T': '#4393C3'
}

def reverse_complement(seq):
    complement = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(complement)[::-1]

def get_logo_matrix(seqs, value_col='Proportion'):
    seqs = [s for s in seqs if all(base in bases for base in s)]
    if not seqs:
        return pd.DataFrame(0.0, index=[], columns=bases)

    k = len(seqs[0])
    matrix = pd.DataFrame(0.0, index=np.arange(k), columns=bases)
    for s in seqs:
        for i, base in enumerate(s):
            matrix.at[i, base] += 1
    matrix = matrix.div(matrix.sum(axis=1), axis=0).fillna(0)
    if value_col == 'BitScore':
        entropy = - (matrix * np.log2(matrix + 1e-10)).sum(axis=1)
        info_content = np.log2(4) - entropy
        matrix = matrix.mul(info_content, axis=0)
    matrix.index = matrix.index.astype(int) + 1
    return matrix

def plot_combined_logos(mono, di, tri, tetra, output_prefix):
    sets = [(mono, 'Mono', 1), (di, 'Di', 2), (tri, 'Tri', 3), (tetra, 'Tetra', 4)]
    widths = [1, 2, 3, 4]
    total_width = sum(widths)
    rel_widths = [w / total_width for w in widths]

    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(2, 4, width_ratios=rel_widths)
    axs = np.empty((2, 4), dtype=object)
    for i in range(2):
        for j in range(4):
            axs[i, j] = fig.add_subplot(gs[i, j])

    for col, (seqs, label, length) in enumerate(sets):
        if not seqs:
            axs[0, col].axis('off')
            axs[1, col].axis('off')
            continue
        prop_mat = get_logo_matrix([s[:length] for s in seqs if len(s) >= length], value_col='Proportion')
        bit_mat = get_logo_matrix([s[:length] for s in seqs if len(s) >= length], value_col='BitScore')
        lm.Logo(prop_mat, ax=axs[0, col], color_scheme=COLOR_SCHEME)
        lm.Logo(bit_mat, ax=axs[1, col], color_scheme=COLOR_SCHEME)
        axs[0, col].set_title(f"{label}-NTA")
        axs[0, col].set_ylabel("Proportion")
        axs[1, col].set_ylabel("Bit Score")
        axs[1, col].set_xlabel("Position")
        positions = [i + 1 for i in range(length)]
        axs[0, col].set_xticks(positions)
        axs[0, col].set_xticklabels([str(i) for i in positions])
        axs[0, col].set_xlim(0.5, length + 0.5)
        axs[1, col].set_xticks(positions)
        axs[1, col].set_xticklabels([str(i) for i in positions])
        axs[1, col].set_xlim(0.5, length + 0.5)

    plt.tight_layout()
    fig.savefig(f"{output_prefix}_nta_combined_logo.svg")
    fig.savefig(f"{output_prefix}_nta_combined_logo.pdf")
    plt.close()


def usage():
    sys.stderr.write("""\nUsage: tailclip_filter.py <input.bam> <noClip.bam> <withClip.bam> <reference.fa> [output_prefix]\n""")
    sys.exit(1)

def process(input_bam, no_clip_bam, with_clip_bam, reference_fa, output_prefix):
    sam = pysam.AlignmentFile(input_bam, "rb")
    no_clip = pysam.AlignmentFile(no_clip_bam, "wb", template=sam)
    with_clip = pysam.AlignmentFile(with_clip_bam, "wb", template=sam)
    ref = FastaFile(reference_fa)

    mono, di, tri, tetra = [], [], [], []
    mono_counts = Counter()
    di_counts = Counter()
    tri_counts = Counter()
    tetra_counts = Counter()

    total, no_clip_count, with_clip_count = 0, 0, 0

    for read in sam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        total += 1
        cig = read.cigartuples or []
        if not cig or cig[-1][0] not in (4, 5):
            no_clip.write(read)
            no_clip_count += 1
            continue

        with_clip.write(read)
        with_clip_count += 1

        if cig[-1][0] != 4:
            continue

        clip_len = cig[-1][1]
        seq = read.query_sequence or ''
        if clip_len < 1 or len(seq) < clip_len:
            continue

        tail = seq[-clip_len:].upper()
        ref_start = read.reference_end
        chrom = read.reference_name

        try:
            ref_seq = ref.fetch(chrom, ref_start, ref_start + clip_len).upper()
            if len(ref_seq) < clip_len:
                continue
        except (ValueError, IndexError):
            continue

        if read.is_reverse:
            tail = reverse_complement(tail)
            ref_seq = reverse_complement(ref_seq)

        if clip_len == 1:
            if tail[-1] != ref_seq[-1]:
                mono_counts[tail[-1]] += 1
                mono.append(tail[-1])

        elif clip_len == 2:
            if tail[-1] != ref_seq[-1] or tail[-2] != ref_seq[-2]:
                di_counts[tail[-2:]] += 1
                di.append(tail[-2:])

        elif clip_len == 3:
            mismatches = sum(tail[-i] != ref_seq[-i] for i in (1, 2, 3))
            if mismatches >= 1:
                tri_counts[tail[-3:]] += 1
                tri.append(tail[-3:])

        elif clip_len >= 4:
            mismatches = sum(tail[-i] != ref_seq[-i] for i in (1, 2, 3, 4))
            if mismatches >= 1:
                tetra_counts[tail[-4:]] += 1
                tetra.append(tail[-4:])

    sys.stderr.write(f"Total mapped reads: {total}\n")
    sys.stderr.write(f"No clip reads     : {no_clip_count}\n")
    sys.stderr.write(f"Clipped reads     : {with_clip_count}\n\n")

    mono_total = sum(mono_counts.values())
    di_total = sum(di_counts.values())
    tri_total = sum(tri_counts.values())
    tetra_total = sum(tetra_counts.values())
    all_total = mono_total + di_total + tri_total + tetra_total

    out_tsv = f"{output_prefix}_nta_counts.tsv"
    with open(out_tsv, "w") as out:
        out.write("Type\tSequence\tCount\n")
        for b in bases:
            out.write(f"Mono\t{b}\t{mono_counts.get(b, 0)}\n")
        for k, v in sorted(di_counts.items()):
            out.write(f"Di\t{k}\t{v}\n")
        for k, v in sorted(tri_counts.items()):
            out.write(f"Tri\t{k}\t{v}\n")
        for k, v in sorted(tetra_counts.items()):
            out.write(f"Tetra\t{k}\t{v}\n")
        out.write(f"\nTotal\tMono\t{mono_total}\n")
        out.write(f"Total\tDi\t{di_total}\n")
        out.write(f"Total\tTri\t{tri_total}\n")
        out.write(f"Total\tTetra\t{tetra_total}\n")
        out.write(f"Total\tAll\t{all_total}\n")

    sys.stderr.write(f"\nCounts written to TSV: {out_tsv}\n")
    plot_combined_logos(mono, di, tri, tetra, output_prefix)

if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        usage()
    input_bam, no_clip_bam, with_clip_bam, ref_fa = sys.argv[1:5]
    prefix = sys.argv[5] if len(sys.argv) == 6 else with_clip_bam.rsplit('.', 1)[0]
    process(input_bam, no_clip_bam, with_clip_bam, ref_fa, prefix)
