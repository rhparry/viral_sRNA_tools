#!/usr/bin/env python3
"""
bam_sRNA_histogram.py
v0.8 Rhys Parry <r.parry@uq.edu.au> University of Queensland

Generates small RNA base composition histograms from BAM files:
  • Per-contig TSV with sense/antisense base and size counts
  • Optional: One figure per contig with stacked bar chart (A/C/G/T)
  • Figures saved as multipage PDF and individual SVGs

Usage:
    python bam_sRNA_histogram.py <input.bam> [-n <read_count>] [--min MIN] [--max MAX] [--plot] [--end 5|3]
"""

import argparse
import pysam
import os
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("bam", help="Input BAM file (indexed)")
    p.add_argument("-n", type=int, default=None, help="Original FASTQ read count for CPM normalization")
    p.add_argument("--min", type=int, default=15, help="Minimum read length to count")
    p.add_argument("--max", type=int, default=50, help="Maximum read length to count")
    p.add_argument("--plot", action="store_true", help="Output plots (PDF and SVG")
    p.add_argument("--end", choices=["5", "3"], default="5", help="Choose whether to use 5' or 3' end for plotting")
    return p.parse_args()

def init_base_dict():
    return {b: defaultdict(int) for b in "ACGT"}

def plot_stacked_bar(df, outprefix, end="5", is_cpm=False):
    colors = {"A": "#d62728", "C": "#2ca02c", "G": "#ff7f0e", "T": "#1f77b4"}
    bases = ["A", "C", "G", "T"]
    
    # Ensure proper ordering by converting Length to int, sorting, then back to string
    df_sorted = df.sort_values('Length').reset_index(drop=True)
    lengths = df_sorted["Length"].astype(str)

    fig_cpm, (ax_pos, ax_neg) = plt.subplots(nrows=2, ncols=1, figsize=(10, 6), sharex=True)
    bottom_pos = [0] * len(df_sorted)
    bottom_neg = [0] * len(df_sorted)

    for b in bases:
        pos = df_sorted[f"Sense_{b}_{'first' if end == '5' else 'last'}"]
        neg = df_sorted[f"Antisense_{b}_{'first' if end == '5' else 'last'}"] * -1
        ax_pos.bar(lengths, pos, bottom=bottom_pos, color=colors[b], label=b)
        ax_neg.bar(lengths, neg, bottom=bottom_neg, color=colors[b])
        bottom_pos = [sum(x) for x in zip(bottom_pos, pos)]
        bottom_neg = [sum(x) for x in zip(bottom_neg, neg)]

    ax_pos.set_title(f"Sense strand (Forward reads, {end}' end)")
    ax_neg.set_title(f"Antisense strand (Reverse reads, {end}' end)")
    ax_neg.set_xlabel("Length (nt)")
    ax_pos.set_ylabel("CPM" if is_cpm else "Count")
    ax_neg.set_ylabel("CPM" if is_cpm else "Count")
    ax_pos.legend(title="Base", bbox_to_anchor=(1.05, 1), loc='upper left')
    ax_pos.axhline(0, color="black", linewidth=0.5)
    ax_neg.axhline(0, color="black", linewidth=0.5)
    fig_cpm.suptitle(os.path.basename(outprefix), fontsize=12)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig_cpm.savefig(outprefix + f"_{end}prime.svg", dpi=300, bbox_inches='tight')
    fig_cpm.savefig(outprefix + f"_{end}prime.pdf", dpi=300, bbox_inches='tight')

    return [fig_cpm]

def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")
    # Calculate CPM modifier: if -n provided, normalize to counts per million
    modifier = 1e6 / args.n if args.n else 1

    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    all_figs = []
    outprefix_base = os.path.basename(args.bam).replace(".bam", "")

    for ref in bam.references:
        # Separate counts for forward-mapping (sense) and reverse-mapping (antisense) reads
        sense_count = defaultdict(int)      # Forward reads (-F 16)
        antisense_count = defaultdict(int)  # Reverse reads (-f 16)
        
        sense_first = init_base_dict()      # 5' end bases for forward reads
        sense_last = init_base_dict()       # 3' end bases for forward reads
        antisense_first = init_base_dict()  # 5' end bases for reverse reads  
        antisense_last = init_base_dict()   # 3' end bases for reverse reads

        for r in bam.fetch(ref):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            seq = r.query_sequence.upper()
            if not seq:
                continue

            L = len(seq)
            if L < args.min or L > args.max:
                continue

            # Follow the bash script logic exactly:
            # For both forward and reverse reads, seq[0] gives us the 5' end of the original RNA
            # seq[-1] gives us the 3' end of the original RNA
            if r.is_reverse:
                # Reverse reads: BAM stores as reverse complement
                # seq[0] is the 5' end of original RNA (what we want for "first")
                # seq[-1] is the 3' end of original RNA (what we want for "last")
                antisense_count[L] += 1
                f5_base = seq[0]   # 5' end of original RNA
                l3_base = seq[-1]  # 3' end of original RNA
                if f5_base in antisense_first:
                    antisense_first[f5_base][L] += 1
                if l3_base in antisense_last:
                    antisense_last[l3_base][L] += 1
            else:
                # Forward reads: BAM stores as-is
                # seq[0] is the 5' end of original RNA
                # seq[-1] is the 3' end of original RNA  
                sense_count[L] += 1
                f5_base = seq[0]   # 5' end of original RNA
                l3_base = seq[-1]  # 3' end of original RNA
                if f5_base in sense_first:
                    sense_first[f5_base][L] += 1
                if l3_base in sense_last:
                    sense_last[l3_base][L] += 1

        rows = []
        for L in range(args.min, args.max + 1):
            row = [L, sense_count[L] * modifier]
            row += [sense_first[b][L] * modifier for b in "ACGT"]
            row += [sense_last[b][L] * modifier for b in "ACGT"]
            row += [100 * sense_first[b][L] / sense_count[L] if sense_count[L] else 0 for b in "ACGT"]
            row += [100 * sense_last[b][L] / sense_count[L] if sense_count[L] else 0 for b in "ACGT"]
            row += [antisense_count[L] * modifier]
            row += [antisense_first[b][L] * modifier for b in "ACGT"]
            row += [antisense_last[b][L] * modifier for b in "ACGT"]
            row += [100 * antisense_first[b][L] / antisense_count[L] if antisense_count[L] else 0 for b in "ACGT"]
            row += [100 * antisense_last[b][L] / antisense_count[L] if antisense_count[L] else 0 for b in "ACGT"]
            rows.append(row)

        hdr = ["Length", "Sense_Count"] + \
              [f"Sense_{b}_first" for b in "ACGT"] + [f"Sense_{b}_last" for b in "ACGT"] + \
              [f"Sense_{b}_first_pct" for b in "ACGT"] + [f"Sense_{b}_last_pct" for b in "ACGT"] + \
              ["Antisense_Count"] + \
              [f"Antisense_{b}_first" for b in "ACGT"] + [f"Antisense_{b}_last" for b in "ACGT"] + \
              [f"Antisense_{b}_first_pct" for b in "ACGT"] + [f"Antisense_{b}_last_pct" for b in "ACGT"]

        df = pd.DataFrame(rows, columns=hdr)
        base_out = os.path.join(output_dir, f"{outprefix_base}_{ref}")
        df.to_csv(base_out + ".tsv", sep="\t", index=False)

        if args.plot:
            figs = plot_stacked_bar(df, base_out, end=args.end, is_cpm=bool(args.n))
            all_figs.extend(figs)

    if args.plot and all_figs:
        pdf_path = os.path.join(output_dir, f"{outprefix_base}_summary_{args.end}prime.pdf")
        with PdfPages(pdf_path) as pdf:
            for fig in all_figs:
                pdf.savefig(fig, dpi=300, bbox_inches='tight')
                plt.close(fig)

if __name__ == '__main__':
    main()
