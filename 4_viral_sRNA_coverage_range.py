#!/usr/bin/env python3
"""
4_viral_sRNA_coverage_range.py
v0.4.2 Rhys Parry <r.parry@uq.edu.au> University of Queensland

Generates strand-specific coverage from a BAM file, normalised to CPM.
Output: Tab-delimited file with columns per chromosome and rows per nucleotide position.
Supports:
  • Coverage from a defined size range (-r MIN MAX)
  • or Coverage from a size list file (-n sizes.txt)
  • X-axis only extends to each chromosome's actual length

Usage:
    python 4_viral_sRNA_coverage_range.py -r 18 30 -m 31707303 input.bam

    # -r MIN MAX : include reads between MIN and MAX nt
    # -m INT     : total mapped reads (to convert to CPM)
    # input.bam  : coordinate-sorted BAM file
"""

import pysam
import argparse
from collections import defaultdict
import os
import re

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("-r", nargs=2, type=int, metavar=('MIN', 'MAX'), help="Size range of reads to extract")
    parser.add_argument("-n", help="File with list of read lengths to include")
    parser.add_argument("-m", type=int, required=True, help="Total mapped reads for CPM normalization")
    return parser.parse_args()

def extract_reads_by_size(bamfile, output_bam, min_size=None, max_size=None, size_list=None):
    with pysam.AlignmentFile(bamfile, "rb") as infile, \
         pysam.AlignmentFile(output_bam, "wb", template=infile) as out:

        for read in infile:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            read_len = read.query_length
            if size_list and read_len in size_list:
                out.write(read)
            elif min_size is not None and max_size is not None and min_size <= read_len <= max_size:
                out.write(read)
    pysam.index(output_bam)

def get_chromosome_coverage(bamfile, chrom, modifier):
    coverage = defaultdict(lambda: [0, 0])  # position: [sense, antisense]
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            strand = 1 if read.is_reverse else 0
            for p in read.get_reference_positions():
                coverage[p + 1][strand] += 1
    for pos in coverage:
        coverage[pos][0] = round(coverage[pos][0] * modifier, 3)
        coverage[pos][1] = round(coverage[pos][1] * modifier, 3)
    return coverage

def calculate_coverage_all_chromosomes(bamfile, total_reads):
    base = os.path.basename(bamfile)
    prefix = re.sub(r"\.bam$", "", base)
    output_bam = prefix + ".filtered.bam"
    output_tsv = prefix + "_coverage_matrix.tsv"

    chromosomes = list(pysam.AlignmentFile(bamfile, "rb").references)
    modifier = 1e6 / total_reads

    if args.n:
        with open(args.n) as f:
            size_list = set(int(x.strip()) for x in f if x.strip().isdigit())
        extract_reads_by_size(bamfile, output_bam, size_list=size_list)
    elif args.r:
        min_size, max_size = args.r
        extract_reads_by_size(bamfile, output_bam, min_size=min_size, max_size=max_size)
    else:
        raise ValueError("You must provide either -n or -r")

    all_coverage = {}
    max_pos_per_chr = {}

    with pysam.AlignmentFile(output_bam, "rb") as bam:
        chr_lengths = dict(zip(bam.references, bam.lengths))

    for chr_name in chromosomes:
        coverage_data = get_chromosome_coverage(output_bam, chr_name, modifier)
        all_coverage[chr_name] = coverage_data
        max_pos_per_chr[chr_name] = chr_lengths[chr_name]

    max_position_overall = max(max_pos_per_chr.values())

    with open(output_tsv, "w") as f:
        headers = []
        for chr_name in chromosomes:
            headers.extend([f"{chr_name}_position", f"{chr_name}_sense", f"{chr_name}_antisense"])
        f.write("\t".join(headers) + "\n")

        for pos in range(1, max_position_overall + 1):
            row_data = []
            for chr_name in chromosomes:
                if pos <= max_pos_per_chr[chr_name]:
                    if pos in all_coverage[chr_name]:
                        sense_cov, antisense_cov = all_coverage[chr_name][pos]
                        row_data.extend([str(pos), str(sense_cov), str(-antisense_cov)])
                    else:
                        row_data.extend([str(pos), '0', '0'])
                else:
                    row_data.extend(['', '', ''])
            f.write('\t'.join(row_data) + '\n')

    os.remove(output_bam)
    if os.path.exists(output_bam + ".bai"):
        os.remove(output_bam + ".bai")

if __name__ == "__main__":
    args = parse_arguments()
    calculate_coverage_all_chromosomes(args.bam, args.m)
