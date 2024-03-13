#v0.1 Rhys H Parry r.parry@uq.edu.au
#Takes a fastq file and outputs the proportion of bias in the 5' and 3' for each read length size.
#Usage: python3 m5prime_3prime_bias.py input.fastq

import sys
from Bio import SeqIO
from collections import defaultdict

def count_nucleotides(fastq_file):
    # Initialize dictionaries to store counts
    start_counts = defaultdict(lambda: defaultdict(int))
    end_counts = defaultdict(lambda: defaultdict(int))
    total_counts = defaultdict(int)

    # Parse FASTQ file
    for record in SeqIO.parse(fastq_file, "fastq"):
        length = len(record.seq)
        start_nt = record.seq[0]
        end_nt = record.seq[-1]

        # Update counts
        start_counts[length][start_nt] += 1
        end_counts[length][end_nt] += 1
        total_counts[length] += 1

    # Calculate proportions
    start_proportions = {length: {nt: count / total_counts[length] for nt, count in nt_counts.items()} for length, nt_counts in start_counts.items()}
    end_proportions = {length: {nt: count / total_counts[length] for nt, count in nt_counts.items()} for length, nt_counts in end_counts.items()}

    return start_proportions, end_proportions

def print_proportions(proportions, label):
    print(label)
    print("Size\tA\tT\tG\tC")
    for length in sorted(proportions.keys()):
        print(f"{length}\t{proportions[length].get('A', 0):.2f}\t{proportions[length].get('T', 0):.2f}\t{proportions[length].get('G', 0):.2f}\t{proportions[length].get('C', 0):.2f}")

def main():
    fastq_file = sys.argv[1]  # Get FASTQ file from command-line argument
    start_proportions, end_proportions = count_nucleotides(fastq_file)
    print_proportions(start_proportions, "Start")
    print_proportions(end_proportions, "End")

if __name__ == "__main__":
    main()
