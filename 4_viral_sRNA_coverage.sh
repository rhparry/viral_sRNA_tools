#!/bin/bash
# Rhys Parry r.parry@uq.edu.au v0.2.3
# Tested with bedtools v2.30.0 and samtools v1.13
# This script calculates coverage for specified read sizes or ranges in a BAM file.
# Normalization to counts per million (CPM) is optional via the -m flag.
# Usage:
#   bash 4_viral_sRNA_coverage_range.sh [-n size] [-r min max] [-m million_reads] input.bam
# Options:
#   -n size         Extract reads with length equal to size.
#   -r min max      Extract reads with length between min and max (inclusive).
#   -m million_reads Normalize coverage to counts per million (CPM).
#   input.bam       Input BAM file (sorted and indexed).

# Function to display usage information
usage() {
    echo "Usage: $0 [-n size] [-r min max] [-m million_reads] input.bam"
    echo "  -n size          Extract reads with length equal to size"
    echo "  -r min max       Extract reads with length between min and max (inclusive)"
    echo "  -m million_reads Normalize coverage to counts per million (CPM)"
    echo "  input.bam        Input BAM file (sorted and indexed)"
    exit 1
}

# Parse command line options
while getopts ":n:r:m:h" opt; do
    case $opt in
        n)
            size=$OPTARG
            ;;
        r)
            min=$OPTARG
            max=${!OPTIND}
            OPTIND=$(( $OPTIND + 1 ))
            ;;
        m)
            million_reads=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument" >&2
            usage
            ;;
    esac
done

# Shift options out of the way
shift $((OPTIND-1))

# Check if input BAM file is provided
if [ -z "$1" ]; then
    echo "Error: No input BAM file provided"
    usage
fi

# Input BAM file
input_bam=$1

# Check if input BAM file is sorted
if ! samtools view -H "$input_bam" | grep -q '^@HD.*SO:coordinate'; then
    echo "Input BAM file is not sorted by coordinate"
    echo "Sorting input BAM file"
    samtools sort -o "${input_bam%.bam}_sorted.bam" "$input_bam"
    input_bam="${input_bam%.bam}_sorted.bam"
fi

# Check if index file for sorted input BAM file exists
if [ ! -f "${input_bam}.bai" ]; then
    echo "Index file for sorted input BAM file not found"
    echo "Indexing sorted input BAM file"
    samtools index "$input_bam"
fi

echo "Input BAM file is sorted and indexed"

# Get input file prefix
input_prefix=$(basename "$input_bam" .bam)

# Output BAM files and extract reads with specified size or range
if [ ! -z "$size" ]; then
    output_bam="${input_bam%.bam}_${size}nt.bam"
    echo "Extracting reads with length $size nt"
    samtools view -h "$input_bam" | awk -v size="$size" 'length($10) == size || $1 ~ /^@/' | samtools view -b -o "$output_bam" -
elif [ ! -z "$min" ] && [ ! -z "$max" ]; then
    output_bam="${input_bam%.bam}_${min}_${max}nt.bam"
    echo "Extracting reads with length between $min and $max nt (inclusive)"
    samtools view -h "$input_bam" | awk -v min="$min" -v max="$max" 'length($10) >= min && length($10) <= max || $1 ~ /^@/' | samtools view -b -o "$output_bam" -
else
    echo "Error: You must specify a size (-n) or range (-r)"
    usage
fi

# Index the output BAM file
echo "Indexing output BAM file"
samtools index "$output_bam"

# Calculate coverage depth for each chromosome
echo "Calculating coverage depth for each chromosome"
chromosomes=$(samtools idxstats "$input_bam" | cut -f 1 | grep -v '*')

for chr in $chromosomes; do
    echo "Processing chromosome $chr"
    bedtools genomecov -ibam "$output_bam" -d -split -strand + | awk -v OFS='\t' '{print $2,$3}' > "${input_prefix}_${chr}_coverage_pos.txt"
    bedtools genomecov -ibam "$output_bam" -d -split -strand - | awk -v OFS='\t' '{print $2,$3}' > "${input_prefix}_${chr}_coverage_neg.txt"

    if [ ! -z "$million_reads" ]; then
        echo "Normalizing coverage to counts per million (CPM)"
        awk -v OFS='\t' -v million_reads="$million_reads" '{print $1, $2 / million_reads}' "${input_prefix}_${chr}_coverage_pos.txt" > "${input_prefix}_${chr}_coverage_pos_normalized.txt"
        awk -v OFS='\t' -v million_reads="$million_reads" '{print $1, $2 / million_reads}' "${input_prefix}_${chr}_coverage_neg.txt" > "${input_prefix}_${chr}_coverage_neg_normalized.txt"
        join "${input_prefix}_${chr}_coverage_pos_normalized.txt" "${input_prefix}_${chr}_coverage_neg_normalized.txt" > "${input_prefix}_${chr}_coverage_normalized.txt"
        rm "${input_prefix}_${chr}_coverage_pos_normalized.txt" "${input_prefix}_${chr}_coverage_neg_normalized.txt"
    else
        join "${input_prefix}_${chr}_coverage_pos.txt" "${input_prefix}_${chr}_coverage_neg.txt" > "${input_prefix}_${chr}_coverage.txt"
    fi

    rm "${input_prefix}_${chr}_coverage_pos.txt" "${input_prefix}_${chr}_coverage_neg.txt"
done

echo "Done"
