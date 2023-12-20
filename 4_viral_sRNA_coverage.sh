#!/bin/bash
#Rhys Parry r.parry@uq.edu.au v0.2.2
#Tested with bedtools v2.30.0 and samtools v1.13
#Takes a BAM file and by default produces 2 BAM files, one for vsiRNAs (21nt) mapped reads, and another for vpiRNAs (24-30nts). 
#It produces three files per chromosome for both vsiRNA and vpiRNA BAM files: filename_chromosomename_coverage_neg.txt filename_chromosomename_coverage_pos.txt and 
#a combined filename_chromosomename_coverage.txt
#You can also use the flags -n for a specific size or -r for a range (-r min max) 
#Usage: bash viral_sRNA_coverage.sh [options] input.bam
#This version of the script extracts 5' position and also outputs normalised coverage.


# Function to display usage information
usage() {
    echo "Usage: $0 [-n size] [-r min max] input.bam"
    echo "  -n size     Extract reads with length equal to size"
    echo "  -r min max  Extract reads with length between min and max (inclusive)"
    echo "  input.bam   Input BAM file"
    exit 1
}

# Parse command line options
while getopts ":n:r:h" opt; do
    case $opt in
        n)
            size=$OPTARG
            ;;
        r)
            min=$OPTARG
            max=${!OPTIND}
            OPTIND=$(( $OPTIND + 1 ))
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

# Check if the input BAM file is sorted
if ! samtools view -H $input_bam | grep -q '^@HD.*SO:coordinate'; then
    echo "Input BAM file is not sorted by coordinate."
    echo "Sorting input BAM file"
    samtools sort -o ${input_bam%.bam}_sorted.bam $input_bam
    input_bam="${input_bam%.bam}_sorted.bam"
fi

# Check if the index file for sorted input BAM file exists
if [ ! -f "${input_bam}.bai" ]; then
    echo "Index file for sorted input BAM file not found."
    echo "Indexing sorted input BAM file."
    samtools index $input_bam
fi

echo "Input BAM file is sorted and indexed"

# Get input file prefix
input_prefix=$(basename $input_bam .bam)

# Output BAM files and extract reads with specified size or range or default settings if no flags are used
if [ ! -z "$size" ]; then
    output_bam="${input_bam%.bam}_${size}nt.bam"
    echo "Extracting reads with length $size nt"
    samtools view -h $input_bam | awk -v size=$size 'length($10) == size || $1 ~ /^@/' | samtools view -b -o $output_bam -
elif [ ! -z "$min" ] && [ ! -z "$max" ]; then
    output_bam="${input_bam%.bam}_${min}_${max}nt.bam"
    echo "Extracting reads with length between $min and $max nt (inclusive)"
    samtools view -h $input_bam | awk -v min=$min -v max=$max 'length($10) >= min && length($10) <= max || $1 ~ /^@/' | samtools view -b -o $output_bam -
else
    output_bam_21nt="${input_bam%.bam}_21nt.bam"
    output_bam_24_30nt="${input_bam%.bam}_24_30nt.bam"

    echo "Extracting reads with default settings (21nt and 24-30nt)"
    samtools view -h $input_bam | awk 'length($10) == 21 || $1 ~ /^@/' | samtools view -b -o $output_bam_21nt -
    samtools view -h $input_bam | awk 'length($10) >= 24 && length($10) <= 30 || $1 ~ /^@/' | samtools view -b -o $output_bam_24_30nt -
fi

echo "Indexing output BAM files."
# Index output BAM files
if [ ! -z "$size" ] || ([ ! -z "$min" ] && [ ! -z "$max" ]); then
    samtools index $output_bam
else
    samtools index $output_bam_21nt
    samtools index $output_bam_24_30nt
fi

echo "Calculating coverage depth for each chromosome."
# Calculate coverage depth for each chromosome and specified size or range of reads or default settings if no flags are used.
chromosomes=$(samtools idxstats $input_bam | cut -f 1 | grep -v '*')

for chr in $chromosomes; do
    echo "Processing chromosome $chr"

    if [ ! -z "$size" ] || ([ ! -z "$min" ] && [ ! -z "$max" ]); then
        bedtools genomecov -ibam $output_bam -d -strand + -scale 1.0 -5 | awk -v OFS='\t' '{print $2,$3}' > ${input_prefix}_${chr}_coverage_pos.txt
        bedtools genomecov -ibam $output_bam -d -strand - -scale 1.0 -5 | awk -v OFS='\t' '{print $2,$3}' > ${input_prefix}_${chr}_coverage_neg.txt
        join ${input_prefix}_${chr}_coverage_pos.txt ${input_prefix}_${chr}_coverage_neg.txt > ${input_prefix}_${chr}_coverage.txt
		awk '{f[NR]=$2;r[NR]=$3;sumf+=$2;sumr+=$3} END {for (i=1;i<=NR;i++) {printf "%d %.6f %.6f\n", i, f[i]/sumf, r[i]/sumr}}' ${input_prefix}_${chr}_coverage.txt > ${input_prefix}_${chr}_cov_norm.txt
    else
        bedtools genomecov -ibam $output_bam_21nt -d -strand + -scale 1.0 -5 | awk -v OFS='\t' '{print $2,$3}' > ${input_prefix}_${chr}_21nt_coverage_pos.txt
        bedtools genomecov -ibam $output_bam_21nt -d -strand - -scale 1.0 -5 | awk -v OFS='\t' '{print $2,$3}' > ${input_prefix}_${chr}_21nt_coverage_neg.txt
        join ${input_prefix}_${chr}_21nt_coverage_pos.txt ${input_prefix}_${chr}_21nt_coverage_neg.txt > ${input_prefix}_${chr}_21nt_coverage.txt
		awk '{f[NR]=$2;r[NR]=$3;sumf+=$2;sumr+=$3} END {for (i=1;i<=NR;i++) {printf "%d %.6f %.6f\n", i, f[i]/sumf, r[i]/sumr}}' ${input_prefix}_${chr}_21nt_coverage.txt > ${input_prefix}_${chr}_21nt_cov_norm.txt

        bedtools genomecov -ibam $output_bam_24_30nt -d -strand + -scale 1.0 -5 | awk -v OFS='\t' '{print $2,$3}' > ${input_prefix}_${chr}_piRNA_coverage_pos.txt
        bedtools genomecov -ibam $output_bam_24_30nt -d -strand - -scale 1.0 -5 | awk -v OFS='\t' '{print $2,$3}' > ${input_prefix}_${chr}_piRNA_coverage_neg.txt
        join ${input_prefix}_${chr}_piRNA_coverage_pos.txt ${input_prefix}_${chr}_piRNA_coverage_neg.txt > ${input_prefix}_${chr}_piRNA_coverage.txt
		awk '{f[NR]=$2;r[NR]=$3;sumf+=$2;sumr+=$3} END {for (i=1;i<=NR;i++) {printf "%d %.6f %.6f\n", i, f[i]/sumf, r[i]/sumr}}' ${input_prefix}_${chr}_piRNA_coverage.txt > ${input_prefix}_${chr}_piRNA_cov_norm.txt
		
    fi

done

echo "Done"
