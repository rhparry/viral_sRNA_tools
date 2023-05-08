#!/bin/bash
#Rhys Parry r.parry@uq.edu.au v0.2 
#Tested with bedtools v2.30.0 and samtools v1.13
#Takes a BAM file (tested with samtools v1.13, bedtools v2.30.0) and it produces 2 BAM files, one for vsiRNAs (21nt) mapped reads, and another for vpiRNAs (25-30nts). It produces three files per chromosome for both vsiRNA and vpiRNA BAM files: coverage_neg.tab coverage_pos.tab and a combined coverage.tab.
#Usage: bash viral_sRNA_coverage.sh input.bam

# Input BAM file
input_bam=$1

# Output BAM files
output_bam_21nt="${input_bam%.bam}_21nt.bam"
output_bam_25_30nt="${input_bam%.bam}_25_30nt.bam"

# Extract reads with length 21nt
samtools view -h $input_bam | awk 'length($10) == 21 || $1 ~ /^@/' | samtools view -b -o $output_bam_21nt -

# Index output BAM files
samtools index $output_bam_21nt

# Extract reads with length between 25 and 30nt
samtools view -h $input_bam | awk 'length($10) >= 25 && length($10) <= 30 || $1 ~ /^@/' | samtools view -b -o $output_bam_25_30nt -

# Index output BAM files
samtools index $output_bam_25_30nt

# Calculate coverage depth for each chromosome
chromosomes=$(samtools idxstats $input_bam | cut -f 1 | grep -v '*')

for chr in $chromosomes; do
    # Calculate coverage depth for 21nt BAM file
    bedtools genomecov -ibam $output_bam_21nt -d -split -strand + | awk -v OFS='\t' '{print $2,$3}' > ${chr}_21nt_coverage_pos.tab
    bedtools genomecov -ibam $output_bam_21nt -d -split -strand - | awk -v OFS='\t' '{print $2,$3}' > ${chr}_21nt_coverage_neg.tab
    join ${chr}_21nt_coverage_pos.tab ${chr}_21nt_coverage_neg.tab > ${chr}_21nt_coverage.tab

    # Calculate coverage depth for 25-30nt BAM file
    bedtools genomecov -ibam $output_bam_25_30nt -d -split -strand + | awk -v OFS='\t' '{print $2,$3}' > ${chr}_piRNA_coverage_pos.tab
    bedtools genomecov -ibam $output_bam_25_30nt -d -split -strand - | awk -v OFS='\t' '{print $2,$3}' > ${chr}_piRNA_coverage_neg.tab
    join ${chr}_piRNA_coverage_pos.tab ${chr}_piRNA_coverage_neg.tab > ${chr}_piRNA_coverage.tab
done
