# viral_sRNA_tools
This provides a range of bash scripts to examine small viral RNA signatures. To ensure compatibility for the BAM file with each script, I recommend running these scripts in the order they have been numbered.

## 1_fastq_histogram.sh
Takes a fastq file (gzipped or otherwise), calculates a histogram of read lengths, and outputs a table of the first nucleotide sequence.

Usage 

>bash 1_fastq_histogram.sh input.fastq.gz
>
Or, if you want to run for every file in the current directory:
>for f in *.fastq; do bash 1_fastq_histogram_updated.sh $f > ${f%.fastq}_fqhisto.txt; done

## 2_mapping_vRNAs.sh
This script maps a fastq file against a reference using bowtie2, suppresses unaligned reads, and outputs a sorted bam file with only mapped reads. Tested with bowtie2 v2.4.5 and samtools v1.13.

Usage

>bash 2_mapping_vRNAs.sh input.fastq reference.fasta

## 3_bam_sRNA_histogram.sh 
This script takes a sorted BAM file (SFTS_S.bam) and provides a histogram output for mapped read lengths and the first nucleotide of both forward and reverse orientation.

This has been tested with samtools v1.13 and does not require an index file.

Usage:

>bash 3_bam_sRNA_histogram.sh input.bam
```
Example output:
Chromosome: KP202165.1_SFTS_virus_HB29_S
Sense Read Sizes:
Size Counts A T G C A% T% G% C%
---- ------ - - - - -- -- -- --
18 8 6 2 0 0 0.027 0.009 0.000 0.000
19 7 0 3 2 2 0.000 0.013 0.009 0.009
20 6 0 6 0 0 0.000 0.027 0.000 0.000
21 18 2 2 12 2 0.009 0.009 0.054 0.009
22 34 4 7 10 13 0.018 0.031 0.045 0.058
23 16 0 11 3 2 0.000 0.049 0.013 0.009
25 11 0 5 6 0 0.000 0.022 0.027 0.000
26 44 9 26 5 4 0.040 0.117 0.022 0.018
27 27 5 4 13 5 0.022 0.018 0.058 0.022
28 18 10 6 2 0 0.045 0.027 0.009 0.000
29 16 3 7 6 0 0.013 0.031 0.027 0.000
30 8 2 5 1 0 0.009 0.022 0.004 0.000
31 3 0 0 0 3 0.000 0.000 0.000 0.013
32 4 4 0 0 0 0.018 0.000 0.000 0.000
33 1 1 0 0 0 0.004 0.000 0.000 0.000
35 1 0 0 0 1 0.000 0.000 0.000 0.004
41 1 0 1 0 0 0.000 0.004 0.000 0.000

Antisense Read Sizes:
Size Counts A T G C A% T% G% C%
---- ------ - - - - -- -- -- --
18 12 2 3 7 0 0.006 0.010 0.023 0.000
19 19 6 11 0 2 0.019 0.036 0.000 0.006
20 37 3 6 13 15 0.010 0.019 0.042 0.049
21 23 3 3 5 12 0.010 0.010 0.016 0.039
22 44 24 7 4 9 0.078 0.023 0.013 0.029
23 11 0 6 0 5 0.000 0.019 0.000 0.016
24 17 2 5 6 4 0.006 0.016 0.019 0.013
25 9 7 0 2 0 0.023 0.000 0.006 0.000
26 35 5 11 3 16 0.016 0.036 0.010 0.052
27 29 11 0 12 6 0.036 0.000 0.039 0.019
28 19 8 4 2 5 0.026 0.013 0.006 0.016
29 19 1 10 0 8 0.003 0.032 0.000 0.026
30 24 1 9 5 9 0.003 0.029 0.016 0.029
31 2 2 0 0 0 0.006 0.000 0.000 0.000
32 7 2 0 0 5 0.006 0.000 0.000 0.016
35 1 0 0 0 1 0.000 0.000 0.000 0.003

```
## 4_viral_sRNA_coverage.sh
Takes a BAM file and subsets the bam file into 2 BAM files, one for vsiRNAs (21nt) mapped reads and another for vpiRNAs (25-30nts). After producing the subsetted BAM files using bedtools, it calculates coverage for each chromosome within the BAM file. It produces three files per chromosome using the following nomenclature: {chromname}_21nt/piRNA_coverage_neg.tab {chromname}_21nt/piRNA_coverage_pos.tab and a {chromsomename}_combined coverage.tab.

Tested with samtools v1.13 and bedtools v2.30.0

Usage:

>bash 4_viral_sRNA_coverage.sh input.bam

## 5_extract_fw_rv_fasta.sh
Takes a sorted BAM file as input and outputs two fasta files, one for the forward reads and one for the reverse reads per chromosome, using the following nomenclature:
{bamfileprefix}_{chromsomename}_fw.fa and {bamfileprefix}_{chromsomename}_rv.fa. By default, it outputs all reads. However, you can specify a [-min=<MIN_SIZE>] [-max=<MAX_SIZE>]

Tested with with samtools v 1.13.

Usage:

>bash 5_extract_fw_rv_fasta.sh -i=<input.bam> [-min=<MIN_SIZE>] [-max=<MAX_SIZE>]
