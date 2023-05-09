# viral_sRNA_bash
This provides a range of bash scripts to examine small viral RNA signatures. To ensure you have compatibility for the BAM file with each script I would recommend running these scripts in the order that they have been numbered.

## 1 fastq_histogram.sh
Takes a fastq file (gzipped or otherwise) and calculates a histogram of read lengths and also outputs a table of the first nucleotide sequence.

Usage 

>bash fastq_histogram.sh input.fastq.gz

## 2 mapping_vRNAs.sh
Takes a fastq file and maps it against a reference using bowtie2 and suppresses unaligned reads and outputs a sorted bam file with only mapped reads. Tested with bowtie2 v2.4.5 and samtools v1.13.

Usage

>bash mapping_vRNAs.sh input.fastq reference.fasta

## 3 sRNA_histogram.sh 
Takes a sorted BAM file (in this case TOSV.bam) and provides a histogram output for mapped read lengths and also the first nucleotide of both forward and reverse orientation.

Tested with samtools v1.13, does not require index file.

Usage:

>bash sRNA_histogram.sh input.bam
```
Example output:
Chromosome: MT032308.1_Toscana_virus_isolate_1500590_segment_L
Sense Read Sizes:
Size  Counts  A  T  G  C
----  ------  -  -  -  -
18 327 147 61 78 41
19 625 246 167 133 79
20 889 285 283 190 131
21 4600 1722 1474 619 785
22 749 281 234 159 75
23 526 173 172 114 67
24 406 181 104 91 30
25 414 176 105 102 31
26 401 175 87 100 39
27 593 238 136 165 54
28 418 198 71 106 43
29 125 37 24 53 11
30 5 2  2 1
31 5 2 1  2
32 1   1
33 1  1
34 1  1

Antisense Read Sizes:
Size  Counts  A  T  G  C
----  ------  -  -  -  -
18 278 133 41 46 58
19 793 324 152 169 148
20 1118 441 210 268 199
21 4131 1306 699 1210 916
22 755 306 155 159 135
23 481 158 121 106 96
24 314 98 70 69 77
25 270 105 70 54 41
26 271 98 55 57 61
27 263 90 59 74 40
28 253 74 85 70 24
29 46 16 9 15 6
30 9 3  3 3
31 2 1  1
```
## 4 viral_sRNA_coverage.sh
Takes a BAM file and subsets the bam file into 2 BAM files, one for vsiRNAs (21nt) mapped reads, and another for vpiRNAs (25-30nts). After producing the subsetted BAM files using bedtools it calculates coverage for each chromosome within the BAM file and produces three files per chromosome using the following nomenclature: {chromsomename}_21nt/piRNA_coverage_neg.tab {chromsomename}_21nt/piRNA_coverage_pos.tab and a {chromsomename}_combined coverage.tab.

Tested with samtools v1.13 and bedtools v2.30.0

Usage:

>bash viral_sRNA_coverage.sh input.bam

## 5 extract_fw_rv_fasta.sh
Takes a sorted BAM file as input and outputs two fasta files, one for the forward reads and one for the reverse reads per chromosome using the following nomenclature:
{bamfileprefix}_{chromsomename}_fw.fa and {bamfileprefix}_{chromsomename}_rv.fa. 

Tested with with samtools v 1.13.

Usage:

>bash extract_fw_rv_fasta.sh input.bam
