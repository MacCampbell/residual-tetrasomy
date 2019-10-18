# magic-eight
For identification of likely tetrasomic homeologs versus likely disomic homeologs in salmonids

## Goals
A plot of overall similarity between homeologous pairs based on protokaryotype notation with classification into "tetrasomic" or "disomic" groups using any high-quality assembly.

## Files
  * Should be designed so you can execute the files in the /magic-eight/ directory
  * 10x.x-xxxxx.xx sequentially ordered files of what to do
  * These will basically run at first by using files in /data/ then putting outputs in /outputs/10x
  * The scripts should receive input from previous scripts and output into their own /outputs/xxx

## File Structure

/data/ genomes by sub directory (in .gitignore as /data/*)

/outputs/ outputs (in .gitignore as /outputs/*)

# Preparing Input Data
## Creating chromosome fasta files
In order to work with NCBI genome files for our purposes, unscaffolded contigs need to be removed. Here is an example of how to prepare one for _Salmo salar_.

1. Download GCA_000233375.4_ICSASG_v2_genomic.fna  via the NCBI site
2. Get the names of the chromosome sequences
  * $ grep "ssa" GCA_000233375.4_ICSASG_v2_genomic.fna > names.txt
  * $ xargs samtools faidx GCA_000233375.4_ICSASG_v2_genomic.fna < names.txt > salmo-salar-chroms.fasta
3. This file should be placed in /data/salmo-salar/
4. You can further split as needed like so:
  * $ samtools faidx salmoChroms.fasta ssa17 > salmo-salar-ssa17.fasta
  * Or with the *pairs.txt files: cut -f 1 -d ' ' t-thymallus-pairs.txt | while read line; do samtools faidx GCA_004348285.1_ASM434828v1_genomic.fna \$line > t-thymallus-\$line.fasta; done;

## Identifying centromeres and further splitting
For Chinook salmon, centromere information was available for Ots04 and Ots12 (Christensen pers. comm.). It was split like so:

  * samtools faidx GCF_002872995.1_Otsh_v1.0_genomic.fna NC_037100.1:1-41977851 > o-tshaw-NC_037100.1p-fasta
  * samtools faidx GCF_002872995.1_Otsh_v1.0_genomic.fna NC_037100.1:55780692-74299616 > o-tshaw-NC_037100.1q-fasta

  * samtools faidx GCF_002872995.1_Otsh_v1.0_genomic.fna NC_037108.1:1-44536220 > o-tshaw-NC_037108.1p.fasta
  * samtools faidx GCF_002872995.1_Otsh_v1.0_genomic.fna NC_037108.1:51430624-77127610 > o-tshaw-NC_037108.1q.fasta


# Running lastz and Generating Figures
Usage notes are noted in each script at the beginning.

__101-lastz.pl__ Dependencies: lastz, perl, Rscript, gnu-parallel. This generates lastz alignments of each PK.

__201-classify-and-test.R__ Generates box plots of overall similarity of PK's while also classifying PK's into inheritance mode categories.

__201.1-classify-and-test-thymallus.R__ Similar to 201-classify-and-test.R except that _Thymallus_ chromosomes were aligned with lastz already named based on protokaryotypes.

__301-combine-pdf.sh__ Takes the output from 201 & 201.1 and merges into two separate .pdf files forming the bases of Figure 4 and Supplemental Document S8 in Blumstein et al..

