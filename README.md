# magic-eight
For characterization of tetrasomic regions in salmonids

## Goals
To outline how to generate plots of similarities between tetrasomically pairing regions of salmonid genomes using any high-quality assembly.

## Files
  * Should be designed so you can execute the files in the /magic-eight/ directory
  * 10x.x-xxxxx.xx Sequentially ordered files of what to do
  * These will basically run at first by using files in /data/ then putting outputs in /outputs/10x
  * 10(x+1) should use /outputs/10x/ as inputs and output in /outputs/10(x+1)/ and so on

## File Structure

/data/ genomes by sub directory (in .gitignore as /data/*)

/outputs/ outputs (in .gitignore as /outputs/*)

# Preparing Input Data
## Creating a chromosome fasta file
In order to work with NCBI genome files for our purposes, unscaffolded contigs need to be removed. We are storing the chromosome-only files on google drive. Here is an example of how to prepare one for _Salmo salar_.

1. Download GCA_000233375.4_ICSASG_v2_genomic.fna  via the NCBI site
2. Get the names of the chromosome sequences
  * $ grep "ssa" GCA_000233375.4_ICSASG_v2_genomic.fna > names.txt
  * $ xargs samtools faidx GCA_000233375.4_ICSASG_v2_genomic.fna < names.txt > salmo-salar-chroms.fasta
3. This file should be placed in /data/salmo-salar/
4. You can further split as needed like so:
  * $ samtools faidx salmoChroms.fasta ssa17 > salmo-salar-ssa17.fasta
  * $ samtools faidx salmoChroms.fasta ssa16 > salmo-salar-ssa16.fasta

## Identifying centromeres and further splitting
1. TBD - presumably with alignment to _Esox lucius_ or self

# Running lastz and Generating Figures

1. Dependencies: lastz, perl, Rscript, gnu-parallel


