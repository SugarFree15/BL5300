# Finalized Pipline

## Downloading RNA-Seq Data from NCBI SRA Database

cd ~/scratch/BL5300/topic_of_choice

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump ERR3001915

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump ERR3001916

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump SRR10054928

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump SRR10054929

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump SRR10054931

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump SRR10054932

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump SRR10320043

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump SRR10320044

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump SRR10320045

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump ERR1727061

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump ERR1727060

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump SRR1610505

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump SRR1610506

~/software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump SRR1610507

> Using the fasterq-dump function in the NCBI SRA-toolkit, I downloaded each replicate of each tissues from a variety of studies in the NCBI Bioproject database, included paired runs, if applicable. NOTE: Not all replicates were used due to time it took to troubleshoot/proccess the samples.

## Downloading A. thaliana Genome and Info Files

mkdir Athal_ncbi/

> Accessed the TAIR10.1 Assembly of the A. thaliana genome through the RefSeq Genome FTP via Cyberduck, coping the genomic (GCF_000001735.4_TAIR10.1_genomic.fna) and info (GCF_000001735.4_TAIR10.1_genomic.gff) files to colossus in my ~/scratch/BL5300/topic_of_choice/Athal_ncbi directory.

## Quality Assessment and Control

fastqc <Sample_ID>

> Ran all the samples through fastqc to ensure quality. No trimming needed.

## Mapping Raw RNA-Seq Reads to the Genome and Counting Transcripts

> NOTE: I was slightly desperate to troubleshoot my syntaxes, so the working directory and some files are named with 'please' with no significance other than hoping this pipeline worked.

cd ~/scratch/BL5300/topic_of_choice/

mkdir please

cd please/

grep -v "unknown_transcript_1" ../Athal_ncbi/GCF_000001735.4_TAIR10.1_genomic.gtf > please_genomic.gtf

mkdir rsem

rsem-prepare-reference --gtf please_genomic.gtf --star ../Athal_ncbi/GCF_000001735.4_TAIR10.1_genomic.fna rsem/please_ncbi

> This indexed the TAIR10.1 genome in both RSEM and STAR at the same time using transcript info from a modified .gtf file that removed a set of unknown transcripts that would through errors from incorrect orientations in the files.

rsem-calculate-expression --no-bam-output --star ../SRR1610505.fastq rsem/please_ncbi SRR1610505

rsem-calculate-expression --no-bam-output --star ../SRR10320043.fastq rsem/please_ncbi 

rsem-calculate-expression --no-bam-output --star --paired-end --forward-prob 0 ../ERR3001915_1.fastq ../ERR3001915_2.fastq rsem/please_ncbi

rsem-calculate-expression --no-bam-output --star --paired-end --forward-prob 0 ../ERR1727060_1.fastq ../ERR1727060_2.fastq rsem/please_ncbi

rsem-calculate-expression --no-bam-output --star --paired-end --forward-prob 0 ../SRR10054931_1.fastq ../SRR10054931_2.fastq rsem/please_ncbi

> These took the raw reads and aligned them via STAR with parameters perfect for RSEM before then counting the transcripts and outputting the needed data with RSEM, taking into account whether or not each sample was a paired-end read.

## Data Configuration Before Analysis

grep -w "AT2G06200" SRR10320043.genes.results >> SRR10320043.GRFs.results

grep -w "AT2G22840" SRR10320043.genes.results >> SRR10320043.GRFs.results

grep -w "AT2G36400" SRR10320043.genes.results >> SRR10320043.GRFs.results

grep -w "AT2G45480" SRR10320043.genes.results >> SRR10320043.GRFs.results

grep -w "AT3G13960" SRR10320043.genes.results >> SRR10320043.GRFs.results

grep -w "AT3G52910" SRR10320043.genes.results >> SRR10320043.GRFs.results

grep -w "AT4G24150" SRR10320043.genes.results >> SRR10320043.GRFs.results

grep -w "AT4G37740" SRR10320043.genes.results >> SRR10320043.GRFs.results

grep -w "AT5G53660" SRR10320043.genes.results >> SRR10320043.GRFs.results


grep -w "AT2G06200" ERR3001915.genes.results >> ERR3001915.GRFs.results

grep -w "AT2G22840" ERR3001915.genes.results >> ERR3001915.GRFs.results

grep -w "AT2G36400" ERR3001915.genes.results >> ERR3001915.GRFs.results

grep -w "AT2G45480" ERR3001915.genes.results >> ERR3001915.GRFs.results

grep -w "AT3G13960" ERR3001915.genes.results >> ERR3001915.GRFs.results

grep -w "AT3G52910" ERR3001915.genes.results >> ERR3001915.GRFs.results

grep -w "AT4G24150" ERR3001915.genes.results >> ERR3001915.GRFs.results

grep -w "AT4G37740" ERR3001915.genes.results >> ERR3001915.GRFs.results

grep -w "AT5G53660" ERR3001915.genes.results >> ERR3001915.GRFs.results


grep -w "AT2G06200" ERR1727060.genes.results >> ERR1727060.GRFs.results

grep -w "AT2G22840" ERR1727060.genes.results >> ERR1727060.GRFs.results

grep -w "AT2G36400" ERR1727060.genes.results >> ERR1727060.GRFs.results

grep -w "AT2G45480" ERR1727060.genes.results >> ERR1727060.GRFs.results

grep -w "AT3G13960" ERR1727060.genes.results >> ERR1727060.GRFs.results

grep -w "AT3G52910" ERR1727060.genes.results >> ERR1727060.GRFs.results

grep -w "AT4G24150" ERR1727060.genes.results >> ERR1727060.GRFs.results

grep -w "AT4G37740" ERR1727060.genes.results >> ERR1727060.GRFs.results

grep -w "AT5G53660" ERR1727060.genes.results >> ERR1727060.GRFs.results


grep -w "AT2G06200" SRR10054931.genes.results >> SRR10054931.GRFs.results

grep -w "AT2G22840" SRR10054931.genes.results >> SRR10054931.GRFs.results

grep -w "AT2G36400" SRR10054931.genes.results >> SRR10054931.GRFs.results

grep -w "AT2G45480" SRR10054931.genes.results >> SRR10054931.GRFs.results

grep -w "AT3G13960" SRR10054931.genes.results >> SRR10054931.GRFs.results

grep -w "AT3G52910" SRR10054931.genes.results >> SRR10054931.GRFs.results

grep -w "AT4G24150" SRR10054931.genes.results >> SRR10054931.GRFs.results

grep -w "AT4G37740" SRR10054931.genes.results >> SRR10054931.GRFs.results

grep -w "AT5G53660" SRR10054931.genes.results >> SRR10054931.GRFs.results


grep -w "AT2G06200" SRR1610505.genes.results >> SRR1610505.GRFs.results

grep -w "AT2G22840" SRR1610505.genes.results >> SRR1610505.GRFs.results

grep -w "AT2G36400" SRR1610505.genes.results >> SRR1610505.GRFs.results

grep -w "AT2G45480" SRR1610505.genes.results >> SRR1610505.GRFs.results

grep -w "AT3G13960" SRR1610505.genes.results >> SRR1610505.GRFs.results

grep -w "AT3G52910" SRR1610505.genes.results >> SRR1610505.GRFs.results

grep -w "AT4G24150" SRR1610505.genes.results >> SRR1610505.GRFs.results

grep -w "AT4G37740" SRR1610505.genes.results >> SRR1610505.GRFs.results

grep -w "AT5G53660" SRR1610505.genes.results >> SRR1610505.GRFs.results


> Using the gene IDs of each GRF in Arabidopsis, I used grep to sort the GRF data from the RSEM results into new files that were then downloaded to my R.project directory via Cyberduck for analysis.

## Data Analysis in RStudio

library(tidyverse)

> These blocks of code format each tissue's data into a proper data frame and renaming the columns from the default "V#" titles to their actual labels, as well as plot the TPM and FPKM of each GRF gene as a point and color-coded in uniform for all.

RootCounts <- read.table("SRR1610505.GRFs.results")

RootCounts <- as.data.frame(RootCounts)

RootCounts <- mutate(RootCounts, tissue = "root")

RootCounts <- rename(RootCounts, "gene_id" = "V1")

RootCounts <- rename(RootCounts, "transcript_id(s)" = "V2")

RootCounts <- rename(RootCounts, "length" = "V3")

RootCounts <- rename(RootCounts, "effective_length" = "V4")

RootCounts <- rename(RootCounts, "expected_count" = "V5")

RootCounts <- rename(RootCounts, "TPM" = "V6")

RootCounts <- rename(RootCounts, "FPKM" = "V7")

ggplot(RootCounts)+
  geom_point(aes(x=gene_id, y=TPM, color=gene_id))+
  ggtitle("GRF Expression in Root Stele Cells")

ggplot(RootCounts)+
  geom_point(aes(x=gene_id, y=FPKM, color=gene_id))+
  ggtitle("GRF Expression in Root Stele Cells")

SeedCounts <- read.table("SRR10054931.GRFs.results")

SeedCounts <- as.data.frame(SeedCounts)

SeedCounts <- mutate(SeedCounts, tissue = "seed")

SeedCounts <- rename(SeedCounts, "gene_id" = "V1")

SeedCounts <- rename(SeedCounts, "transcript_id(s)" = "V2")

SeedCounts <- rename(SeedCounts, "length" = "V3")

SeedCounts <- rename(SeedCounts, "effective_length" = "V4")

SeedCounts <- rename(SeedCounts, "expected_count" = "V5")

SeedCounts <- rename(SeedCounts, "TPM" = "V6")

SeedCounts <- rename(SeedCounts, "FPKM" = "V7")

ggplot(SeedCounts)+
  geom_point(aes(x=gene_id, y=TPM, color=gene_id))+
  ggtitle("GRF Expression in Seed 12 Days Post-pollination")

ggplot(SeedCounts)+
  geom_point(aes(x=gene_id, y=FPKM, color=gene_id))+
  ggtitle("GRF Expression in Seed 12 Days Post-pollination")

CallusCounts <- read.table("ERR1727060.GRFs.results")

CallusCounts <- as.data.frame(CallusCounts)

CallusCounts <- mutate(CallusCounts, tissue = "callus")

CallusCounts <- rename(CallusCounts, "gene_id" = "V1")

CallusCounts <- rename(CallusCounts, "transcript_id(s)" = "V2")

CallusCounts <- rename(CallusCounts, "length" = "V3")

CallusCounts <- rename(CallusCounts, "effective_length" = "V4")

CallusCounts <- rename(CallusCounts, "expected_count" = "V5")

CallusCounts <- rename(CallusCounts, "TPM" = "V6")

CallusCounts <- rename(CallusCounts, "FPKM" = "V7")

ggplot(CallusCounts)+
  geom_point(aes(x=gene_id, y=TPM, color=gene_id))+
  ggtitle("GRF Expression in Calli")

ggplot(CallusCounts)+
  geom_point(aes(x=gene_id, y=FPKM, color=gene_id))+
  ggtitle("GRF Expression in Calli")

WholeCounts <- read.table("ERR3001915.GRFs.results")

WholeCounts <- as.data.frame(WholeCounts)

WholeCounts <- mutate(WholeCounts, tissue = "full organism")

WholeCounts <- rename(WholeCounts, "gene_id" = "V1")

WholeCounts <- rename(WholeCounts, "transcript_id(s)" = "V2")

WholeCounts <- rename(WholeCounts, "length" = "V3")

WholeCounts <- rename(WholeCounts, "effective_length" = "V4")

WholeCounts <- rename(WholeCounts, "expected_count" = "V5")

WholeCounts <- rename(WholeCounts, "TPM" = "V6")

WholeCounts <- rename(WholeCounts, "FPKM" = "V7")

ggplot(WholeCounts)+
  geom_point(aes(x=gene_id, y=TPM, color=gene_id))+
  ggtitle("GRF Expression Summary of Full Plant")

ggplot(WholeCounts)+
  geom_point(aes(x=gene_id, y=FPKM, color=gene_id))+
  ggtitle("GRF Expression Summary of Full Plant")

ApexCounts <- read.table("SRR10320043.GRFs.results")

ApexCounts <- as.data.frame(ApexCounts)

ApexCounts <- mutate(ApexCounts, tissue = "stem apex")

ApexCounts <- rename(ApexCounts, "gene_id" = "V1")

ApexCounts <- rename(ApexCounts, "transcript_id(s)" = "V2")

ApexCounts <- rename(ApexCounts, "length" = "V3")

ApexCounts <- rename(ApexCounts, "effective_length" = "V4")

ApexCounts <- rename(ApexCounts, "expected_count" = "V5")

ApexCounts <- rename(ApexCounts, "TPM" = "V6")

ApexCounts <- rename(ApexCounts, "FPKM" = "V7")

ggplot(ApexCounts)+
  geom_point(aes(x=gene_id, y=TPM, color=gene_id))+
  ggtitle("GRF Expression in Apical Meristem")

ggplot(ApexCounts)+
  geom_point(aes(x=gene_id, y=FPKM, color=gene_id))+
  ggtitle("GRF Expression in Apical Meristem")

> From here, I created a data frame with all of the data collected and created summary plots of TPM/FPKM of each gene per tissue and box plots of the overall TPM/FPKM of the study.

consensus <- rbind(RootCounts, CallusCounts, SeedCounts, WholeCounts, ApexCounts)

ggplot(consensus)+
  geom_point(aes(x=tissue, y=TPM, color=gene_id))+
  ggtitle("GRF Expression Across Tissues")

ggplot(consensus)+
  geom_point(aes(x=tissue, y=FPKM, color=gene_id))+
  ggtitle("GRF Expression Across Tissues")

ggplot(consensus)+
  geom_boxplot(aes(x=gene_id, y=TPM, fill=gene_id))+
  ggtitle("Average Expression Levels of GRFs in Total")

ggplot(consensus)+
  geom_boxplot(aes(x=gene_id, y=FPKM, fill=gene_id))+
  ggtitle("Average Expression Levels of GRFs in Total")


# Errored Pipelines

## Original Reference Files

> Arabidopsis_thaliana.TAIR10.37.gtf was downloaded from the Ensembl FTP via Cyberduck.

wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

wget https://www.arabidopsis.org/download_files/Sequences/Araport11_blastsets/Araport11_genes.201606.cdna.fasta.gz

> There appeared to be some issues with formatting of Chromosome 1, but could have been user/syntax error.

## Failed Mapping and Counting

> NOTE: Though there is overlap in file/directory names, each command was run in a unique directory.

bowtie2-build Araport11_genes.201606.cdna.fasta TAIR10

bowtie2-build GCF_000001735.4_TAIR10.1_cds_from_genomic.fna TAIR10.1

> bowtie2 outputs gave a large difficulty in running through RSEM for counting.

STAR --runMode genomeGenerate --genomeDir TAIR10 --genomeFastaFiles TAIR10_chr_all.fas --sjdbGTFfile Arabidopsis_thaliana.TAIR10.37.gtf

> Once this was counted began, this reference flagged at error due to some formatting issues with the first chromosome, according to the error messages.

STAR --runMode genomeGenerate --genomeDir TAIR10.1 --genomeFastaFiles GCF_000001735.4_TAIR10.1_genomic.fna --sjdbGTFfile GCF_000001735.4_TAIR10.1_genomic.gtf

> This worked properly, but the .gtf needed a few genes removed due to incorrect orientation data.

> There was an attempt at using salmon, but the permissions on colossus would not allow it to run and the software was too old for the security on my personal computer to allow the package to download.

STAR --genomeDir ../A.thaliana_Genome/TAIR10.1 --readFilesIn <Read1.fastq> <Read2.fastq> --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix <Sample_ID>

> This STAR mapping would run successfully for each sample, but the output would give formatting issues to RSEM, failing to count with the command below.

rsem-calculate-expression --bam --no-bam-output --paired-end --forward-prob 0 <Sample_ID>Aligned.toTranscriptome.out.bam rsem/TAIR10.1 <Sample_ID>





