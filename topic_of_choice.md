# Finalized Pipline

## Downloading RNA-Seq Data from NCBI Bioprojects

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

> Accessed the TAIR10.1 Assembly of the A. thaliana genome through the RefSeq Genome FTP via Cyberduck, coping the genomic (GCF_000001735.4_TAIR10.1_genomic.fna) and info (GCF_000001735.4_TAIR10.1_genomic.gff) files to colossus in my ~/scratch/BL5300/topic_of_choice/Athal_ncbi directory.

## Mapping Raw RNA-Seq Reads to the Genome and Counting Transcripts

> NOTE: I was slightly desperate to troubleshoot my syntaxes, so the working directory and some files are named with 'please' with no significance other than hoping this pipeline worked.

cd ~/scratch/BL5300/topic_of_choice/please

grep -v "unknown_transcript_1" GCF_000001735.4_TAIR10.1_genomic.gtf > please_genomic.gtf

rsem-prepare-reference --gtf please_genomic.gtf --star ../Athal_ncbi/GCF_000001735.4_TAIR10.1_genomic.fna rsem/please_ncbi

> This indexed the TAIR10.1 genome in both RSEM and STAR at the same time using transcript info from a modified .gtf file that removed a set of unknown transcripts that would through errors from incorrect orientations in the files.

rsem-calculate-expression --no-bam-output --star ../SRR1610505.fastq rsem/please_ncbi SRR1610505

rsem-calculate-expression --no-bam-output --star ../.fastq rsem/please_ncbi 

rsem-calculate-expression --no-bam-output --star --paired-end --forward-prob 0 ../ERR3001915_1.fastq ../ERR3001915_2.fastq rsem/please_ncbi

rsem-calculate-expression --no-bam-output --star --paired-end --forward-prob 0 ../ERR1727060_1.fastq ../ERR1727060_2.fastq rsem/please_ncbi

rsem-calculate-expression --no-bam-output --star --paired-end --forward-prob 0 ../SRR10054931_1.fastq ../SRR10054931_2.fastq rsem/please_ncbi

> These took the raw reads and aligned them via STAR with parameters perfect for RSEM before then counting the transcripts and outputting the needed data with RSEM, taking into account whether or not each sample was a paired-end read.

# Trial Run

> Started by downloading the GCF_000001735.4_TAIR10.1_genomic.fna A. thaliana genome file onto my computer, along with its cooresponding GFF file, and used cyberduck to place it on colossus.

bowtie2-build ../Athal_ncbi/GCF_000001735.4_TAIR10.1_genomic.fna Athal_gen

