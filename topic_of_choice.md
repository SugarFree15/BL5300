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

> Using the fasterq-dump function in the NCBI SRA-toolkit, I downloaded each replicate of each tissues from a variety of studies in the NCBI Bioproject database, included paired runs, if applicable.

## Downloading A. thaliana Genome and Info Files

> Accessed the TAIR10.1 Assembly of the A. thaliana genome through the RefSeq Genome FTP via Cyberduck, coping the genomic (GCF_000001735.4_TAIR10.1_genomic.fna) and info (GCF_000001735.4_TAIR10.1_genomic.gff) files to colossus in my ~/scratch/BL5300/topic_of_choice/Athal_ncbi directory.

## Mapping Raw RNA-Seq Reads to the Genome

cd ~/scratch/BL5300/topic_of_choice
mkdir star-rsem
cd star-rsem
STAR --runMode genomeGenerate --genomeDir ../Athal_ncbi/ --genomeFastaFiles GCF_000001735.4_TAIR10.1_genomic.fna --sjdbGTFfile GCF_000001735.4_TAIR10.1_genomic.gff

> Indexed the TAIR10 Assembly for use with the STAR aligner.

STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../ERR3001915_1.fastq ../ERR3001915_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "ERR3001915"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../ERR3001916_1.fastq ../ERR3001916_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "ERR3001916"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../SRR10054928_1.fastq ../SRR10054928_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "SRR10054928"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../SRR10054929_1.fastq ../SRR10054929_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "SRR10054929"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../SRR10054931_1.fastq ../SRR10054931_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "SRR10054931"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../SRR10054932_1.fastq ../SRR10054932_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "SRR10054932"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../SRR10320043_1.fastq ../SRR10320043_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "SRR10320043"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../SRR10320044_1.fastq ../SRR10320044_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "SRR10320044"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../SRR10320045_1.fastq ../SRR10320045_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "SRR10320045"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../ERR1727061_1.fastq ../ERR1727061_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "ERR1727061"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../ERR1727060_1.fastq ../ERR1727060_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "ERR1727060"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../ERR1727059_1.fastq ../ERR1727059_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "ERR1727059"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../SRR1610505_1.fastq ../SRR1610505_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "SRR1610505"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../SRR1610506_1.fastq ../SRR1610506_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "SRR1610506"
STAR --genomeDir ../Athal_ncbi/ --readFilesIn ../SRR1610507_1.fastq ../SRR1610507_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "SRR1610507"


> Mapped each of my datasets to the TAIR10 indexed genome in preparation for RSEM via STAR (WARNING: this step takes excessively long).

## Counting Transcript Reads


# Trial Run

> Started by downloading the GCF_000001735.4_TAIR10.1_genomic.fna A. thaliana genome file onto my computer, along with its cooresponding GFF file, and used cyberduck to place it on colossus.

bowtie2-build ../Athal_ncbi/GCF_000001735.4_TAIR10.1_genomic.fna Athal_gen

