# Downloading RNA-Seq Data from NCBI Bioprojects

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

# Downloading A. thaliana Genome and .gtf Files

wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < TAIR10_chr_all.fas > TAIR10_chr_all_tmp.fas
tr "\t" "\n" < TAIR10_chr_all_tmp.fas > TAIR10_chr_all.fasta

> Downloaded the TAIR10 Arabidopsis genome and spliced the sequences of the .fasta into a single line each with awk and tr, making a properly formatted .fasta of chromosomes and organelles.

wget -c ftp://ftp.ensemblgenomes.org/pub/release-37/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.37.gtf.gz

> Downloaded the gene information .gtf file for TAIR10.

# Mapping Raw RNA-Seq Reads to the Genome

STAR --runMode genomeGenerate --genomeDir A.thaliana_Genome --genomeFastaFiles TAIR10_chr_all.fasta --sjdbGTFfile Arabidopsis_thaliana.TAIR10.37.gtf

> Indexed the TAIR10 Assembly for use with the STAR aligner.

STAR --genomeDir A.thaliana_Genome/ --readFilesIn ERR3001915_1.fastq ERR3001915_2.fastq --outSAMtype SAM --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --outFileNamePrefix "ERR3001915"


> Mapped each of my datasets to the TAIR10 indexed genome with outputs into .sam files.
