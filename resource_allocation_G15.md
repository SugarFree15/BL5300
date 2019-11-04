# Downloading the Raw Reads

gdown.pl https://drive.google.com/file/d/1AAi_-gKOmu9lEOgTJTrAtg4pHHqxessZ/edit G15_S3_L001_R1_001.fast.gz

gdown.pl https://drive.google.com/file/d/1p5ht1Tm9ZuJmqMroZWh3MbjM0qCpTAuo/edit G15_S3_L001_R1_001.fast.gz

> These commands pull the zipped .fastq files for the forward and reverse reads of this sequencing run off of our Google Drive and onto Colossus.

# Quality Filtering

fastqc G15_S3_L001_R1_001.fastq.gz

fastqc G15_S3_L001_R2_001.fastq.gz

> These run the untrimmed reads through FastQC to assess their initial quality.

cutadapt -q 20,20 -m 50 --max-n 0 -o G15_S3_L001_R1_001.cutadapt.fastq -p G15_S3_L001_R2_001.cutadapt.fastq G15_S3_L001_R1_001.fastq.gz G15_S3_L001_R2_001.fastq.gz

> This trims the raw reads from both runs according to a minimum Fredd score of 20, minimum final length of 50 bp, and no "N" calls allowed, depositing the new sequences in new .fastq files.

fastqc G15_S3_L001_R1_001.cutadapt.fastq

fastqc G15_S3_L001_R2_001.cutadapt.fastq

> These verify the improved quality of the trimmed reads using FastQC.

# De Novo Assembly

conda activate de_novo

spades.py -k 21,51,71,91,111,127 --careful --pe1-1 G15_S3_L001_R1_001.cutadapt.fastq --pe1-2 G15_S3_L001_R2_001.cutadapt.fastq -o G15_spades_output

> These first create a region we can perform SPAdes in, then runs SPAdes to make de bruijn graphs at various k-mer lengths (21, 51, 71, 91, 111, and 127) and create a concensus de novo assembly in a new directory.

# Assembly Quality Assessment

quast G15_spades_output/contigs.fasta

> This runs our concensus contigs through quast to assess their quality and output an .html file of the report. 

# Genome Annotation

conda activate annotation

> This activates my annotation environment containing prokka.

mkdir annotation

cp G15_spades_output/contigs.fasta annotation

> These make a directory for my annotation data and then copies my assembled contigs into this new directory.

awk '/^>/{print ">G15_" ++i; next}{print}' < contigs.fasta > contigs_names.fasta

> Using awk, the names of the contigs were changed to simply a G15 prefix followed by an ID number.

prokka --outdir G15 --prefix G15 contigs_names.fasta

> With the renamed contigs, I annotated them using the prokka pipeline, outputting the results into the directory G15.

# Genome Analysis

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < G15.ffn > G15.fasta

> This takes the .ffn file and makes all the sequences only on single lines in .fasta format.

grep -A1 "resistance" G15.fasta > G15_resistance.fasta

> This created a .fasta file of all nucleotide sequences of all 24 "resistance" proteins found in the organism.

grep -A1 "transporter" G15.fasta > G15_transporter.fasta

> This created a .fasta file of all nucleotide sequences of all 82 "transporter" proteins found in the organism.

grep -A1 "efflux" G15.fasta > G15_efflux.fasta

> This created a .fasta file of all nucleotide sequences of all 15 "efflux" proteins found in the organism.
