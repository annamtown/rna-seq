## This is a shell script for running RNA-seq data from Listeria monocytogenes:
## D-glucose sample: https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6281667
## D-allose sample: https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6281666
## reference genome: https://www.ebi.ac.uk/ena/data/view/AL591824
## For BINF8441: Bioinformatic Statistics Spring 2018 Final Project
## Anna Townsend
## May 1 2018

# reference genome fasta file was downloaded from ENA and copied using scp onto server to run this script
# name of reference genome: AL591824.fasta

# Download sequence reads
fastq-dump SRR6281667    #EGD-e with D-glucose

fastq-dump SRR6281666    #EGD-e with D-allose

# Index reference genome
# -f indicates fasta format
# -c indicates reference sequence given on cmd line (as <seq_in>)
bowtie2-build -f AL591824.fasta ref_seq

# run Bowtie2 to align RNA-seq files to indexed reference file
bowtie2 -x ref_seq -U SRR6281667.fastq -S SRR6281667.sam

#anna@GuanYin:~/rna-seq$ bowtie2 -x ref_seq -U SRR6281667.fastq -S SRR6281667.sam
#25052599 reads; of these:
#  25052599 (100.00%) were unpaired; of these:
#    481129 (1.92%) aligned 0 times
#    23871576 (95.29%) aligned exactly 1 time
#    699894 (2.79%) aligned >1 times
#98.08% overall alignment rate

bowtie2 -x ref_seq -U SRR6281666.fastq -S SRR6281666.sam

#anna@GuanYin:~/rna-seq$ bowtie2 -x ref_seq -U SRR6281666.fastq -S SRR6281666.sam
#25222151 reads; of these:
#  25222151 (100.00%) were unpaired; of these:
#    719810 (2.85%) aligned 0 times
#    24167759 (95.82%) aligned exactly 1 time
#    334582 (1.33%) aligned >1 times
#97.15% overall alignment rate

# Use Tophat for creation of accepted_hits.bam file
# for SRR6281666
tophat ref_seq SRR6281666.fasta

# for SRR6281667
tophat -o tophat_out_1 ref_seq SRR6281667.fasta


# run cufflinks for isoform-based summarization
# for SRR6281666
cufflinks /home/anna/rna-seq/tophat_out/accepted_hits.bam

# for SRR6281667
cufflinks /home/anna/rna-seq/tophat_out_1/accepted_hits.bam

# Create text file containing paths to transcript.gtf files for each sample
# nano assembly_GTF_list.txt

# merge transcript.gtf files from each lineage
cuffmerge assembly_GTF_list.txt

# discover novel transcript
# cuffcompare /home/anna/rna-seq/merged_asm/merged.gtf

# compare gene expressions of two samples
cuffdiff merged.gtf /home/anna/rna-seq/tophat_out/accepted_hits.bam /home/anna/rna-seq/tophat_out_1/accepted_hits.bam




################################################################################
# Commands I tried to use, but did not run successfully
################################################################################
# Convert sam file to bam file using Samtools view
#samtools view -S -b SRR6281667.sam > SRR6281667.bam

#samtools view -S -b SRR6281666.sam > SRR6281666.bam

# Fix header in bam file
# samtools view -bh SRR6281667.bam | samtools bam2fq -

# sort the bam file
# samtools sort -@ 4 -O bam -T SRR6281667.bam -o SRR6281667.sort.bam

# samtools sort -@ 4 -O bam -T SRR6281666.bam -o SRR6281666.sort.bam

# convert bam file to bed file to use in DEGseq
#bedtools bamtobed -i SRR6281667.bam > SRR6281667.bed

#bedtools bamtobed -i SRR6281666.bam > SRR6281666.bed

# Gene annotation of reference sequence using Prokka
#prokka --prefix refgenome AL591824.fasta

# the following commmands are to use in R
## To install DEGseq package:
#> source("https://bioconductor.org/biocLite.R")
#> biocLite("DEGseq")
#> library(DEGseq)

#> glucose_667 <- system.file("extdata", "SRR6281667.bed.txt", package="DEGseq")
#> allose_666 <- system.file("extdata", "SRR6281666.bed.txt", package="DEGseq")
#> refFlat <- system.file("extdata", "refgenome.gff.txt", package="DEGseq")
#> mapResultBatch1 <- c(glucose_667) ## only use the data from kidneyR1L1 and liverR1L2
#> mapResultBatch2 <- c(allose_666)
#> outputDir <- file.path(tempdir(), "/Users/tennisluver/Desktop")
#> DEGseq(mapResultBatch1, mapResultBatch2, fileFormat="bed", refFlat=refFlat, outputDir=outputDir, method="LRT")
