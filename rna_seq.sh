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
