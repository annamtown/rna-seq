## This is a shell script for running RNA-seq data from Listeria monocytogenes:
## D-glucose sample: https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6281667
## D-allose sample: https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6281666
## reference genome: https://www.ebi.ac.uk/ena/data/view/AL591824
## For BINF8441: Bioinformatic Statistics Spring 2018 Final Project
## Anna Townsend
## May 1 2018

# reference genome fasta file was downloaded from ENA and copied using scp into directory 

# Download sequence reads
fastq-dump --fasta SRR6281667    #EGD-e with D-glucose

fastq-dump --fasta SRR6281666    #EGD-e with D-allose

# run Bowtie2
# bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
# <bt2-idx>  Index filename prefix (minus trailing .X.bt2)
# <m1>       Files with #1 mates, paired with files in <m2>
#           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
# <m2>       Files with #2 mates, paired with files in <m1>
#           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
# <r>        Files with unpaired reads
#           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
# <sam>      File for SAM output (default: stdout)
#
# <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
# specified many times  E.g. '-U file1.fq,file2.fq -U file3.fq'
