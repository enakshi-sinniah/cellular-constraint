
# Permutation test using regioneR

# HPC = Bunya
# wd = "/scratch/user/uqesinni/TRIAGEi_analysis"
# RDM: /QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/

#################################################################################################################
######################## RegioneR Permutation Tests splitting coding vs. non-coding #############################
#################################################################################################################

# Generate bed file for whole genome

# wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes

hg38 <- fread('hg38.chrom.sizes')
hg38$V3 <- hg38$V2
hg38$V2 <- "1"
write.table(hg38, 'Whole_genome_hs_hg38.bed', sep = '\t', quote = F, col.names = F, row.names = F)

###################################################################

# 1) Get bed files with hg19 coding and non-coding sequence split

# scp uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/NCBI_RefSeqselect_singlegene_hg19_input.bedgraph /scratch/user/uqesinni/TRIAGEi_analysis/ 
# scp uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/Whole_genome_hs_hg19.bed /scratch/user/uqesinni/TRIAGEi_analysis/ 
# scp uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/phastcons_top0.9_hg19.bed /scratch/user/uqesinni/TRIAGEi_analysis/ 
# scp uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/triagebp_P_hg19.bed /scratch/user/uqesinni/TRIAGEi_analysis/ 

# Use bedtools to generate required input files

bedtools subtract -a Whole_genome_hs_hg19.bed -b NCBI_RefSeqselect_singlegene_hg19_input.bedgraph > Noncoding_genome_hs_hg19.bed

bedtools intersect -a triagebp_P_hg19.bed -b Noncoding_genome_hs_hg19.bed > triagebp_P_noncoding_hg19.bed

bedtools intersect -a triagebp_P_hg19.bed -b Noncoding_genome_hs_hg19.bed -v > triagebp_P_coding_hg19.bed


# Check on generated files
library(regioneR)
library(data.table)
library(ggplot2)
library(BSgenome)
installed.genomes()

pcgs <- fread('NCBI_RefSeqselect_singlegene_hg19_input.bedgraph') #21436 NM_ transcripts

hg19 <- fread('Whole_genome_hs_hg19.bed')
noncoding <- fread('Noncoding_genome_hs_hg19.bed')

phastcons <- fread('phastcons_top0.9_hg19.bed')

all_triage <- fread('triagebp_P_hg19.bed')
coding_triage <- fread('triagebp_P_coding_hg19.bed')
noncoding_triage<- fread('triagebp_P_noncoding_hg19.bed')


################################################################################################################

# 2) Run overlapPermTest on phastcons top >0.9

# Submit slurm job on Bunya
# Rscript
# 240103_Permute_1000_phastcons_noncoding.R

library(regioneR)
library(data.table)
library(ggplot2)
library(BSgenome)

file1 <- read.table("triagebp_P_noncoding_hg19.bed", header = FALSE, stringsAsFactors = FALSE)
file2 <- read.table("phastcons_top0.9_hg19.bed", header = FALSE, stringsAsFactors = FALSE)

# Convert the BED data to GRanges objects
gr1 <- toGRanges(read.table("triagebp_P_noncoding_hg19.bed",sep="\t"))
gr2 <- toGRanges(read.table("phastcons_top0.9_hg19.bed",sep="\t"))

# numOverlaps(gr1, gr2, count.once=TRUE)
# [1] 62773

result <- overlapPermTest(A = gr1, B = gr2, ntimes = 1000, mc.cores=20, genome="BSgenome.Hsapiens.UCSC.hg19", verbose=TRUE, count.once=TRUE)

result 

pdf(file="240103_phastcons0.9_triagebp_noncoding_1000_perm.pdf")
plot(result)
dev.off()

p_val <- result$pvalue

p_val

result <- write.table(result, 'phastconsTRIAGEbp_noncoding_permute100.txt')

# scp  uqesinni@bunya.rcc.uq.edu.au:/scratch/user/uqesinni/TRIAGEi_analysis/240103_phastcons0.9_triagebp_noncoding_1000_perm.pdf .


# Job IDs:
# 240103_Permute_1000_phastcons_noncoding.SH
# 240103_Permute_1000_phastcons_noncoding.R
# Took ~1hr to run using select=1, ncpus=20

#################################################################################
# Slurm batch script
#240103_Permute_1000_phastcons_noncoding.SH

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --job-name=Permute_1000
#SBATCH --time=3:00:00
#SBATCH --partition=general
#SBATCH --account=a_palpant

MY_DIR=/scratch/user/uqesinni/TRIAGEi_analysis
cd $MY_DIR

module load miniconda3/4.9.2

conda activate r_env

/home/uqesinni/.conda/envs/r_env/bin/Rscript 240103_Permute_1000_phastcons_noncoding.R


# sbatch 240103_Permute_1000_phastcons_noncoding.SH
# squeue -u uqesinni

#####################################################################################################
# Run overlap permutation test on zooHARs and hsSVs 

# zooHARs and hsSVs downloaded manually from Keough et al. (2023) paper: https://www.science.org/doi/full/10.1126/science.abm1696
# Supplementary table 1 & 3

# scp uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/phastCons_100way_hg19.bedgraph /scratch/user/uqesinni/TRIAGEi_analysis/

library(regioneR)
library(data.table)
library(ggplot2)
library(BSgenome)
installed.genomes()

hars <- fread('zooHARs_hg38.txt')
dim(hars)
# [1] 312   8
summary(hars$end - hars$start)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 50.00   76.75  117.50  157.61  187.25  789.00 

hsSVs <- fread('hsSVs_hg38.txt')
dim(hsSVs)
# [1] 17789     7
summary(hsSVs$end - hsSVs$start)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   50     294     324    1043    1015   64998

hars_bed <- subset(hars, select=c(chrom, start, end, simple_name))
write.table(hars_bed, 'zooHARs_hg38_TRIAGEinput.bed', quote = F, col.names = F, row.names = F, sep = '\t') 

hsSVs_bed <- subset(hsSVs, select=c(chrom, start, end, name))
write.table(hsSVs_bed, 'hsSVs_hg38_TRIAGEinput.bed', quote = F, col.names = F, row.names = F, sep = '\t') 

###############################
# Get hg38 P TRIAGEbp file
library(inflection)

x <- fread('TRIAGEi_ReferenceH3K27me3_epimap_hg38.bedgraph')
x <- x[order(-x$V4),]
xy <- cbind(x, "Rank"=1:nrow(x)) 
ede(xy$Rank, xy$V4,index=0)
#         j1 j2 chi
# EDE 727985  1 NaN

z <- xy[1:727985,]
z <- subset(z, select=c(V1, V2, V3, V4))
write.table(z, 'triagebp_P_hg38.bed', quote = F, col.names = F, row.names = F, sep = '\t') 

#################################################################################################

# Run overlap permutation test for zooHARs vs TRIAGEbp P 

gr1 <- toGRanges(read.table("zooHARs_hg38_TRIAGEinput.bed",sep="\t"))
gr2 <- toGRanges(read.table("triagebp_P_hg38.bed",sep="\t"))

# run 10,000 permutations
result <- overlapPermTest(A = gr1, B = gr2, ntimes = 10000, mc.cores=20, genome="BSgenome.Hsapiens.UCSC.hg38", verbose=TRUE, count.once=TRUE)

result 
# P-value: 0.0001999800019998
# Z-score: 5.4543
# Number of iterations: 10000
# Alternative: greater
# Evaluation of the original region set: 15
# Evaluation function: numOverlaps
# Randomization function: randomizeRegions


pdf(file="231017_zooHAR_triagebp_10000_perm.pdf")
plot(result)
dev.off()

## Bunya slurm job details. NOTE: ran in less than 10mins

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --job-name=Permute_1000
#SBATCH --time=1:00:00
#SBATCH --partition=general
#SBATCH --account=a_palpant

MY_DIR=/scratch/user/uqesinni/TRIAGEi_analysis
cd $MY_DIR

module load miniconda3/4.9.2

conda activate r_env

/home/uqesinni/.conda/envs/r_env/bin/Rscript 231017_Permute_1000_zooHAR.R

# sbatch 231017_Permute_1000_zooHAR.SH
# squeue -u uqesinni


#################################################################################################

# Run overlap permutation test for hsSVs vs TRIAGEbp P 

gr1 <- toGRanges(read.table("hsSVs_hg38_TRIAGEinput.bed",sep="\t"))
gr2 <- toGRanges(read.table("triagebp_P_hg38.bed",sep="\t"))

# run 10,000 permutations
result <- overlapPermTest(A = gr1, B = gr2, ntimes = 10000, mc.cores=20, genome="BSgenome.Hsapiens.UCSC.hg38", verbose=TRUE, count.once=TRUE)

result 
# P-value: 9.99900009999e-05
# Z-score: -7.8533
# Number of iterations: 10000
# Alternative: less
# Evaluation of the original region set: 124
# Evaluation function: numOverlaps
# Randomization function: randomizeRegions


pdf(file="231018_hsSVs_triagebp_10000_perm.pdf")
plot(result)
dev.off()

## Bunya slurm job details. NOTE: ran in less than 10mins

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --job-name=Permute_1000
#SBATCH --time=1:00:00
#SBATCH --partition=general
#SBATCH --account=a_palpant

MY_DIR=/scratch/user/uqesinni/TRIAGEi_analysis
cd $MY_DIR

module load miniconda3/4.9.2

conda activate r_env

/home/uqesinni/.conda/envs/r_env/bin/Rscript 231018_Permute_1000_hsSVs.R

# sbatch 231018_Permute_1000_hsSVs.R
# squeue -u uqesinni









