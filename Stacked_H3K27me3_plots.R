
# Generating stacked H3K27me3 plots

# HPC = Bunya
# wd = "/scratch/user/uqesinni/TRIAGEi_analysis"
# RDM: /QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/

#################################################################################################################
################################### Stacked H3K27me3 plot across 833 biosamples #################################
#################################################################################################################

# 1) Get Epimap file with merged H3K27me3 peak data from all 833 biosamples

# scp uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q4879/Epimap/H3K27me3/epimap_H3K27me3_merged_hg19.bed /scratch/user/uqesinni/TRIAGEi_analysis/

library(data.table)
library(dplyr)
library(ggplot2)

epimap <- fread('/QRISdata/Q4879/Epimap/H3K27me3/epimap_H3K27me3_merged_hg19.bed')

#checks
head(epimap)
#      V1     V2     V3     V4    V5       V6
# 1: chr1  10175  10475 peak_1 0.001 BSS00215
# 2: chr1  10175  10475 peak_1 0.001 BSS00216
# 3: chr1  10175  10475 peak_1 0.001 BSS01446
# 4: chr1  87050  87300 peak_1 0.000 BSS00343
# 5: chr1  87075  87300 peak_1 0.000 BSS00341
# 6: chr1 795025 802600 peak_2 0.029 BSS00343
length(unique(epimap$V6))
# [1] 833
nrow(epimap)
# [1] 67907856
summary(epimap$V5)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00200 0.00500 0.01271 0.01300 1.00000

# figuring out what V5 is --> looks like it's normalized domain breadth within each unique biosample
# subs <- epimap[epimap$V6 == 'BSS00215', ]
# subs <- subs[order(-subs$V5), ]
# subs <- subs[order(-subs$V7), ]
# subs$V8 <- ifelse(subs$V5 >= quantile(subs$V5, 0.95), "Y", "N")

#### Create new column with top 5% broadest domains marked within each biosample
epimap1 <- epimap %>%
  group_by(V6) %>%
  mutate(V7 = ifelse(V5 >= quantile(V5, 0.95), "Y", "N"))

epimap1 <- as.data.table(epimap1)

write.table(epimap1, 'epimap_H3K27me3_merged_broadannot_hg19.bed', quote = F, col.names = F, row.names = F, sep = '\t') 

##############################################################################################################################

# 2) Get Epimap biosample metadata for ALL biosamples

# wget https://personal.broadinstitute.org/cboix/epimap/metadata/main_metadata_table.tsv

library(data.table)
library(ggplot2)

meta <- fread('main_metadata_table.tsv')

#checks
length(unique(meta$id))
# [1] 916
length(unique(meta$GROUP))
# [1] 33
table(meta$GROUP)

#        Adipose Blood & T-cell           Bone          Brain         Cancer 
#             13             41              5             54            122 
#      Digestive      Endocrine    Endothelial     Epithelial       ES-deriv 
#             83             29             11             57             18 
#            ESC            Eye          Heart   HSC & B-cell           iPSC 
#              9              5             43             44             15 
#         Kidney          Liver           Lung Lymphoblastoid        Mesench 
#             62             12             45             35              4 
#         Muscle         Myosat       Neurosph          Other       Pancreas 
#             61              6              4              8              8 
# Placenta & EEM            PNS   Reproductive     Sm. Muscle         Spleen 
#             31             10              7              8              9 
#        Stromal         Thymus        Urinary 
#             42             11              4 


# Merge H3K27me3 833 biosamples with metadata
epimap <- fread('epimap_H3K27me3_merged_broadannot_hg19.bed')

all_biosamples <- as.data.table(unique(epimap$V6))
meta_merge <- merge(meta, all_biosamples, by.x = "id", by.y = "V1")

write.table(meta_merge, 'Epimap_833_biosamples_metadata.bed', quote = F, col.names = T, row.names = F, sep = '\t') 

##############################################################################################################################

# 3) Bedtools intersect to get H3K27me3 peaks that overlap with locus of interest

# MYH7
# chr14:23,881,947-23,904,870

bedtools intersect -a epimap_H3K27me3_merged_broadannot_hg19.bed -b <(echo -e "chr14\t23881947\t23904870") -u > overlapping_peaks_MYH7.bed

##############################################################################################################################

# 4) Add missing biosample information and metadata

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)

epimap <- fread('epimap_H3K27me3_merged_broadannot_hg19.bed')

loi <- fread('overlapping_peaks_u.bed') #locus of interest
# loi <- fread('overlapping_peaks_MYOD1.bed') #locus of interest
# loi <- fread('overlapping_peaks_GAPDH.bed') #locus of interest
# loi <- fread('overlapping_peaks_MYH2.bed') #locus of interest
# loi <- fread('overlapping_peaks_H19.bed') #locus of interest


# Get domain breadth column and reorder
loi$V8 <- loi$V3-loi$V2
loi <- loi[order(-loi$V8), ]

# Get missing biosamples added on
all_biosamples <- unique(epimap$V6)
loi_biosamples <- unique(loi$V6)
missing_bs <- setdiff(all_biosamples, loi_biosamples)
length(missing_bs)
# [1] 34

new_rows <- data.frame(
  V1 = "chr14",
  V2 = NA_real_,
  V3 = NA_real_,
  V4 = NA_character_,
  V5 = NA_real_,
  V6 = missing_bs,
  V7 = NA_character_,
  V8 = NA_integer_,
  stringsAsFactors = FALSE
)

# Add on missing biosamples to loi dataframe
loi <- rbind(loi, new_rows)

length(unique(loi$V6))
# [1] 833

# Get biosample metadata (celltype group) added on

meta <- fread('Epimap_833_biosamples_metadata.bed')

# Merge to get metadata info
loi <- merge(meta[, c("id","GROUP")], loi, by.x = "id", by.y = "V6")

length(unique(loi$GROUP))
# [1] 33


##############################################################################################################################

# 5) Generate stacked H3K27me3 plot for locus of interest

# MYH7
plot <- ggplot(loi, aes(x = V2, xend = V3, y = reorder(id, V8), yend = reorder(id, V8), color = V7)) +
  geom_segment(size = 0.5) +  
  scale_color_manual(values = c("N" = "#2377B9", "Y" = "#E86666"), na.value = "white") +
  labs(x = "Genomic Range", y = "Biosample") +
  theme_minimal()
ggsave("240101_stacked_H3K27me3_MYH7.pdf", plot, width = 8, height = 6)


##############################################################################################################################

# 6) Biosample distribution plot

table(loi$V7, loi$GROUP)

tbl <- as.data.frame(table(loi$V7, loi$GROUP))

tbl$Proportion <- with(tbl, Freq / ave(Freq, Var2, FUN = sum))

# Stacked prop plot
plot <- ggplot(tbl, aes(fill=Var1, y=Proportion, x=Var2)) + 
    geom_bar(position="fill", stat="identity") +
    labs(x = "Biosample", y = "Proportion") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c('N' = 'darkgrey', 'Y' = 'darkred'))
ggsave("240101_stacked_props_MYH7.pdf", plot, width = 8, height = 6)


# Boxplot
plot <- ggplot(loi, aes(x = GROUP, y = V8, fill = GROUP)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "GROUP", y = "V8") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("240101_stacked_boxplot_MYH7.pdf", plot, width = 8, height = 6)


# Violin plot
plot <- ggplot(loi, aes(x = GROUP, y = V8, fill = GROUP)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, show.legend = FALSE) +
  labs(x = "GROUP", y = "V8") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("240101_stacked_violinplot_MYH7.pdf", plot, width = 8, height = 6)


#################################################################################################################
######################################### Generating Figure Panels ############################################
#################################################################################################################

# Generate bed files for example loci

# GAPDH
# chr12:6643585-6647537
bedtools intersect -a epimap_H3K27me3_merged_broadannot_hg19.bed -b <(echo -e "chr12\t6643585\t6647537") -u > overlapping_peaks_GAPDH.bed

# MYH2
# chr17:10424465-10452940
bedtools intersect -a epimap_H3K27me3_merged_broadannot_hg19.bed -b <(echo -e "chr17\t10424465\t10452940") -u > overlapping_peaks_MYH2.bed

# MYOD1
# chr11:17741110-17743678
bedtools intersect -a epimap_H3K27me3_merged_broadannot_hg19.bed -b <(echo -e "chr11\t17741110\t17743678") -u > overlapping_peaks_MYOD1.bed

# H19
# chr11:2016406-2019065
bedtools intersect -a epimap_H3K27me3_merged_broadannot_hg19.bed -b <(echo -e "chr11\t2016406\t2019065") -u > overlapping_peaks_H19.bed


##############################################################################################################################

# Generate stacked H3K27me3 plots for example loci

#MYOD1
plot <- ggplot(loi, aes(x = V2, xend = V3, y = reorder(id, V8), yend = reorder(id, V8), color = V7)) +
  geom_segment(size = 0.5) +  
  scale_color_manual(values = c("N" = "#2377B9", "Y" = "#E86666"), na.value = "white") +
  geom_vline(xintercept = c(17705000, 17805000, 17741110, 17743678), linetype = "dashed", color = "black", size = 0.5) + #add lines to mark gene location
  labs(x = "Genomic Range", y = "Biosample") +
  theme_minimal()
ggsave("240101_stacked_H3K27me3_MYOD1.pdf", plot, width = 8, height = 6)


# GAPDH
plot <- ggplot(loi, aes(x = V2, xend = V3, y = reorder(id, V8), yend = reorder(id, V8), color = V7)) +
  geom_segment(size = 0.5) +  
  scale_color_manual(values = c("N" = "#2377B9", "Y" = "#E86666"), na.value = "white") +
  geom_vline(xintercept = c(6643200, 6644900, 6643585, 6647537), linetype = "dashed", color = "black", size = 0.5) + #add lines to mark gene location
  labs(x = "Genomic Range", y = "Biosample") +
  theme_minimal()
ggsave("240101_stacked_H3K27me3_GAPDH.pdf", plot, width = 8, height = 6)


# H19
plot <- ggplot(loi, aes(x = V2, xend = V3, y = reorder(id, V8), yend = reorder(id, V8), color = V7)) +
  geom_segment(size = 0.5) +  
  scale_color_manual(values = c("N" = "#2377B9", "Y" = "#E86666"), na.value = "white") +
  geom_vline(xintercept = c(1990000, 2095000, 2016406, 2019065), linetype = "dashed", color = "black", size = 0.5) + #add lines to mark gene location
  labs(x = "Genomic Range", y = "Biosample") +
  theme_minimal()
ggsave("240101_stacked_H3K27me3_H19.pdf", plot, width = 8, height = 6)


# MYH2
plot <- ggplot(loi, aes(x = V2, xend = V3, y = reorder(id, V8), yend = reorder(id, V8), color = V7)) +
  geom_segment(size = 0.5) +  
  scale_color_manual(values = c("N" = "#2377B9", "Y" = "#E86666"), na.value = "white") +
  geom_vline(xintercept = c(10400000, 10455000, 10424465, 10452940), linetype = "dashed", color = "black", size = 0.5) + #add lines to mark gene location
  labs(x = "Genomic Range", y = "Biosample") +
  theme_minimal()
ggsave("240101_stacked_H3K27me3_MYH2.pdf", plot, width = 8, height = 6)


##############################################################################################################################
##############################################################################################################################
# scp files over from bunya scratch to RDM
# scp -r /scratch/user/uqesinni/TRIAGEi_analysis/*.bed uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/

# scp -r /scratch/user/uqesinni/TRIAGEi_analysis/*.pdf uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/

# scp -r /scratch/user/uqesinni/TRIAGEi_analysis/*.txt uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/

# # transfer to local
# scp -r uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/240101* .



