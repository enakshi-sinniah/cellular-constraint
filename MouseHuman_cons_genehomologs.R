# Human-Mouse CCS Conservation

# HPC = Bunya
# wd = "/scratch/user/uqesinni/TRIAGEi_analysis"
# RDM: /QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/

#################################################################################################################
########################################## Human-Mouse CCS Conservation #########################################
#################################################################################################################

# 1) Homologous genes between hg19 and mm10 

library(data.table)
library(dplyr)

homos <- fread("/QRISdata/Q4879/TRIAGEbp/homologs_hg19_mm10.txt")

#formatting
colnames(homos) <- c("Gene_id", "mm10_chr", "mm10_start", "mm10_end", "hg19_chr","hg19_start","hg19_end")
homos$mm10_chr <- paste("chr", homos$mm10_chr, sep="")
homos$hg19_chr <- paste("chr", homos$hg19_chr, sep="")
typeof(homos$mm10_end) #check if integer
# [1] "integer"

# remove duplicates
length(which(duplicated(homos$Gene_id)))
# [1] 4254
homos <- homos[!duplicated(homos$Gene_id), ]

# remove weird chr entries
homos <- homos[!(grepl("^chrGL|^chrMT", homos$hg19_chr))]
homos <- homos[!(grepl("^chrGL|^chrMT|^chrJH", homos$mm10_chr))]

# ignore case finding duplicates
homos$Gene_id <- tolower(homos$Gene_id)
homos <- homos[!duplicated(homos$Gene_id), ]

#order 
homos <- homos[order(homos$Gene_id), ]

#remove weird 'rik' genes
homos <- homos[!(grepl("rik$|rik2$", homos$Gene_id))]

dim(homos)
# [1] 21314     7

write.table(homos, 'homologous_genes_mm10hg19_meta_clean.txt', quote = F, col.names = T, row.names = F, sep = '\t')

# Generate BED file with human/mouse chr coordinates for homologous genes
hum_homos <- subset(homos, select=c(hg19_chr, hg19_start, hg19_end, Gene_id))
mou_homos <- subset(homos, select=c(mm10_chr, mm10_start, mm10_end, Gene_id))

write.table(hum_homos, 'homologous_genes_hg19_TRIAGEinput.bed', quote = F, col.names = F, row.names = F, sep = '\t') #no colnames
write.table(mou_homos, 'homologous_genes_mm10_TRIAGEinput.bed', quote = F, col.names = F, row.names = F, sep = '\t') 


# Get TRIAGEbp RTSs 
python TRIAGEi_main.py -i homologous_genes_hg19_TRIAGEinput.bed -o homologous_genes_hg19_rts.txt -r TRIAGEi_ReferenceH3K27me3_epimap_hg19.bedgraph

python TRIAGEi_main.py -i homologous_genes_mm10_TRIAGEinput.bed -o homologous_genes_mm10_rts.txt -r scaled_fold_change_RTS_mm10.bedgraph


################################

library(data.table)

homos <- fread("homologous_genes_mm10hg19_meta_clean.txt")
dim(homos)
# [1] 21314     7

hum_rts <- fread('homologous_genes_hg19_rts.txt')
dim(hum_rts)
# [1] 18728     2

hum_merged <- merge(hum_rts, homos, by.x = "V1" , by.y = "Gene_id", all.y = TRUE)
colnames(hum_merged)[colnames(hum_merged) %in% c("V1", "V2")] <- c("Gene_id", "human_RTS")
hum_merged$human_RTS <- ifelse(is.na(hum_merged$human_RTS), 0, hum_merged$human_RTS) #change NAs i.e. non-RTS regions to 0
summary(hum_merged$human_RTS)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.0000002 0.0001020 0.0096310 0.0024918 1.0000000

mou_rts <- fread('homologous_genes_mm10_rts.txt')
dim(mou_rts)
# [1] 20916     2
mou_merged <- merge(mou_rts, hum_merged, by.x="V1", by.y="Gene_id", all.y=TRUE)
colnames(mou_merged)[colnames(mou_merged) %in% c("V1", "V2")] <- c("Gene_id", "mouse_RTS")
mou_merged$human_RTS <- ifelse(is.na(mou_merged$human_RTS), 0, mou_merged$human_RTS) #change NAs i.e. non-RTS regions to 0
summary(mou_merged$mouse_RTS)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
 # 0.0000  0.0008  0.0018  0.0096  0.0059  0.9939     398

write.table(mou_merged, 'homologous_genes_mm10hg19_meta_cleaned.txt', quote = F, col.names = T, row.names = F, sep = '\t')


##########################################################################################################

# 2) Plot correlation 

library(ggplot2)
library(ggpubr)

comb <- fread('homologous_genes_mm10hg19_meta_cleaned.txt')

# remove weird gene names contributing to noise
comb <- comb[!(grepl("^gm|^bc|^klk|^al|^mir|^iq|^n-", comb$Gene_id))]

comb_corr <- cor.test(comb$mouse_RTS, comb$human_RTS, method = "spearman")
# 	Spearman's rank correlation rho

# data:  comb$mouse_RTS and comb$human_RTS
# S = 6.0496e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.678346  

ggscatter(comb, x = "mouse_RTS", y = "human_RTS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman")

# Density plot

library(dplyr)
library(viridis)
library(ggpointdensity)
library(ggpubr)
library(ggrastr)

ggplot(comb, aes(x = mouse_RTS, y = human_RTS)) +
  geom_pointdensity(size=0.7, adjust= 10) + #adjust value changes n_neighbours for density
  scale_color_viridis() +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +
  labs(x = "mouse RTS", y = "human RTS") +
  ggtitle(paste("Spearman Correlation =", round(comb_corr$estimate, 2), " (N =", length(comb$human_RTS), ")")) +
  theme_minimal()

# save rasterized pdf
p <- ggplot(comb, aes(x = mouse_RTS, y = human_RTS)) +
  geom_point(size = 1) +  
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +  
  labs(x = "mouse RTS", y = "human RTS") + 
  ggtitle(paste("Spearman Correlation =", round(comb_corr$estimate, 2), " (N =", length(comb$human_RTS), ")"))  

rastered_plot <- rasterize(p, width = 8, height = 6)
ggsave("240205_HomologousGenes_ALL_density_rts_correlation_cleaned.pdf", plot = p)

# scp uqesinni@bunya.rcc.uq.edu.au:/scratch/user/uqesinni/TRIAGEi_analysis/240205_HomologousGenes_ALL_density_rts_correlation_cleaned.pdf .

###############################################################################
# 3) Calculate inflection point 

library(inflection)

########
# Inflection calculation
comb$bin_human <- ifelse(comb$human_RTS > 0, 'Y', 'N')
comb$bin_mouse <- ifelse(comb$mouse_RTS > 0, 'Y', 'N')

#Human inflection
human <- comb[comb$bin_human == "Y", ]
human <- human[order(human$human_RTS, decreasing = T),] 
dim(human)
# [1] 16192    11
human$human_rank <- 1:16192 #formatting
#calculate inflection point
ede(human$human_rank, human$human_RTS, index=0) 
#      j1    j2  chi
# EDE 857 16192 8524.5
# Add priority info into homos metadata
comb <- comb[order(comb$human_RTS, decreasing = T),] 
comb$human_inflection <- ifelse(seq_len(nrow(comb)) <= 857, 'P', 'NP')

#Mouse inflection
mouse <- comb[comb$bin_mouse == "Y", ]
mouse <- mouse[order(mouse$mouse_RTS, decreasing = T),] 
dim(mouse)
# [1] 17556    13
mouse$mouse_rank <- 1:17556 #formatting
#calculate inflection point
ede(mouse$mouse_rank, mouse$mouse_RTS, index=0) 
#      j1    j2     chi
# EDE 743 17556 9149.5
# Add priority info into homos metadata
comb <- comb[order(comb$mouse_RTS, decreasing = T),] 
comb$mouse_inflection <- ifelse(seq_len(nrow(comb)) <= 743, 'P', 'NP')

#####################################################################################

# Correlation plot for Sig(P) in either species

library(ggplot2)
library(ggpubr)

both_inf <- comb[comb$mouse_inflection == "P" | comb$human_inflection == "P", ]
dim(both_inf)
# [1] 1081   13

corr <- cor.test(both_inf$mouse_RTS, both_inf$human_RTS, method = "spearman")
# 	Spearman's rank correlation rho

# data:  both_inf$mouse_RTS and both_inf$human_RTS
# S = 315539657, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.4789995  

ggscatter(both_inf, x = "mouse_RTS", y = "human_RTS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman")


# Check normality of variable1
# hist(both_inf$mouse_RTS)
# qqnorm(both_inf$mouse_RTS)
# shapiro.test(both_inf$mouse_RTS)


length(which(both_inf$mouse_inflection == 'NP' & both_inf$human_inflection == 'P'))
# [1] 338 		# Significant in human but not mouse
length(which(both_inf$mouse_inflection == 'P' & both_inf$human_inflection == 'NP'))
# [1] 224		# Significant in mouse but not human

####################################################

# Checks

hum_specific <- both_inf[both_inf$mouse_inflection == "NP" & both_inf$human_inflection == "P", ]
hum_specific <- hum_specific[order(hum_specific$human_RTS, decreasing = T),] 

mou_specific <- both_inf[both_inf$mouse_inflection == "P" & both_inf$human_inflection == "NP", ]
mou_specific <- mou_specific[order(mou_specific$mouse_RTS, decreasing = T),] 


both_inf$outliers <- (both_inf$human_RTS - both_inf$mouse_RTS)
both_inf <- both_inf[order(both_inf$outliers, decreasing = T),] 
# found weird genes starting with "^gm|^bc|^klk" contributing to noise- need to remove


##############################
# Species-specific Percentages

table(both_inf$mouse_inflection)
#  NP   P 
# 338 743 
table(both_inf$human_inflection)
#  NP   P 
# 224 857
dim(both_inf)
# [1] 1081   14

length(which(both_inf$mouse_inflection == 'NP' & both_inf$human_inflection == 'NP'))
# [1] 0
length(which(both_inf$mouse_inflection == 'P' & both_inf$human_inflection == 'P'))
# [1] 519

519/1081
# [1] 0.480111

hum_p <- both_inf[both_inf$human_inflection == "P", ]
dim(hum_p)
# [1] 857  14
table(hum_p$mouse_inflection)
#  NP   P 
# 338 519

519/857
# [1] 0.6056009






