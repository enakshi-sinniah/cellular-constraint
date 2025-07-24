# Human-Mouse CCS Conservation

# HPC = Bunya
# wd = "/scratch/user/uqesinni/TRIAGEi_analysis"
# RDM: /QRISdata/Q3856/01_TRIAGE_intergenic/04-Delta2_90days_tempdata/

#################################################################################################################
########################################## Human-Mouse CCS Conservation #########################################
#################################################################################################################

################################################
# Get CCS for mm10 100bp syntenic blocks

python TRIAGEi_main.py -i syntenic_consensus_mm10_100bpLiftover.bed -o syntenic_consensus_mm10_100bp_ccs.txt -r scaled_fold_change_RTS_mm10.bedgraph

########
#Get merged file with hg19-mm10 syntenic 100bp blocks and RTS values

library(data.table)

### Human hg19
hum_rts <- fread('syntenic_consensus_hg19_100bp_rts.txt')
hh <- fread('syntenic_consensus_hg19_100bp.bed')
#get non-rts region_ids back and assign 0
hum_merged <- merge(hum_rts, hh, by.x = "V1" , by.y = "V4", all.y = TRUE)
colnames(hum_merged) <- c("region_id", "human_RTS", "hg19_chr", "hg19_start", "hg19_end" )
hum_merged$human_RTS <- ifelse(is.na(hum_merged$human_RTS), 0, hum_merged$human_RTS) #change NAs i.e. non-RTS regions to 0

### Mouse mm10
mou_rts <- fread('syntenic_consensus_mm10_100bp_ccs.txt')
mm <- fread('syntenic_consensus_mm10_100bpLiftover.bed')
#get non-rts region_ids back and assign 0
mou_merged <- merge(mou_rts, mm, by.x = "V1" , by.y = "V4", all.y = TRUE)
colnames(mou_merged) <- c("region_id", "mouse_RTS", "mm10_chr", "mm10_start", "mm10_end" )
mou_merged$mouse_RTS <- ifelse(is.na(mou_merged$mouse_RTS), 0, mou_merged$mouse_RTS) #change NAs i.e. non-RTS regions to 0

# Merge human and mouse metadata
comb <- merge(mou_merged, hum_merged, by = "region_id", all.x = TRUE)
dim(comb)
# [1] 2692660       9

length(which(comb$mouse_RTS == 0 & comb$human_RTS == 0))
# [1] 28451

comb$bin_mouse <- ifelse(comb$mouse_RTS > 0, 'Y', 'N')
comb$bin_human <- ifelse(comb$human_RTS > 0, 'Y', 'N')

# comb[order(-comb$human_RTS), ]

write.table(comb, '240117_combined_RTS_syntenic_100bp_meta.txt', quote = F, col.names = T, row.names = F, sep = '\t') 

##########################################################################################################

# Plot correlation 

library(ggplot2)
# install.packages("ggpubr")
library(ggpubr)
library(data.table)
library(ggrastr)

comb <- fread('240117_combined_RTS_syntenic_100bp_meta.txt')

comb_corr <- cor.test(comb$mouse_RTS, comb$human_RTS, method = "spearman")
comb_corr
# 	Spearman's rank correlation rho

# data:  comb$mouse_RTS and comb$human_RTS
# S = 1.3227e+18, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5934864 

# ggscatter(comb, x = "mouse_RTS", y = "human_RTS", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "spearman")


ggplot(comb, aes(x = mouse_RTS, y = human_RTS)) +
  geom_point(size = 1) +  
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +  
  labs(x = "mouse RTS", y = "human RTS") + 
  ggtitle(paste("Spearman Correlation =", round(comb_corr$estimate, 2), " (N =", length(comb$human_RTS), ")")) 
  ggsave("240117_Syntenic_100bp_ALL_rts_correlation.pdf", width = 8, height = 6)


# save rasterized pdf
p <- ggplot(comb, aes(x = mouse_RTS, y = human_RTS)) +
  geom_point(size = 0.7, color= "#6D6E71") +  
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +  
  labs(x = "mouse RTS", y = "human RTS") + 
  ggtitle(paste("Spearman Correlation =", round(comb_corr$estimate, 2), " (N =", length(comb$human_RTS), ")")) 

rastered_plot <- rasterize(p, width = 8, height = 6)
ggsave("240117_Syntenic_100bp_ALL_rts_correlation.pdf", plot = rastered_plot)

# scp uqesinni@bunya.rcc.uq.edu.au:/scratch/user/uqesinni/TRIAGEi_analysis/240117_Syntenic_100bp_ALL_rts_correlation.pdf .

# scp uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/07-Bunya_scratch_backup/240117_combined_RTS_syntenic_100bp_meta.txt .
 

###############

# Plot the above with density coloured

# if (!requireNamespace("devtools", quietly = TRUE))
#     install.packages("devtools")
# devtools::install_github("LKremer/ggpointdensity")

library(data.table)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)
library(ggpubr)
library(ggrastr)

dim(comb)
# [1] 2692660      11

ggplot(comb, aes(x = mouse_RTS, y = human_RTS)) +
  geom_pointdensity(size=0.7, adjust= 1000) + #adjust value changes n_neighbours for density
  scale_color_viridis() +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +
  labs(x = "mouse RTS", y = "human RTS") +
  ggtitle(paste("Spearman Correlation =", round(comb_corr$estimate, 2), " (N =", length(comb$human_RTS), ")")) +
  theme_minimal()


# save rasterized pdf
p <- ggplot(comb, aes(x = mouse_RTS, y = human_RTS)) +
  geom_pointdensity(size=0.7, adjust= 1000) + #adjust value changes n_neighbours for density
  scale_color_viridis() +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +
  labs(x = "mouse RTS", y = "human RTS") +
  ggtitle(paste("Spearman Correlation =", round(comb_corr$estimate, 2), " (N =", length(comb$human_RTS), ")")) +
  theme_minimal()

rastered_plot <- rasterize(p, width = 8, height = 6)
ggsave("240117_Syntenic_100bp_ALL_density_rts_correlation.pdf", plot = rastered_plot)

# scp uqesinni@bunya.rcc.uq.edu.au:/scratch/user/uqesinni/TRIAGEi_analysis/240117_Syntenic_100bp_ALL_density_rts_correlation.pdf .



###############################################################################
# Calculate inflection point 

library(inflection)

#### Human
human <- comb[comb$bin_human == "Y", ]
human <- human[order(human$human_RTS, decreasing = T),] 
dim(human) #1503358

human$human_rank <- 1:1503358 #formatting

plot(human$human_rank, human$human_RTS, type= "l", xlab="Rank_order", ylab="RTS_score") #plot all RTS SNPs

#calculate inflection point
ede(human$human_rank, human$human_RTS, index=0) 
#        j1 j2 chi
# EDE 79450 58 NaN

# Add priority info into comb metadata
comb <- comb[order(comb$human_RTS, decreasing = T),] 
comb$human_inflection <- ifelse(seq_len(nrow(comb)) <= 79450, 'P', 'NP')

#### Mouse
mouse <- comb[comb$bin_mouse == "Y", ]
mouse <- mouse[order(mouse$mouse_RTS, decreasing = T),] 
dim(mouse) #2653559

mouse$mouse_rank <- 1:2653559 #formatting

plot(mouse$mouse_rank, mouse$mouse_RTS, type= "l", xlab="Rank_order", ylab="RTS_score") #plot all RTS SNPs

#calculate inflection point
ede(mouse$mouse_rank, mouse$mouse_RTS, index=0) 
#         j1 j2 chi
# EDE 82188  10 NaN

# Add priority info into comb metadata
comb <- comb[order(comb$mouse_RTS, decreasing = T),] 
comb$mouse_inflection <- ifelse(seq_len(nrow(comb)) <= 82188, 'P', 'NP')

#######################

# Keep only regions that are P(sig) in both species

inflec <- comb[comb$mouse_inflection == "P" & comb$human_inflection == "P", ]
dim(inflec)
# [1] 47238      13
corr <- cor.test(inflec$mouse_RTS, inflec$human_RTS, method = "spearman")
# 	Spearman's rank correlation rho

# data:  inflec$mouse_RTS and inflec$human_RTS
# S = 5.3661e+12, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6945544 

summary(inflec$human_RTS)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05423 0.08783 0.13755 0.18507 0.21506 1.00000 
summary(inflec$mouse_RTS)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03679 0.05891 0.09540 0.14181 0.16322 1.00000 

ggscatter(inflec, x = "mouse_RTS", y = "human_RTS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman")


ggplot(inflec, aes(x = mouse_RTS, y = human_RTS)) +
  geom_point(size = 1) +  
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +  
  labs(x = "mouse RTS", y = "human RTS") + 
  ggtitle(paste("Spearman Correlation =", round(corr$estimate, 2), " (N =", length(inflec$human_RTS), ")")) 
  ggsave("250529_Syntenic_100bp_Priority_rts_correlation.pdf", width = 8, height = 6)

# scp uqesinni@bunya.rcc.uq.edu.au:/scratch/project/palpant_scratch/ESinniah/CellConstraint_project/250529_Syntenic_100bp_Priority_rts_correlation.pdf .


##############################################################

# Species specific comparisons

both_inf <- comb[comb$mouse_inflection == "P" | comb$human_inflection == "P", ]
dim(both_inf)
# [1] 114400    13

# How many significant regions have RTS in one species but not the other
table(both_inf$bin_mouse)
   #   N      Y 
   # 107 114293
table(both_inf$bin_human)
  #    N      Y 
  # 1755 112645

# How many significant regions are P(sig) in one species but not the other
table(both_inf$mouse_inflection)
#    NP     P 
# 32212 8218
table(both_inf$human_inflection)
#    NP     P 
# 34950 79450 

corr <- cor.test(both_inf$mouse_RTS, both_inf$human_RTS, method = "spearman")
# 	Spearman's rank correlation rho

# data:  both_inf$mouse_RTS and both_inf$human_RTS
# S = 1.732e+14, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3059183



mouse_specific <- inflec[inflec$bin_mouse == 'Y'  & inflec$bin_human == 'N', ]






