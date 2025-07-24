
###################################################################################################################
################################## Evolutionary Constraint Head-to-head- FETs Comparisons ###############################################
###################################################################################################################

# AIM: Compare significantly constrained gene lists- run FETs, dig out specific examples.

###############################################################################
############# Compare gene enrichments ########################################
###############################################################################

# Step 1) Get genes overlapping and non-overlapping CCS/Phast

library(data.table)

ccs_gene <- fread('240701_CCS_pergene_PCG_hg38.txt')
phast_gene <- fread('240701_Phastcons_consscores_pergene_PCG_hg38.txt')

ccs_gene$ccs_rank <- 1:18075 
ccs_gene$ccs_inflection <- ifelse(seq_len(nrow(ccs_gene)) <= 857, 'P', 'NP')

phast_gene$phast_rank <- 1:20052 
phast_gene$phast_inflection <- ifelse(seq_len(nrow(phast_gene)) <= 943, 'P', 'NP')

comb <- merge(phast_gene, ccs_gene, by = "Gene_id", all.x= TRUE)
dim(comb)
# [1] 20052     7

# Fill empty rows
comb[, ccs_inflection := fifelse(is.na(ccs_inflection) | ccs_inflection == "", "NP", ccs_inflection)]
comb[, phast_inflection := fifelse(is.na(phast_inflection) | phast_inflection == "", "NP", phast_inflection)]

########################################################################

# Subset overlapping sets

# Both CCS+Phastcons significant
both_sig <- comb[phast_inflection == "P" & ccs_inflection == "P"]
dim(both_sig)
# [1] 55  7
write.table(both_sig, '240701_Both_CCSPhast_sig_pergene.txt', quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

# CCS only sig
ccs_sig <- comb[phast_inflection == "NP" & ccs_inflection == "P"]
dim(ccs_sig)
# [1] 802   7
write.table(ccs_sig, '240701_CCS_only_sig_pergene.txt', quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

# Phastcons only sig
phast_sig <- comb[phast_inflection == "P" & ccs_inflection == "NP"]
dim(phast_sig)
# [1] 888   7
write.table(phast_sig, '240701_Phast_only_sig_pergene.txt', quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

###########################################################################

# Step 2) Run FET for CCS/Phastcons

###### Create contingency matrix
# FET 1: CCS vs. Phastcons
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
all_genes <- comb
ccs_pos <- all_genes[all_genes$ccs_inflection == "P",] 
dim(ccs_pos) 
# [1] 857   7
phast_pos <- all_genes[all_genes$phast_inflection == "P",] 
dim(phast_pos)
# [1] 943   7
temp <- ccs_pos[which(ccs_pos$Gene_id %in% phast_pos$Gene_id),] #subset Q1
dim(temp)
# [1] 55  7
cmx[1,1] <- 55 #assign Q1
cmx[1,2] <- 943-55 #assign Q2

ccs_neg <- all_genes[all_genes$ccs_inflection == "NP",] 
dim(ccs_neg)
# [1] 19195     7
phast_neg <- all_genes[all_genes$phast_inflection == "NP",] 
dim(phast_neg)
# [1] 19109   7

temp <- phast_neg[which(phast_neg$Gene_id %in% ccs_pos$Gene_id),] #assign Q3
nrow(temp)
# [1] 802
cmx[2,1] <- 802 #assign Q3
cmx[2,2] <- 19109 - 802 #subset Q4

colnames(cmx) <- c("CCS_sig", "CCS_ns") 
rownames(cmx) <- c("Phastcons_sig", "Phastcons_ns")
cmx
#               CCS_sig CCS_ns
# Phastcons_sig      55    888
# Phastcons_ns      802  18307

###### Run Fishers exact test

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                 CCS_sig     CCS_ns
# Phastcons_sig  40.30276   902.6972
# Phastcons_ns  816.69724 18292.3028
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 5.4823, df = 1, p-value = 0.01921

test$p.value
# [1] 0.01921005

###### Barplot of proportions

library(ggplot2)
library(reshape2)

ggplot(all_genes) +
  aes(x = ccs_inflection, fill = phast_inflection) +
  geom_bar(position = "fill")
ggsave('240701_Barplot_CCSsig_vs_Phastconssig.pdf', width=12, height=10, dpi=600)


###########################################################################
# Step 2) Run FET comparing CCS/Phastcons enrichment of VETFs
############################################################################

library(data.table)
vetfs <- fread('VETFs_list.txt', header = FALSE)  
vetfs <- as.data.frame(vetfs)

# Prepare testing lists
comb[, short_id := gsub("_.*", "", Gene_id)]
dim(comb)
# [1] 20052     8
comb2 <- comb[!is.na(short_id) & short_id != ""]
dim(comb2)
# [1] 19449     8

write.table(comb2, 'SigGenes_CCSandPhastcons_FETinput.txt', quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

length(which(vetfs$V1 %in% comb2$short_id))
# [1] 629
dim(vetfs)
# [1] 634   1

vetfs$VETF <- "VETF"

comb3 <- merge(comb2, vetfs, by.x = "short_id",by.y = "V1", all.x= TRUE)
# Fill empty rows
comb3[, VETF := fifelse(is.na(VETF), "Non-VETF", VETF)]

write.table(comb3, 'AnnotGenes_CCSandPhastcons_FETinput.txt', quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

#####################################
# 2.1) FET for CCS vs. VETFs

'%!in%' <- function(x,y)!('%in%'(x,y))

all_genes <- comb3
ccs_pos <- all_genes[all_genes$ccs_inflection == "P",] 
ccs_neg <- all_genes[all_genes$ccs_inflection == "NP",] 
vetf_pos <- all_genes[all_genes$VETF == "VETF",] 
vetf_neg <- all_genes[all_genes$VETF == "Non-VETF",]

Q1 <- length(which(ccs_pos$short_id %in% vetf_pos$short_id))
Q2 <- length(which(ccs_pos$short_id %!in% vetf_pos$short_id))
Q3 <- length(which(ccs_neg$short_id %in% vetf_pos$short_id))
Q4 <- length(which(ccs_neg$short_id %!in% vetf_pos$short_id))

cmx <- matrix(nrow = 2, ncol = 2)
cmx[1,1] <- Q1
cmx[2,1] <- Q2
cmx[1,2] <- Q3
cmx[2,2] <- Q4

colnames(cmx) <- c("CCS_sig", "CCS_ns") 
rownames(cmx) <- c("VETF", "Non-VETF")

###### Run Fishers exact test

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#            CCS_sig     CCS_ns
# VETF      27.17723   602.8228
# Non-VETF 811.82277 18007.1772
#NOTE: expected frequencies have values over 5 so run a chisquare

cmx
#          CCS_sig CCS_ns
# VETF         231    379
# Non-VETF     608  16361

#Run chi-square test
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 1642.9, df = 1, p-value < 2.2e-16

test$p.value
# [1] 0

###### Barplot of proportions

ggplot(all_genes) +
  aes(x = ccs_inflection, fill = VETF) +
  geom_bar(position = "fill")
ggsave('240701_Barplot_CCSsig_vs_VETF.pdf', width=12, height=10, dpi=600)

#####################################
# 2.2) FET for Phastcons vs. VETFs

all_genes <- comb3
phast_pos <- all_genes[all_genes$phast_inflection == "P",] 
phast_neg <- all_genes[all_genes$phast_inflection == "NP",] 
vetf_pos <- all_genes[all_genes$VETF == "VETF",] 
vetf_neg <- all_genes[all_genes$VETF == "Non-VETF",]

Q1 <- length(which(phast_pos$short_id %in% vetf_pos$short_id))
Q2 <- length(which(phast_pos$short_id %!in% vetf_pos$short_id))
Q3 <- length(which(phast_neg$short_id %in% vetf_pos$short_id))
Q4 <- length(which(phast_neg$short_id %!in% vetf_pos$short_id))

cmx <- matrix(nrow = 2, ncol = 2)
cmx[1,1] <- Q1
cmx[2,1] <- Q2
cmx[1,2] <- Q3
cmx[2,2] <- Q4

colnames(cmx) <- c("Phastcons_sig", "Phastcons_ns") 
rownames(cmx) <- c("VETF", "Non-VETF")
cmx
#          Phastcons_sig Phastcons_ns
# VETF                68          562
# Non-VETF           771        18048

###### Run Fishers exact test

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#          Phastcons_sig Phastcons_ns
# VETF          27.17723     602.8228
# Non-VETF     811.82277   18007.1772
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 1642.9, df = 1, p-value < 9.097e-16

test$p.value
# [1] 9.096527e-16

###### Barplot of proportions

ggplot(all_genes) +
  aes(x = phast_inflection, fill = VETF) +
  geom_bar(position = "fill")
ggsave('240701_Barplot_Phastsig_vs_VETF.pdf', width=12, height=10, dpi=600)

###### Joint Barplot of proportions

# Reshape the data into long format
all_genes_long <- melt(all_genes, id.vars = c("short_id", "VETF"), measure.vars = c("phast_inflection", "ccs_inflection"),
                       variable.name = "inflection_type", value.name = "inflection_value")

# Joint barplot
ggplot(all_genes_long, aes(x = inflection_value, fill = VETF)) +
  geom_bar(position = "fill", stat = "count") +
  facet_wrap(~ inflection_type, scales = "free_x") +
  labs(x = "Inflection Type", y = "Proportion", fill = "VETF") +
  theme_minimal()
ggsave('240701_Barplot_Joint_CCSPhast_vs_VETF.pdf', width=12, height=10, dpi=600)






###########################################################################
# Step 3) Run FET comparing CCS/Phastcons enrichment of GOs
############################################################################

# NOTE: Ran below for multiple lists of GO terms (see results below)

library(data.table)

comb3 <- fread('AnnotGenes_CCSandPhastcons_FETinput.txt')
set1 <- fread('GO_0003700_list.txt')

set1$GO_0003700 <- "Yes"
comb3 <- merge(comb3, set1, by.x = "short_id",by.y = "V1", all.x= TRUE)
# Fill empty rows
comb3[, GO_0003700 := fifelse(is.na(GO_0003700), "No", GO_0003700)]
dim(comb3)
# [1] 19449    13

write.table(comb3, 'AnnotGenes_CCSandPhastcons_FETinput.txt', quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

# 2.1) FET for CCS vs. GO:0003700

'%!in%' <- function(x,y)!('%in%'(x,y))
all_genes <- comb3
ccs_pos <- all_genes[all_genes$ccs_inflection == "P",] 
ccs_neg <- all_genes[all_genes$ccs_inflection == "NP",] 
GO_pos <- all_genes[all_genes$GO_0003700 == "Yes",] 
GO_neg <- all_genes[all_genes$GO_0003700 == "No",]

Q1 <- length(which(ccs_pos$short_id %in% GO_pos$short_id))
Q2 <- length(which(ccs_pos$short_id %!in% GO_pos$short_id))
Q3 <- length(which(ccs_neg$short_id %in% GO_pos$short_id))
Q4 <- length(which(ccs_neg$short_id %!in% GO_pos$short_id))

cmx <- matrix(nrow = 2, ncol = 2)
cmx[1,1] <- Q1
cmx[2,1] <- Q2
cmx[1,2] <- Q3
cmx[2,2] <- Q4

colnames(cmx) <- c("CCS_sig", "CCS_ns") 
rownames(cmx) <- c("GO_yes", "GO_no")
cmx
#        CCS_sig CCS_ns
# GO_yes     420   4000
# GO_no      419  14610

###### Run Fishers exact test

#Run chi-square test
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 1789.8, df = 1, p-value < 2.2e-16

test$p.value
# [1]  9.246393e-83

###### Barplot of proportions
library(ggplot2)

ggplot(all_genes) +
  aes(x = ccs_inflection, fill = GO_0003700) +
  geom_bar(position = "fill")

###############################################################
# 2.2) FET for Phastcons vs. GO:0003700

'%!in%' <- function(x,y)!('%in%'(x,y))
all_genes <- comb3
phast_pos <- all_genes[all_genes$phast_inflection == "P",] 
phast_neg <- all_genes[all_genes$phast_inflection == "NP",] 
GO_pos <- all_genes[all_genes$GO_0003700 == "Yes",] 
GO_neg <- all_genes[all_genes$GO_0003700 == "No",]

Q1 <- length(which(phast_pos$short_id %in% GO_pos$short_id))
Q2 <- length(which(phast_pos$short_id %!in% GO_pos$short_id))
Q3 <- length(which(phast_neg$short_id %in% GO_pos$short_id))
Q4 <- length(which(phast_neg$short_id %!in% GO_pos$short_id))

# length(which(ccs_pos$short_id %in% phast_pos$short_id))

cmx <- matrix(nrow = 2, ncol = 2)
cmx[1,1] <- Q1
cmx[2,1] <- Q2
cmx[1,2] <- Q3
cmx[2,2] <- Q4

colnames(cmx) <- c("phast_sig", "phast_ns") 
rownames(cmx) <- c("GO_yes", "GO_no")
cmx
#        phast_sig phast_ns
# GO_yes       243     4177
# GO_no        596    14433

###### Run Fishers exact test

#Run chi-square test
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 1789.8, df = 1, p-value 1.271e-05

test$p.value
# [1] 1.271482e-05

###### Barplot of proportions
library(ggplot2)

ggplot(all_genes) +
  aes(x = phast_inflection, fill = GO_0003700) +
  geom_bar(position = "fill")

###### Joint Barplot of proportions

# Reshape the data into long format
all_genes_long <- melt(all_genes, id.vars = c("short_id", "GO_0003700"), measure.vars = c("phast_inflection", "ccs_inflection"),
                       variable.name = "inflection_type", value.name = "inflection_value")

# Joint barplot
ggplot(all_genes_long, aes(x = inflection_value, fill = GO_0003700)) +
  geom_bar(position = "fill", stat = "count") +
  facet_wrap(~ inflection_type, scales = "free_x") +
  labs(x = "Inflection Type", y = "Proportion", fill = "GO_0003700") +
  theme_minimal()
ggsave('240701_Barplot_Joint_CCSPhast_vs_GO_0003700.pdf', width=12, height=10, dpi=600)




###############################################################################
#################### RESULTS ##################################################
###############################################################################
####### Repeated above for multiple lists of GO terms

#### 1) GO:0003700	DNA-binding transcription factor activity

# set1 <- fread('GO_0003700_list.txt')
# 240701_Barplot_Joint_CCSPhast_vs_GO_0003700.pdf

# 1.1) CCS results
cmx
#        CCS_sig CCS_ns
# GO_yes     295    627
# GO_no      544  17983

test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 1789.8, df = 1, p-value < 2.2e-16
test$p.value
# [1] 0

# 1.2) Phastcons results
cmx
#        phast_sig phast_ns
# GO_yes        98      824
# GO_no        741    17786

test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 91.917, df = 1, p-value < 2.2e-16
test$p.value
# [1] 9.039295e-22

########################################
#### 2) GO_0030154	cell differentiation

# set1 <- fread('GO_0030154_list.txt')
# 240701_Barplot_Joint_CCSPhast_vs_GO_0030154.pdf

# 2.1) CCS results
cmx
#        CCS_sig CCS_ns
# GO_yes     420   4000
# GO_no      419  14610
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 371.41, df = 1, p-value < 2.2e-16
test$p.value
# [1] 9.246393e-83

# 2.2) Phastcons results
cmx
#        phast_sig phast_ns
# GO_yes       243     4177
# GO_no        596    14433
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 19.053, df = 1, p-value = 1.271e-05
test$p.value
# [1] 1.271482e-05

###############################################
#### 3) GO:0001708    cell fate specification

# set1 <- fread('GO_0001708_list.txt')
# 240701_Barplot_Joint_CCSPhast_vs_GO_0001708.pdf

# 3.1) CCS results
cmx
#        CCS_sig CCS_ns
# GO_yes      54     58
# GO_no      785  18552
test <- chisq.test(cmx)
test
# Warning message:
# In chisq.test(cmx) : Chi-squared approximation may be incorrect
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 515.31, df = 1, p-value < 2.2e-16
test$p.value
# [1] 4.425855e-114

# 3.2) Phastcons results
cmx
#        phast_sig phast_ns
# GO_yes        10      102
# GO_no        829    18508
test <- chisq.test(cmx)
test
# Warning message:
# In chisq.test(cmx) : Chi-squared approximation may be incorrect
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 4.7417, df = 1, p-value = 0.02944
test$p.value
# [1] 0.0294408

###############################################
#### 4) GO:0045165    cell fate commitment

# set1 <- fread('GO_0045165_list.txt')
# 240701_Barplot_Joint_CCSPhast_vs_GO_0045165.pdf


# 4.1) CCS results
cmx
#        CCS_sig CCS_ns
# GO_yes     118    180
# GO_no      721  18430
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 904.09, df = 1, p-value < 2.2e-16
test$p.value
# [1] 1.268254e-198

# 4.2) Phastcons results
cmx
#        phast_sig phast_ns
# GO_yes        24      274
# GO_no        815    18336
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 9.355, df = 1, p-value = 0.002224
test$p.value
# [1] 0.00222373

###############################################
#### 5) GO:0009887	animal organ morphogenesis

# set1 <- fread('GO_0009887_list.txt')
# 240701_Barplot_Joint_CCSPhast_vs_GO_0009887.pdf

# 5.1) CCS results
cmx
#        CCS_sig CCS_ns
# GO_yes     217    926
# GO_no      622  17684
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 629.48, df = 1, p-value < 2.2e-16
test$p.value
# [1] 6.500775e-139

# 5.2) Phastcons results
cmx
#        phast_sig phast_ns
# GO_yes        67     1076
# GO_no        772    17534
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 6.6563, df = 1, p-value = 0.009881
test$p.value
# [1] 0.009880514

###############################################
#### 6) GO:0048729	tissue morphogenesis 

# set1 <- fread('GO_0048729_list.txt')
# 240701_Barplot_Joint_CCSPhast_vs_GO_0048729.pdf

# 6.1) CCS results
cmx
#        CCS_sig CCS_ns
# GO_yes     115    586
# GO_no      724  18024
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 254.54, df = 1, p-value < 2.2e-16
test$p.value
# [1] 2.661518e-57

# 6.2) Phastcons results
cmx
#        phast_sig phast_ns
# GO_yes        40      661
# GO_no        799    17949
test <- chisq.test(cmx)
test
# 	Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 3.0742, df = 1, p-value = 0.07955
test$p.value
# [1] 0.07954566





