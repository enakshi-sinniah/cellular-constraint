
# Analysis of Sarropoulos (Nature, 2019) LncRNA data
#Sarrapoulos et al. 2019 paper
#Developmental dynamics of lncRNAs across mammalian organs and species

# 1) Download full list of lncRNAs expressed during mammalian organ development

#https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1341-x/MediaObjects/41586_2019_1341_MOESM3_ESM.xlsx
# Saved above as: "Sarrapoulos_2019_supp.xlsx"
#Kept only Table S1: "Developmentally_dynamics_lncRNAs_full_dataset.txt"

library(data.table)

dd_all <- fread('Developmentally_dynamics_lncRNAs_full_dataset.txt')
nrow(dd_all) #31678 all lncrnas

dd <- dd_all[dd_all$Dynamic ==  TRUE] #5887 developmentally dynamic lncrnas
nn <- dd_all[dd_all$Dynamic ==  FALSE] #25791 other lncrnas

# 2) Generate RTS scores for full list of lncrnas

# Create input bedfile
dd_all$chr <- sub("^", "chr", dd_all$chr)
table(dd_all$chr)
dd_subs <- dd_all[, .(dd_all$chr, dd_all$Gene.left, dd_all$Gene.right, dd_all$Human_ID)]
write.table(dd_subs, 'DD_lncRNAs_bedtoolsinput.bed', quote = F, col.names = F, row.names = F, sep = '\t') 

#Get CCS
python /scratch/90days/uqesinni/TRIAGEi_analysis/TRIAGEi_main.py -i /scratch/90days/uqesinni/TRIAGEi_analysis/DD_lncRNAs_bedtoolsinput.bed -o /scratch/90days/uqesinni/TRIAGEi_analysis/DD_lncRNAs_all_v2rts.txt -r /scratch/90days/uqesinni/TRIAGEi_analysis/TRIAGEi_ReferenceH3K27me3_epimap_hg19.bedgraph

wc -l DD_lncRNAs_all_v2rts.txt #25232 i.e. lost 6446 lncrnas

#Format output file for inflection calc
library(inflection)
library(data.table)
xy <- fread('DD_lncRNAs_all_v2rts.txt')
xy <- xy[order(-V2),] #sort in ascending order based on v2rts
xy <- cbind(xy, "V3"=1:nrow(xy)) #add new column with rank
colnames(xy) <- c('Human_ID','RTS','RTS_Rank')

# 3) Calculate inflection point

plot(xy$V3, xy$V2, type= "l",xlab="Rank_order", ylab="RTS_score")
ede(xy$V3, xy$V2,index=0)
# j1 j2 chi
# EDE 1276  1 NaN
abline(v=1276,col="red") #1276 lncrnas have high RTS past the inflection point

# Merge rts file to get lncrna info
merged <- merge(xy, dd_all, by="Human_ID")
nrow(merged) #25232

write.table(merged, 'DD_lncRNAs_all_v2rts_withinfo.txt', quote = T, col.names = T, row.names = F, sep = '\t') 

# 4) Binned Fishers Exact Test (FET) measure enrichment of Developmentally Dynamic LncRNAs in RTS priority list

merged <- merged[order(RTS_Rank)] #sort by rank
all <- merged[, .(merged$Human_ID)] #keep only Human_ID column
nrow(all) #25232

priority <- all[1:1276] #priority lncrnas

ddyn <- dd_all[Dynamic %in% c('TRUE')] #developmentally dynamic lncrna positive list
ddyn <- ddyn[, .(ddyn$Human_ID)]
nrow(ddyn) #5887

# Setup FET parameters

'%!in%' <- function(x,y)!('%in%'(x,y))

m_table <- ddyn
N_table <- all
k_table <- priority

# binning
bin <- rep(c(1:10), each = length(k_table$V1)/10)
bin <- c(bin, rep(11, length(k_table$V1)-length(bin)))
k_table$bin <- bin

pvals <- c()
for (b in unique(k_table$bin)[1:11]){
  m <- m_table$V1
  N <- N_table$V1
  n <- N_table$V1[which(N_table$V1 %!in% m_table$V1)]
  k <- k_table[which(k_table$bin == b),]$V1
  
  x <- length(which(k %in% m))
  
  A <- x
  B <- length(k) - x
  C <- length(m) - x
  D <- length(n) - (length(k) - x)
  
  mat <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)
  result <- fisher.test(mat, alternative='greater')
  
  pvals <- c(pvals, result$p.value)
  print(pvals)
}

plt <- data.frame('bin' = c(1:11), 'pval' = pvals)
plt$padj <- p.adjust(plt$pval, method = 'fdr')

library(ggplot2)
ggplot(plt, aes(bin, -log10(padj))) + geom_point(size=1) + geom_line(size=0.8)
ggsave('FET_DDlncs_vs_TRIAGElncs_binned.pdf', width=12, height=10, dpi=600)

# 5) Overall Fishers Exact Test (FET) to test for enrichment of Developmentally Dynamic LncRNAs in RTS priority list

library(data.table)
'%!in%' <- function(x,y)!('%in%'(x,y))

dd_all <- fread('DD_lncRNAs_all_v2rts_withinfo.txt', header = T)
dd <- dd_all[dd_all$Dynamic ==  TRUE] #5887 developmentally dynamic lncrnas
nrow(dd_all) #25232
nrow(dd) #5215

dd_all <- dd_all[order(RTS_Rank),]
dd_all$Priority <- "NA" #add priority column
dd_all[dd_all$RTS_Rank <= 1276,]$Priority <- "P" 
dd_all[dd_all$RTS_Rank > 1276,]$Priority <- "NP"


#Overall Barplot
library(ggplot2)

ggplot(dd_all) +
  aes(x = 1, fill = Priority) +
  geom_bar(position = "fill")
ggsave('Barplot_all_TRIAGElncs_P_vs_NP.pdf', width=12, height=10, dpi=600)


#Generate contingency table

cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 

temp <- dd_all[dd_all$Priority == "P",] #extract all RTS priority lncrnas
dim(temp) 
# 1276   31

temp <- temp[which(temp$Human_ID %in% dd$Human_ID),] #subset Q1: all priority lncRNAs which are in DD positive gene list
dim(temp)
# 431  31
cmx[1,1] <- 431 #assign Q1 to contingency matrix

cmx[2,1] <- 1276-431 #subset Q2: all RTS priority lncRNAs which are NOT in lnctard positive gene list 
cmx
#     [,1] [,2]
# [1,]  431   NA
# [2,]  845   NA

temp <- dd_all[dd_all$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(temp)
# [1] 23956

temp <- temp[which(temp$Human_ID %in% dd$Human_ID),]  #subset Q3: all RTS non-priority lncRNAs which are in DD positive gene list
nrow(temp)
# [1] 4784
cmx[1,2] <- 4786 #assign Q3

cmx[2,2] <- 23956 - 4786 #subset Q4: all RTS non-priority lncRNAs which are NOT in lnctard positive gene list  (i.e. all others)
cmx
#      [,1]  [,2]
# [1,]  431  4786
# [2,]  845 19170

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("DD_yes", "DD_no")
cmx
#         P    NP
# DD_yes 431  4786
# DD_no  845 19170


#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                P        NP
# DD_yes  263.8274  4953.173
# DD_no  1012.1726 19002.827
#NOTE: expected frequencies show no values below 5 so run a chi-square

#Run chi-square test
test <- chisq.test(cmx)
test
#   Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 139.81, df = 1, p-value < 2.2e-16

test$p.value
#[1] 2.927392e-32


#Run FET

test2 <- fisher.test(cmx)
#   Fisher's Exact Test for Count Data
# data:  cmx
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.807190 2.307094
# sample estimates:
# odds ratio 
#   2.042988
test2$p.value
# [1] 4.631395e-29


# 5A) Barplot visualization

library(ggplot2)

ggplot(dd_all) +
  aes(x = Dynamic, fill = Priority) +
  geom_bar(position = "fill")
ggsave('Barplot_DevDyn_vs_TRIAGElncs.pdf', width=12, height=10, dpi=600)

#Flipped
ggplot(dd_all) +
  aes(x = Priority, fill = Dynamic) +
  geom_bar(position = "fill")
ggsave('Barplot_DevDyn_vs_TRIAGElncs_flipped.pdf', width=12, height=10, dpi=600)


# 6) Wilcox test and boxplot for DevDyb vs. CCS 

wilcox.test(dd_all[dd_all$Dynamic == "TRUE",]$RTS, dd_all[dd_all$Dynamic == "FALSE",]$RTS)

#   Wilcoxon rank sum test with continuity correction

# data:  dd_all[dd_all$Dynamic == "TRUE", ]$RTS and dd_all[dd_all$Dynamic == "FALSE", ]$RTS
# W = 58243078, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0


p <- ggplot(dd_all, aes(x = Dynamic, y = log10(abs(RTS)))) + geom_boxplot()
pdf("Boxplot_DevDynlncs_vs_RTS.pdf")
print(p) #save pdf
dev.off()

#_________________________________________________________________________________________________________________

# Reference: LncRNAdb functional lncrnas

library(data.table)
'%!in%' <- function(x,y)!('%in%'(x,y))
library(ggplot2)

dd_all <- fread('DD_lncRNAs_all_v2rts_withinfo.txt', header = T)
dd <- dd_all[dd_all$Dynamic ==  TRUE] #5887 developmentally dynamic lncrnas
nrow(dd_all) #25232
nrow(dd) #5215

dd_all <- dd_all[order(RTS_Rank),]
dd_all$Priority <- "NA" #add priority column
dd_all[dd_all$RTS_Rank <= 1276,]$Priority <- "P" 
dd_all[dd_all$RTS_Rank > 1276,]$Priority <- "NP"

lncdb <- dd_all[dd_all$lncRNAdb ==  TRUE] #73 lncRNAdb functional lncrnas
not <- dd_all[dd_all$lncRNAdb ==  FALSE] #25159 other lncrnas


#Generate contingency table
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 

temp <- dd_all[dd_all$Priority == "P",] #extract all RTS priority lncrnas
dim(temp) 
# 1276   31

temp <- temp[which(temp$Human_ID %in% lncdb$Human_ID),] #subset Q1
dim(temp)
# 17  31
cmx[1,1] <- 17 #assign Q1 to contingency matrix

cmx[2,1] <- 1276-17 #subset Q2
cmx
#      [,1] [,2]
# [1,]   17   NA
# [2,] 1259   NA

temp <- dd_all[dd_all$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(temp)
# [1] 23956

temp <- temp[which(temp$Human_ID %in% lncdb$Human_ID),]  #subset Q3
nrow(temp)
# [1] 56
cmx[1,2] <- 56 #assign Q3

cmx[2,2] <- 23956 - 56 #subset Q4
cmx
#      [,1]  [,2]
# [1,]   17    56
# [2,] 1259 23900

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("lncdb_yes", "lncdb_no")
cmx
#             P    NP
# lncdb_yes   17    56
# lncdb_no  1259 23900


#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                     P          NP
# lncdb_yes    3.691661    69.30834
# lncdb_no  1272.308339 23886.69166
# Warning message:
# In chisq.test(cmx) : Chi-squared approximation may be incorrect
#NOTE: expected frequencies have values below 5 so run a FET

#Run chi-square test
test <- fisher.test(cmx)
test
#   Fisher's Exact Test for Count Data

# data:  cmx
# p-value = 9.978e-08
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   3.129092 10.095868
# sample estimates:
# odds ratio 
#   5.761932 

test$p.value
#[1] 9.977771e-08


ggplot(dd_all) +
  aes(x = lncRNAdb, fill = Priority) +
  geom_bar(position = "fill")
ggsave('Barplot_lncRNAdb_vs_TRIAGEpriority.pdf', width=12, height=10, dpi=600)


#______________________________________________________________________________________________________________

# Reference: Evolutionary conservation of lncrnas

library(data.table)
'%!in%' <- function(x,y)!('%in%'(x,y))
library(ggplot2)

dd_all <- fread('DD_lncRNAs_all_v2rts_withinfo.txt', header = T)
dd <- dd_all[dd_all$Dynamic ==  TRUE] #5887 developmentally dynamic lncrnas
nrow(dd_all) #25232
nrow(dd) #5215

dd_all <- dd_all[order(RTS_Rank),]
dd_all$Priority <- "NA" #add priority column
dd_all[dd_all$RTS_Rank <= 1276,]$Priority <- "P" 
dd_all[dd_all$RTS_Rank > 1276,]$Priority <- "NP"

table(dd_all$Age_orthoMCL)
#    180Mya              25Mya             300Mya              90Mya 
#       598               4241                354               2281 
# ambiguous     human-specific multimember_family 
#       221              16782                755 

hspec <- dd_all[dd_all$Age_orthoMCL ==  "human-specific"]
mya_25 <- dd_all[dd_all$Age_orthoMCL ==  "25Mya"]
mya_90 <- dd_all[dd_all$Age_orthoMCL ==  "90Mya"]
mya_180 <- dd_all[dd_all$Age_orthoMCL ==  "180Mya"] 
mya_300 <- dd_all[dd_all$Age_orthoMCL ==  "300Mya"] 


ggplot(dd_all) +
  aes(x = Age_orthoMCL, fill = Priority) +
  geom_bar(position = "fill")
ggsave('Barplot_EvolutionCons_lncRNAdb_vs_TRIAGEpriority.pdf', width=12, height=10, dpi=600)


#_______________________________________________________________________________________________________________

# FET for all lncrnas vs. functional lncRNAdb

library(data.table)
'%!in%' <- function(x,y)!('%in%'(x,y))
library(ggplot2)

dd_all <- fread('DD_lncRNAs_all_v2rts_withinfo.txt', header = T)
dd <- dd_all[dd_all$Dynamic ==  TRUE] #5887 developmentally dynamic lncrnas
nrow(dd_all) #25232
nrow(dd) #5215

dd_all <- dd_all[order(RTS_Rank),]
dd_all$Priority <- "NA" #add priority column
dd_all[dd_all$RTS_Rank <= 1276,]$Priority <- "P" 
dd_all[dd_all$RTS_Rank > 1276,]$Priority <- "NP"

lncdb <- dd_all[dd_all$lncRNAdb ==  TRUE] #73 lncRNAdb functional lncrnas
not <- dd_all[dd_all$lncRNAdb ==  FALSE] #25159 other lncrnas

########################################################################################

# FET 5) 90Mya vs. 300Mya
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(mya_90, mya_300) #2879
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 200  31

temp <- priority[which(priority$Human_ID %in% mya_90$Human_ID),] #subset Q1
dim(temp)
# 150  31
cmx[1,1] <- 150 #assign Q1 to contingency matrix
cmx[2,1] <- 200-150 #subset Q2

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 2435
temp <- NP[which(NP$Human_ID %in% mya_90$Human_ID),]  #subset Q3
nrow(temp)
# [1] 2131
cmx[1,2] <- 2131 #assign Q3
cmx[2,2] <- 2435 - 2131 #subset Q4

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("90Mya", "300Mya")
cmx
#          P   NP
# 90Mya  150 2131
# 300Mya  50  304

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                P        NP
# 90Mya  173.13093 2107.8691
# 300Mya  26.86907  327.1309
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
#   Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 23.828, df = 1, p-value = 1.053e-06

test$p.value
#[1] 1.053343e-06


###################################################################################################################

# Reference: Evolutionary conservation of lncrnas

library(data.table)
'%!in%' <- function(x,y)!('%in%'(x,y))
library(ggplot2)

dd_all <- fread('DD_lncRNAs_all_v2rts_withinfo.txt', header = T)
dd <- dd_all[dd_all$Dynamic ==  TRUE] #5887 developmentally dynamic lncrnas
nrow(dd_all) #25232
nrow(dd) #5215

dd_all <- dd_all[order(RTS_Rank),]
dd_all$Priority <- "NA" #add priority column
dd_all[dd_all$RTS_Rank <= 1276,]$Priority <- "P" 
dd_all[dd_all$RTS_Rank > 1276,]$Priority <- "NP"

table(dd_all$Age_orthoMCL)
#    180Mya              25Mya             300Mya              90Mya 
#       598               4241                354               2281 
# ambiguous     human-specific multimember_family 
#       221              16782                755 

hspec <- dd_all[dd_all$Age_orthoMCL ==  "human-specific"]
mya_25 <- dd_all[dd_all$Age_orthoMCL ==  "25Mya"]
mya_90 <- dd_all[dd_all$Age_orthoMCL ==  "90Mya"]
mya_180 <- dd_all[dd_all$Age_orthoMCL ==  "180Mya"] 
mya_300 <- dd_all[dd_all$Age_orthoMCL ==  "300Mya"] 


ggplot(dd_all) +
  aes(x = Age_orthoMCL, fill = Priority) +
  geom_bar(position = "fill")
ggsave('Barplot_EvolutionCons_lncRNAdb_vs_TRIAGEpriority.pdf', width=12, height=10, dpi=600)


#FETs for all combinations run in ('Evolutionary_conserved_lncs_FETS.R')


# FETs for Evolutionary conservation

# Sarrapoulos data

# FET 1) Human-specific vs. 25Mya
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(hspec, mya_25) #21023
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 937  31

temp <- priority[which(priority$Human_ID %in% hspec$Human_ID),] #subset Q1
dim(temp)
# 699  31
cmx[1,1] <- 699 #assign Q1 to contingency matrix
cmx[2,1] <- 937-699 #subset Q2
cmx
#    [,1] [,2]
# [1,]  699   NA
# [2,]  238   NA

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 20086
temp <- NP[which(NP$Human_ID %in% hspec$Human_ID),]  #subset Q3
nrow(temp)
# [1] 16083
cmx[1,2] <- 16083 #assign Q3
cmx[2,2] <- 20086 - 16083 #subset Q4
cmx
#  [,1]  [,2]
# [1,]  699 16083
# [2,]  238  4003

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("human_specific", "25Mya")
cmx
#                 P    NP
# human_specific 699 16083
# 25Mya          238  4003


#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                 P        NP
# human_specific 747.9776 16034.022
# 25Mya          189.0224  4051.978
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
#   Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 16.301, df = 1, p-value = 5.403e-05

test$p.value
#[1] 5.402709e-05





cmx <- matrix(nrow = 2, ncol = 2)
> cmx[1,1] <- 64
> cmx[1,2] <- 340
> cmx[2,1] <- 1867
> cmx[2,2] <- 37286




#############################################################
# FET 2) 25Mya vs. 90Mya
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(mya_25, mya_90) #6522
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 388  31

temp <- priority[which(priority$Human_ID %in% mya_25$Human_ID),] #subset Q1
dim(temp)
# 238  31
cmx[1,1] <- 238 #assign Q1 to contingency matrix
cmx[2,1] <- 388-238 #subset Q2
cmx
#  [,1] [,2]
# [1,]  238   NA
# [2,]  150   NA

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 6134
temp <- NP[which(NP$Human_ID %in% mya_25$Human_ID),]  #subset Q3
nrow(temp)
# [1] 4003
cmx[1,2] <- 4003 #assign Q3
cmx[2,2] <- 6134 - 4003 #subset Q4
cmx
# [,1] [,2]
# [1,]  238 4003
# [2,]  150 2131

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("25Mya", "90Mya")
cmx
#         P   NP
# 25Mya 238 4003
# 90Mya 150 2131


#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#           P       NP
# 25Mya 252.3011 3988.699
# 90Mya 135.6989 2145.301
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
#   Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 2.2951, df = 1, p-value = 0.1298

test$p.value
#[1] 0.1297823

#############################################################
# FET 3) 25Mya vs. 180Mya
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(mya_25, mya_180)
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 298  31

temp <- priority[which(priority$Human_ID %in% mya_25$Human_ID),] #subset Q1
dim(temp)
# 238  31
cmx[1,1] <- 238 #assign Q1 to contingency matrix
cmx[2,1] <- 298-238 #subset Q2

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 4541
temp <- NP[which(NP$Human_ID %in% mya_25$Human_ID),]  #subset Q3
nrow(temp)
# [1] 4003
cmx[1,2] <- 4003 #assign Q3
cmx[2,2] <- 4541 - 4003 #subset Q4

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("25Mya", "180Mya")
cmx
#          P   NP
# 25Mya  238 4003
# 180Mya  60  538

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                P        NP
# 25Mya  261.17338 3979.8266
# 180Mya  36.82662  561.1734
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
#   Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 16.973, df = 1, p-value = 3.791e-05

test$p.value
#[1]   3.791232e-05

#############################################################
# FET 4) 90Mya vs. 180Mya
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(mya_90, mya_180) #2879
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 210  31

temp <- priority[which(priority$Human_ID %in% mya_90$Human_ID),] #subset Q1
dim(temp)
# 150  31
cmx[1,1] <- 150 #assign Q1 to contingency matrix
cmx[2,1] <- 210-150 #subset Q2

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 2669
temp <- NP[which(NP$Human_ID %in% mya_90$Human_ID),]  #subset Q3
nrow(temp)
# [1] 2131
cmx[1,2] <- 2131 #assign Q3
cmx[2,2] <- 2669 - 2131 #subset Q4

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("90Mya", "180Mya")
cmx
#         P   NP
# 25Mya 238 4003
# 90Mya 150 2131

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                P        NP
# 90Mya  166.38069 2114.6193
# 180Mya  43.61931  554.3807
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
#   Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 7.8717, df = 1, p-value = 0.005021

test$p.value
#[1]  0.005021405

#############################################################
# FET 5) 90Mya vs. 300Mya
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(mya_90, mya_300) #2879
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 200  31

temp <- priority[which(priority$Human_ID %in% mya_90$Human_ID),] #subset Q1
dim(temp)
# 150  31
cmx[1,1] <- 150 #assign Q1 to contingency matrix
cmx[2,1] <- 200-150 #subset Q2

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 2435
temp <- NP[which(NP$Human_ID %in% mya_90$Human_ID),]  #subset Q3
nrow(temp)
# [1] 2131
cmx[1,2] <- 2131 #assign Q3
cmx[2,2] <- 2435 - 2131 #subset Q4

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("90Mya", "300Mya")
cmx
#          P   NP
# 90Mya  150 2131
# 300Mya  50  304

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                P        NP
# 90Mya  173.13093 2107.8691
# 300Mya  26.86907  327.1309
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
#   Pearson's Chi-squared test with Yates' continuity correction

# data:  cmx
# X-squared = 23.828, df = 1, p-value = 1.053e-06

test$p.value
#[1] 1.053343e-06








#########################################################################################################################

# 28/09/21

# SyntDB data

# Step 4B) Run FETs
# 4B) Protein-coding SyntDB

# FET 1) Human-specific vs. Apes
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(hspec, apes)
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 30  27

temp <- priority[which(priority$Ensembl_reference_Transcript %in% hspec$Ensembl_reference_Transcript),] #subset Q1
dim(temp)
# 27 27
cmx[1,1] <- 27 #assign Q1 to contingency matrix
cmx[2,1] <- 30-27 #subset Q2

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 1518
temp <- NP[which(NP$Ensembl_reference_Transcript %in% hspec$Ensembl_reference_Transcript),]  #subset Q3
nrow(temp)
# [1] 1407
cmx[1,2] <- 1407 #assign Q3
cmx[2,2] <- 1518 - 1407 #subset Q4

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("human_specific", "Apes")
cmx
#                 P   NP
# human_specific 27 1407
# Apes            3  111

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                        P        NP
# human_specific 27.790698 1406.2093
# Apes            2.209302  111.7907
#NOTE: expected frequencies have values below 5 so run a fischers exact test

#Run chi-square test
test <- fisher.test(cmx)
test
#   Fisher's Exact Test for Count Data
# data:  cmx
# p-value = 0.4807
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.2134078 3.7145500
# sample estimates:
# odds ratio 
#  0.7102069 

###########################################################################################
# FET 2) Human-specific vs. Conserved
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(hspec, cons) #10225
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 94  27

temp <- priority[which(priority$Ensembl_reference_Transcript %in% hspec$Ensembl_reference_Transcript),] #subset Q1
dim(temp)
# 27 27
cmx[1,1] <- 27 #assign Q1 to contingency matrix
cmx[2,1] <- 94-27 #subset Q2

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 3201
temp <- NP[which(NP$Ensembl_reference_Transcript %in% hspec$Ensembl_reference_Transcript),]  #subset Q3
nrow(temp)
# [1] 1407
cmx[1,2] <- 1407 #assign Q3
cmx[2,2] <- 3201 - 1407 #subset Q4

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("human_specific", "Conserved")
cmx
#                 P   NP
# human_specific 27 1407
# Conserved      67 1794

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                       P       NP
# human_specific 40.90926 1393.091
# Conserved      53.09074 1807.909
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
#   Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 8.0106, df = 1, p-value = 0.00465



# Run FET

test2 <- fisher.test(cmx)
#   Fisher's Exact Test for Count Data
# data:  cmx
# p-value = 0.003075
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.3142051 0.8195808
# sample estimates:
# odds ratio 
#  0.5139309 

test2$p.value
# [1] 0.003074816


################################################################
################################################################

# Step 4B) Run FETs
# 4B) Non-coding SyntDB

# FET 1) Human-specific vs. Apes
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(hspec, apes) #3170
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 118  27

temp <- priority[which(priority$Ensembl_reference_Transcript %in% hspec$Ensembl_reference_Transcript),] #subset Q1
dim(temp)
# 95 27
cmx[1,1] <- 95 #assign Q1 to contingency matrix
cmx[2,1] <- 118-95 #subset Q2

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 3052
temp <- NP[which(NP$Ensembl_reference_Transcript %in% hspec$Ensembl_reference_Transcript),]  #subset Q3
nrow(temp)
# [1] 2694
cmx[1,2] <- 2694 #assign Q3
cmx[2,2] <- 3052 - 2694 #subset Q4

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("human_specific", "Apes")
cmx
#                 P   NP
# human_specific 95 2694
# Apes           23  358

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                        P        NP
# human_specific 103.81767 2685.1823
# Apes            14.18233  366.8177
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
#   Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 5.7589, df = 1, p-value = 0.01641

test$p.value
#[1] 0.0164052

###########################################################################################
# FET 2) Human-specific vs. Conserved
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(hspec, cons) #10225
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 515  27

temp <- priority[which(priority$Ensembl_reference_Transcript %in% hspec$Ensembl_reference_Transcript),] #subset Q1
dim(temp)
# 95 27
cmx[1,1] <- 95 #assign Q1 to contingency matrix
cmx[2,1] <- 515-95 #subset Q2

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 9710
temp <- NP[which(NP$Ensembl_reference_Transcript %in% hspec$Ensembl_reference_Transcript),]  #subset Q3
nrow(temp)
# [1] 2694
cmx[1,2] <- 2694 #assign Q3
cmx[2,2] <- 9710 - 2694 #subset Q4

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("human_specific", "Conserved")
cmx
#                  P   NP
# human_specific  95 2694
# Conserved      420 7016

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                       P       NP
# human_specific 140.4729 2648.527
# Conserved      374.5271 7061.473
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
#   Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 20.849, df = 1, p-value = 4.971e-06

test$p.value
#[1] 4.970589e-06

###########################################################################################
# FET 3) Human-specific vs. Ultraconserved
cmx <- matrix(nrow = 2, ncol = 2) #create contingency table 
comp1 <- rbind(hspec, ultra)
priority <- comp1[comp1$Priority == "P",] #extract all RTS priority lncrnas
dim(priority) 
# [1] 107  27

temp <- priority[which(priority$Ensembl_reference_Transcript %in% hspec$Ensembl_reference_Transcript),] #subset Q1
dim(temp)
# 95 27
cmx[1,1] <- 95 #assign Q1 to contingency matrix
cmx[2,1] <- 107-95 #subset Q2

NP <- comp1[comp1$Priority == "NP",] #extract all RTS non-priority lncrnas
nrow(NP)
# [1] 3504
temp <- NP[which(NP$Ensembl_reference_Transcript %in% hspec$Ensembl_reference_Transcript),]  #subset Q3
nrow(temp)
# [1] 2694
cmx[1,2] <- 2694 #assign Q3
cmx[2,2] <- 3504 - 2694 #subset Q4

colnames(cmx) <- c("P", "NP") 
rownames(cmx) <- c("human_specific", "Ultraconserved")
cmx
#                 P   NP
# human_specific 95 2694
# Ultraconserved 12  810

#Check expected frequencies to know whether to run Chi-square test or fischer's exact test
chisq.test(cmx)$expected
#                       P        NP
# human_specific 82.64276 2706.3572
# Ultraconserved 24.35724  797.6428
#NOTE: expected frequencies have values over 5 so run a chisquare

#Run chi-square test
test <- chisq.test(cmx)
test
#     Pearson's Chi-squared test with Yates' continuity correction
# data:  cmx
# X-squared = 7.7016, df = 1, p-value = 0.005517

test$p.value
#[1] 0.005517155

