###################################################################################################################
################################## Histone Modification Comparisons ###############################################
###################################################################################################################

##############################################
######### ALL FUNCTIONAL ELEMENTS ############
##############################################

# Step 1) Get element ranks based on HM breadth 
# see '240620_Histone_mods_comparisons_genes.R' for pseudo-CCS breadth score generation
# 240620_Histone_mods_comparisons_lncrnas.R
# 240620_Histone_mods_comparisons_enhancers.R
# 240620_Histone_mods_comparisons_ZooCons.R

library(data.table)
library(dplyr)
library(ggplot2)

'%!in%' <- function(x,y)!('%in%'(x,y))

# 1) Run for TFs
files <- list.files(pattern = "^PCGs_BCscore_.*\\.txt$")
positive_genes <- fread('TF_transcripion_factors_genelist.txt', header = FALSE)  
positive_genes <- as.data.frame(positive_genes)

# 2) Run for VETFs
files <- list.files(pattern = "^PCGs_BCscore_.*\\.txt$")
positive_genes <- fread('VETFs_list.txt', header = FALSE)  
positive_genes <- as.data.frame(positive_genes)

# 3) Run for GO terms 
files <- list.files(pattern = "^PCGs_BCscore_.*\\.txt$")
positive_genes <- fread('GO_0048729_list.txt', header = FALSE)  
positive_genes <- as.data.frame(positive_genes)

# GO_0030154_list.txt
# GO_0001708_list.txt
# GO_0009887_list.txt
# GO_0048729_list.txt

# 3) Run for developmentally dynamic lncRNAs
files <- list.files(pattern = "^LncSars_BCscore_.*\\.txt$")
positive_genes <- fread('LncRNAs_dev_dynamic.txt', header = FALSE)  
positive_genes <- as.data.frame(positive_genes)

# 4) Run for evolutionarily conserved lncRNAs
files <- list.files(pattern = "^LncSars_BCscore_.*\\.txt$")
positive_genes <- fread('LncRNAs_evo_conserved.txt', header = FALSE)  
positive_genes <- as.data.frame(positive_genes)

# 4) Run for human-specific lncRNAs
files <- list.files(pattern = "^LncSars_BCscore_.*\\.txt$")
positive_genes <- fread('LncRNAs_human_specific.txt', header = FALSE)  
positive_genes <- as.data.frame(positive_genes)

# 4) Run for disease-associated lncRNAs
files <- list.files(pattern = "^Lncipeds_BCscore_.*\\.txt$")
positive_genes <- fread('LncRNAs_disease_associated.txt', header = FALSE)  
positive_genes <- as.data.frame(positive_genes)

# 5) Run for cell-type specific enhancers 
files <- list.files(pattern = "^EnhancerAtlas_BCscore_.*\\.txt$")
positive_genes <- fread('EnhancerAtlas_CTspecific_enhancers.txt', header = FALSE)
positive_genes <- as.data.frame(positive_genes)

# 6) Run for cell-type specific enhancers 
files <- list.files(pattern = "^dbSUPER_BCscore_.*\\.txt$")
positive_genes <- fread('dbSUPER_superenhancers.txt', header = FALSE)
positive_genes <- as.data.frame(positive_genes)



# 6) Run for Zoonomia alignment Phylop conserved sequences


# 7) Run for HAR sequences



###########################################################################
# Step 2) Process input HM files

### Function to process all input HM files
process_file <- function(file) {
  data <- fread(file)
  colnames(data) <- c("gene", "score")

  # Sort data by score in descending order
  data <- data[order(-data$score)]

  # Create 100 bins
  data$bin <- ntile(data$score, 100) # create 100 bins
  # data$bin <- ceiling(seq_along(data$gene) / 100) # create bins of 100 genes
  # data$single_bin <- 1:nrow(data) # create single gene bins

  return(data)
}

processed_data <- lapply(files, process_file)
names(processed_data) <- files
# df <- processed_data[["PCGs_BCscore_H3K4me3.txt"]]

###########################################################################
# Step 3) Run FET enrichments across cumulative gene bins

### Function to cumulatively bin genes and run looped FET enrichment
all_results <- data.frame()

# Loop through each file in the processed data list
for (file_name in names(processed_data)) {
  
  df <- processed_data[[file_name]] # Read in the data

  # Initialize vectors to store p-values, OR and proportions for the current file
  p_vals <- c()
  odds_ratio <- c()
  proportions <- c()
  
  # Loop through each bin in the data frame
  for (b in 1:length(unique(df$bin))) {

  	# Cumulatively add the bins
  	cum_bins <- unique(df$bin)[1:b]

    # Identify the genes in the positive gene list
    m <- positive_genes$V1
    N <- df$gene
    
    # Identify genes not in the positive gene list
    n <- df$gene[which(df$gene %!in% positive_genes$V1)]
    
    # Identify the genes in the current bin
    k <- df[which(df$bin %in% cum_bins), ]$gene
    
    # Count the number of positive genes in the current bin
    x <- length(which(k %in% m))
    
    # Construct the contingency table
    A <- x
    B <- length(k) - x
    C <- length(m) - x
    D <- length(n) - (length(k) - x)
    
    mat <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)
    
    # Perform Fisher's Exact Test
    result <- fisher.test(mat, alternative = 'greater')
    
    # Store the p-value
    p_vals <- c(p_vals, result$p.value)

    # Store the odd's ratio
    odds_ratio <- c(odds_ratio, result$estimate)

    # Store the proportion of positive genes in bin
    proportions <- c(proportions, A/(A+C))
  }
  
  # Create results dataframe
  temp_df <- data.frame('bin' = c(1:100), 'p-value' = p_vals,'odds_ratio' = odds_ratio, 'proportion' = proportions, 'file' = file_name)
  
  # Append to the aggregated data frame
  all_results <- rbind(all_results, temp_df)
}

### Format the output

# Adjust p-values for multiple testing across all files
all_results$padj <- p.adjust(all_results$p.value, method = 'fdr')

# Round the OR for categorical 
all_results$OR_rounded <- round(all_results$odds_ratio)

# Save results tables
# write.table(all_results, "240827_CumBinned_FET_results_VETFs.txt")
# write.table(all_results, "240827_CumBinned_FET_results_DevDynLncRNAs.txt")
# write.table(all_results, "240827_CumBinned_FET_results_VETFs.txt")
# write.table(all_results, "240911_CumBinned_FET_results_CTSenhancers.txt")

###########################################################################
# Step 4) Plot enrichments

### FINAL VETF ENRICHMENT FIGURE X PLOTS:

# P-vals
pval <- ggplot(all_results, aes(bin, -log10(padj), color = file)) +
  geom_point() +
  geom_line() +
  labs(x = 'Bin', y = '-log10(Adjusted P-value)', title = 'Tissue morph Cumulative bins') +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("240827_HMcomp_superenhancers_Cumulative_bins_FET_pvals.pdf", plot = pval, width = 10, height = 6)

# ggsave("240827_HMcomp_VETF_Cumulative_bins_FET_pvals.pdf", plot = pval, width = 10, height = 6)
# ggsave("240827_HMcomp_DevDynLncRNAs_Cumulative_bins_FET_pvals.pdf", plot = pval, width = 10, height = 6)

# Odd's ratio
odds <- ggplot(all_results, aes(bin, odds_ratio, color = file)) +
  geom_point() +
  geom_line() +
  labs(x = 'Bin', y = 'Odds ratio', title = 'Tissue morph Cumulative bins') +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("240827_HMcomp_superenhancers_Cumulative_bins_FET_oddsratio.pdf", plot = odds, width = 10, height = 6)

# ggsave("240827_HMcomp_VETFs_Cumulative_bins_FET_oddsratio.pdf", plot = odds, width = 10, height = 6)
# ggsave("240827_HMcomp_DevDynLncRNAs_Cumulative_bins_FET_oddsratio.pdf", plot = odds, width = 10, height = 6)


# scp uqesinni@bunya.rcc.uq.edu.au:/scratch/project/palpant_scratch/ESinniah/CellConstraint_project/ 


# scp uqesinni@bunya.rcc.uq.edu.au:/scratch/project/palpant_scratch/ESinniah/CellConstraint_project/240827_HMcomp_superenhancers_Cumulative_bins_FET_oddsratio.pdf . 


# Get statistics for manuscript

subs <- subset(all_results, bin == 1 & file == "PCGs_BCscore_H3K27me3.txt")

subs <- subset(all_results, file == "dbSUPER_BCscore_H3K27me3.txt")

subs <- subset(all_results, file == "dbSUPER_BCscore_H3K27ac.txt")
subs <- subs[order(-subs$padj), ]


