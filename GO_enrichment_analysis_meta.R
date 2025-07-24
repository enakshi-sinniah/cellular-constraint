###################################################################################################################
################################## Histone Modification Comparisons- Big GO enrichment analysis ###############################################
###################################################################################################################

########################
##### GO analysis ###### 
########################

### Step 1) Get list of all GO terms and their metadata

library(GO.db)
go <- keys(GO.db, keytype="GOID") # list of all go terms
GO_all <- select(GO.db, columns=c("GOID","TERM","ONTOLOGY","DEFINITION"), keys=go, keytype="GOID")
dim(GO_all)
# [1] 42888     2
write.table(GO_all, 'GO_terms_meta.txt', quote = F, col.names = T, row.names = F, sep = '\t') 

###################################################################
### Step 2) Find GO terms of interest based on keywords

library(data.table)
GO_all <- fread('GO_terms_meta.txt')

# Search for keywords in GO metadata file
keep_terms <- grepl("cell fate specification|cell development|development|differentiation", GO_all$TERM, ignore.case = TRUE)

remove_terms <- !grepl("regulation", GO_all$TERM, ignore.case = TRUE)

subs <- GO_all[keep_terms & remove_terms, ]
subs <- as.data.table(subs)

# subs <- GO_all[grep("cell fate specification|cell development|morphogenesis|development|differentiation", GO_all$TERM, ignore.case = TRUE), ]
# subs <- GO_all[grep("cell fate specification", GO_all$TERM, ignore.case = TRUE), ]

# Get GO associated genes for list of GO terms of interest
library(org.Hs.eg.db)
go_terms <- subs$GOID
all_go_terms <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys= go_terms, columns= c("ENSEMBL","GENENAME","SYMBOL"))

all_go_terms <- as.data.table(all_go_terms)
gots <- all_go_terms[!is.na(SYMBOL), .(GOALL, SYMBOL)]
gots <- unique(gots)

# Filter to only keep GO terms with atleast 100 genes associated
counts <- gots[, .N, by = GOALL]
keep <- counts[N >= 100, GOALL]
gots <- gots[GOALL %in% keep, .(GOALL, SYMBOL)]

# Manual check for filter of unwanted or redundant GO terms
gots_ids <- as.character(unique(gots$GOALL))
subs[, GOID := as.character(GOID)]
temp <- subs[GOID %in% gots_ids, .(GOID, TERM)]
temp <- temp[order(TERM)]

###################################################################
### Step 3) Run GO enrichment analysis for GO terms of interest

# Process input HM files

files <- list.files(pattern = "^PCGs_BCscore_.*\\.txt$") 

### Function to process all input HM files
process_file <- function(file) {
  data <- fread(file)
  colnames(data) <- c("gene", "score")
  # Sort data by score in descending order
  data <- data[order(-data$score)]
  # Annotate top 200 genes
  data$Top_200 <- ifelse(1:nrow(data) <= 200, 'Y', 'N')
  return(data)
}
processed_data <- lapply(files, process_file)
names(processed_data) <- files
# df <- processed_data[["PCGs_BCscore_H3K27me3.txt"]]

# Run FET enrichments

### Function to run looped FET enrichment
all_results <- data.frame()

'%!in%' <- function(x,y)!('%in%'(x,y))

# Loop through each file in the processed data list
for (goall in unique(gots$GOALL)) {
  # Filter genes for the current GOALL
  go_genes <- as.data.frame(gots[GOALL == goall, SYMBOL])
  colnames(go_genes) <- "V1"
  # go_genes <- as.data.frame(gots[GOALL == "GO:0001501", SYMBOL])

  # Loop through each file 
  for (file_name in names(processed_data)) {
    
    df <- processed_data[[file_name]] # Read in the data
    
    # Initialize vectors to store results
    p_vals <- c()
    odds_ratio <- c()
    proportions <- c()

    # Create contingency matrix
    N <- df$gene
    GO_pos <- go_genes$V1[go_genes$V1 %in% df$gene]
    GO_neg <- df$gene[!df$gene %in% go_genes$V1]
    Top_pos <- df$gene[df$Top_200 == 'Y']
    Top_neg <- df$gene[df$Top_200 == 'N']
    
    A <- length(which(Top_pos %in% GO_pos))
    B <- length(which(Top_pos %!in% GO_pos))
    C <- length(which(Top_neg %in% GO_pos))
    D <- length(which(Top_neg %!in% GO_pos))

    mat <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)
    
    # Perform Fisher's Exact Test
    result <- fisher.test(mat, alternative = 'greater')
    
    # Store results
    p_val <- result$p.value
    odds_ratio_val <- result$estimate
    proportion <- A / (A + C)
    
    # Create results dataframe for the current GOALL and file
    temp_df <- data.frame('GOALL' = goall, 'p-value' = p_val, 'odds_ratio' = odds_ratio_val, 'proportion' = proportion, 'file' = file_name)
    
    # Append to the aggregated data frame
    all_results <- rbind(all_results, temp_df)
    rownames(all_results) <- NULL
  }
}

# Format the output
all_results$HM_id <- sub(".*re_(.*)\\.txt", "\\1", all_results$file)
all_results$padj <- p.adjust(all_results$p.value, method = 'fdr')
all_results$OR_rounded <- round(all_results$odds_ratio)
all_results$log_pvalue <- -log10(all_results$p.value)
all_results$significant <- ifelse(all_results$p.value <= 0.01, 'S', 'NS')

# Get GOid names
GO_all <- fread('GO_terms_meta.txt')
all_results <- merge(all_results, GO_all[, c("GOID", "TERM")], by.x = "GOALL", by.y="GOID", keep.all.x = TRUE)

length(which(all_results$significant == "S" & all_results$HM_id == "H3K27me3"))
# 119 # H3K27me3
# 57 # H3K4me3

length(unique(all_results$GOALL))
# [1] 147

########################################################################################
### Step 4) Plot heatmap of enrichments

# Heatmap of meta GO-term enrichment 

library(ggplot2)

# Alphabetically order Y-axis top down
all_results$TERM <- factor(all_results$TERM, levels = rev(sort(unique(all_results$TERM))))

# Create the heatmap
heaty <- ggplot(all_results, aes(x = HM_id, y = TERM, fill = log_pvalue)) +
  geom_tile(color = "white") + 
  scale_fill_gradientn(colors = c("White", "#C7433A")) + 
  labs(title = "GO enrichment heatmap",
       x = "HM",
       y = "GO TERM",
       fill = "-log10(p-value)") +
  theme_minimal() +
  geom_text(data = all_results[all_results$significant == 'S', ],
            aes(label = "*"), color = "black", size = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("240916_HMcomp_GOterms_meta_heatmap.pdf", plot = heaty, width = 10, height = 6)

# scp uqesinni@bunya.rcc.uq.edu.au:/scratch/project/palpant_scratch/ESinniah/CellConstraint_project/240916_HMcomp_GOterms_meta_heatmap.pdf . 






















