library(cluster)
library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(parallel)
library(inflection)

# read in gene-drug-trait trios with their approval status 
ddir ="/QRISdata/Q3164/TRIAGE/drug_king"
dt <- fread("target_indication.tsv")

# read in RTS scores for genes 
#rdir = directory with RTS scores for triage genes  
rtsg <- fread(paste0(rdir, "EnsemblID_all_genes_v2rts.txt"))

# read in predicted drug success probability from King manuscript
result <- readRDS("Model_result_success_probability.rds")
result <- arrange(result, desc(p.mean))

# join data with gene-drug-trait trios to get approval status 
resultd <- merge(result, dt, by = c("ensembl_id", "MSH"))
resultd <- arrange(resultd, desc(resultd$lApprovedUS.EU))
# remove duplicate gene associations, preferentially 
# keeping approved drugs and high success probability
resultd2 <- resultd[!duplicated(resultd$ensembl_id),]
# assign RTS scores to the genes
resultd2 <- merge(resultd, rtsg, by.x = "ensembl_id", by.y = "V1")

# arrange by RTS (decreasing)
resultd2 <- arrange(resultd2, desc(V2))
resultd2$index <- 1:nrow(resultd2)

# bin by RTS 
bins <- seq(1, nrow(resultd2), length.out = 21)[-c(1,21)] %>% round()

# Get average predicted success probability in RTS bin for approved and non-approved drugs
std <- function(x) sd(x)/sqrt(length(x))
prob_bin <- function(bin) {

	resultd2 <- resultd2[resultd2$index < bin,]

	resultd21 <- resultd2[resultd2$lApprovedUS.EU == TRUE,]
	resultd22 <- resultd2[resultd2$lApprovedUS.EU != TRUE,]
	df1 <- resultd21 %>%  summarise(mean = mean(p.mean), sem = std(p.mean), nobjects = bin)
	df2 <- resultd22 %>%  summarise(mean = mean(p.mean), sem = std(p.mean), nobjects = bin)
	df <- rbind(df1, df2)
	df$status <- c("Approved",  "Not_Approved")
	return(df)
}

results <- lapply(bins, prob_bin) %>% bind_rows()

# Visualise results
p <- ggplot(results, aes(x = nobjects, y = mean, color = status))+
geom_line(linewidth = 1.5)+
geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2, position=position_dodge(0.05))+
scale_color_manual(values = c("#C14949", "#948B89"))+
geom_point(aes(fill=factor(status)), size=3, shape=21, stroke=0)+
scale_fill_manual(values = c("#A73939", "#766D6B"))+
theme_bw()+
ylab("Average Success Probability")+
xlab("Gene Index (Ordered by RTS")

pdf("Success_probability_nodupgenedis_binned.pdf", width = 6, height = 8)
print(p)
dev.off()


# Quantify AUC for predicting approval status using predicted success probability
# across RTS bins 
AUC_separation <- function(bin) {

	temp <- resultd2[resultd2$index < bin,]

	roc_obj <- roc(temp$lApprovedUS.EU, temp$p.mean, auc = TRUE)

	return(roc_obj$auc)
}


results <- sapply(bins, AUC_separation)
results <- data.frame(auc = results, bin = bins)

# Visualise results
p <- ggplot(results, aes(x = bin, y = auc))+
geom_line(linewidth = 1.5)+
geom_point(size=3, shape=21, stroke=0)+
theme_bw()+
ylab("AUC (Separating Approved vs Non-Approved Drugs")+
xlab("Gene Index (Ordered by RTS")

pdf("AUC_across_RTSbins.pdf", width = 6, height = 4)
print(p)
dev.off()
