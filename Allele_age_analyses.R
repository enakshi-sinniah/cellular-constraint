library(cluster)
library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(parallel)
library(inflection)

### DATA PRE-PROCESSING:

### Compile all individual chromosome files (atlas.chr*.csv.gz) 
# and merge with triage RTS score

### Load in pre-processed data:
#V2 = RTS score column 
allele_age <- readRDS("Allele_age_compiled_triage.RDS")

### Retrieve minor allele frequencies for all common variants 
allsnp <- fread("snp150Common.txt.gz")
# remove excess information 
maf <- allsnp[,c(2,3,4,5,10,25)]
# only keep variants that have a minor allele frequency and are present in the allele age database
maf <- maf[maf$V5 %in% allele_age$V1 & maf$V25 != "",]
# format minor allele frequency 
maf2 <- str_split_fixed(maf$V25, ",", 3)
maf$maf <- as.numeric(maf2[,2])
colnames(maf) <- c("CHR", "POS1", "POS2", "RSID", "Alleles", "MAF", "maf")

# remove variants with no minor allele freq/non common variants 
allele_age <- allele_age[allele_age$V1 %in% maf$RSID,]
# join datasets
comb <- inner_join(allele_age, maf, by = c("V1" = "RSID"))

# remove variants with more > 1 age estimation (preferentially keeping Combined)
comb <- arrange(comb, DataSource)
comb <- comb[!duplicated(comb$V1),]

# identiy priority variants:
comb <- arrange(comb, desc(V2))
comb$index <- 1:nrow(comb)
# get inflection point of priority:
ede(comb$index, comb$V2, index = 1)
# assign priority vs non-priority
comb$priority <- ifelse(comb$index < 45707, "P", "NP")


##############################################################################
############# ANALYSE DATA ###################################################
##############################################################################

# Break data into bins based on minor allele frequency to account for MAF effect on allele age:
MAF_bins <- seq(min(comb$maf), max(comb$maf), length.out = 30)

bin_age <- function(pos) {
	# restrict data to variants in the desired MAF window
	combt <- comb[comb$maf > MAF_bins[pos] & comb$maf < MAF_bins[pos+1],]
	# compare average mutational allele age in priority vs non-priority 
	df <- combt %>% group_by(priority) %>% summarise(mean = mean(AgeMean_Mut), sd = sd(AgeMean_Mut), bin = pos, maf = mean(maf))
	return(df)
}


results <- lapply(1:100, bin_age) %>% bind_rows()
results$priority <- factor(results$priority, levels = c("P", "NP"))

# Visalise results: 
p <- ggplot(results, aes(x = maf, y = mean, color = priority))+
geom_smooth()+
scale_color_manual(values = c("#C14949", "#948B89"))+
geom_point(aes(fill=factor(priority)), size=3, shape=21, stroke=0)+
scale_fill_manual(values = c("#A73939", "#766D6B"))+
theme_bw()+
ylab("Average Allele Age")+
xlab("Allele Frequency")

pdf("Allele_age_acrossMAF_comparison.pdf", width = 8, height = 7)
print(p)
dev.off()


## Run bootstrapping to check significance of difference:
bootstrap <- function(pos) {
	# get MAF window
	combt <- comb[comb$maf > MAF_bins[pos] & comb$maf < MAF_bins[pos+1],]
	# how many priority variants in window
	sn <- sum(combt$priority == "P")
	all <- vector()
	print(pos)
	# bootstrap 200 times
	for (i in 1:200) {
		# sample # of priority variants from non priority & get mean value 
		all[i] <- mean(combt[combt$V1 %in% sample(combt[combt$priority == "NP"]$V1, sn),]$AgeMean_Mut)
	}
	df <- data.frame(meanexp = mean(all), sdexp = sd(all), mean = mean(combt[combt$priority == "P",]$AgeMean_Mut),
					 bin = pos, maf = mean(combt$maf), sample = sn)
	return(df)


}

MAF_bins <- seq(min(comb$maf), max(comb$maf), length.out = 51)
results <- mclapply(1:50, bootstrap, mc.cores = 10) %>% bind_rows()

# get z scores + p-values from bootstrapping:
results$zscore <- (results$mean - results$meanexp)/results$sdexp
results$pval <- pnorm(results$mean, mean = results$meanexp, sd = results$sdexp, lower.tail = FALSE)
results$pval <- ifelse(results$zscore > 0, pnorm(results$mean, mean = results$meanexp, sd = results$sdexp, lower.tail = FALSE), 
	pnorm(results$mean, mean = results$meanexp, sd = results$sdexp, lower.tail = TRUE))

# Visualise Results
p <- ggplot(results, aes(x = maf, y = zscore))+
geom_smooth()+
geom_point()+
theme_bw()+
geom_hline(yintercept=1.65, linetype="dashed", color = "red", size=0.5)+ # at 0.05
xlab("Bootstrapped Z score")+
ylab("Allele Frequency")


pdf("Allele_age_acrossMAF_significance.pdf", width = 8, height = 7)
print(p)
dev.off()






