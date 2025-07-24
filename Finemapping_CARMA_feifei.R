
# CARMA Fine-mapping results from Fei-Fei


# HPC = Bunya
# wd = "/scratch/user/uqesinni/TRIAGEi_analysis"
# RDM: /QRISdata/Q3856/01_TRIAGE_intergenic/06-External_Collaborators/


#################################################################################################################
############################### CARMA Fine-mapping Results from Fei-Fei #########################################
#################################################################################################################

# 1) Get number of PIP >0.9 variants for 28 traits in three groups

# scp -r uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/06-External_Collaborators/finemapping /scratch/user/uqesinni/TRIAGEi_analysis/


library(dplyr)
library(stringr)

# testing:
# a <- readLines('/scratch/user/uqesinni/TRIAGEi_analysis/finemapping//CARMA_Allergic.log')
# values <- as.numeric(unlist(strsplit(a[5], "\\s+")))
# values <- values[!is.na(values)]

# file_df <- data.frame(
#  	Bina = values[1],
#     Conti = values[2],
#     Nofun = values[3],
#     Filename = basename(a))
# yy <- bind_rows(result_df, file_df)


folder_path <- "/scratch/user/uqesinni/TRIAGEi_analysis/finemapping/"
log_files <- list.files(path = folder_path, pattern = "\\.log$", full.names = TRUE) #28 traits
result_df <- data.frame()

for (x in log_files) {
	log_content <- readLines(x) # Read the contents of the log file
  
	values <- as.numeric(unlist(strsplit(log_content[5], "\\s+")))
	values <- values[!is.na(values)] # Extract the values for "bina," "conti," and "nofun"

	file_df <- data.frame(
    Bina = values[1],
    Conti = values[2],
    Nofun = values[3],
    Trait = basename(x)  # Extract the filename
  )
    result_df <- bind_rows(result_df, file_df)
}

print(result_df)

result_df$Trait <- gsub(".log", "", result_df$Trait)
result_df <- result_df[-c(22, 28), ] # remove rows with NAs

write.table(result_df, '231012_CARMA_finemapping_results.txt', quote = F, col.names = T, row.names = F, sep = '\t') 

# Paired t-test

bina_vs_nofun <- t.test(result_df$Bina, result_df$Nofun, paired = TRUE)  # 0.008243
# 	Paired t-test
# data:  result_df$Bina and result_df$Nofun
# t = 2.8695, df = 25, p-value = 0.008243
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#   5.753938 35.015293
# sample estimates:
# mean difference 
#        20.38462 

conti_vs_nofun <- t.test(result_df$Conti, result_df$Nofun, paired = TRUE) # 0.009503
#Paired t-test

# data:  result_df$Conti and result_df$Nofun
# t = 2.8092, df = 25, p-value = 0.009503
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#   6.558365 42.595481
# sample estimates:
# mean difference 
#        24.57692 

bina_vs_conti <- t.test(result_df$Bina, result_df$Conti, paired = TRUE) # 0.02001
# 	Paired t-test

# data:  result_df$Bina and result_df$Conti
# t = -2.485, df = 25, p-value = 0.02001
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#  -7.6668982 -0.7177172
# sample estimates:
# mean difference 
#       -4.192308 

###########################################
# Barplot
library(ggplot2)
library(tidyr)
library(data.table)

result_df <- fread('231012_CARMA_finemapping_results.txt')

result_df_long <- pivot_longer(result_df, cols = c(Conti,Bina, Nofun), names_to = "Variable", values_to = "Count")

result <- ggplot(result_df_long, aes(x = Trait, y = Count, fill = Variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  labs(x = "Trait", y = "Number of PIP>0.9 variants", fill = "Variable") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file="231023_CARMA_finemapping_alltraits_PIP_barplot.pdf")
plot(result)
dev.off()

#############################
## Bubble plot

# # Normalize the Count values within each Trait
# result_df$Sum <- rowSums(result_df[, c("Bina", "Conti", "Nofun")])

# # Calculate the percentages within each Trait
# result_df$Bina_Percentage <- (result_df$Bina / result_df$Sum) * 100
# result_df$Conti_Percentage <- (result_df$Conti / result_df$Sum) * 100
# result_df$Nofun_Percentage <- (result_df$Nofun / result_df$Sum) * 100

# perc_long <- pivot_longer(result_df, cols = c(Bina_Percentage,Conti_Percentage, Nofun_Percentage), names_to = "Variable", values_to = "Count")

# Normalize the data to "Nofun"
norm <- result_df
norm$Bina <- norm$Bina / norm$Nofun
norm$Conti <- norm$Conti / norm$Nofun
norm$Nofun <- norm$Nofun / norm$Nofun

norm_long <- pivot_longer(norm, cols = c(Conti,Bina, Nofun), names_to = "Variable", values_to = "Count")

bubble_plot <- ggplot(norm_long, aes(x = Trait, y = Variable, size = Count)) +
  geom_point() +
  scale_size_continuous(range = c(5, 15)) +  
  labs(x = "Trait", y = "Variable", size = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(bubble_plot)


# library(viridis)
# library(wesanderson)
# pal <- wes_palette("Zissou1", 100, type = "continuous")


bubble_plot <- ggplot(norm_long, aes(x = Trait, y = Variable, size = Count, color = Count)) +
  geom_point() +
  scale_size_continuous(range = c(5, 15)) +  
  scale_color_gradient(low = "sky blue", high = "salmon2") +  
  labs(x = "Trait", y = "Variable", size = "Count", color = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(bubble_plot)

pdf(file="231023_CARMA_finemapping_alltraits_PIP_bubbleplot.pdf")
plot(bubble_plot)
dev.off()


# scp uqesinni@bunya.rcc.uq.edu.au:/scratch/user/uqesinni/TRIAGEi_analysis/231023_CARMA_finemapping_alltraits_PIP_barplot.pdf .

##################################################

# Boxplot

library(dplyr)

sem_data <- norm_long %>%
  group_by(Variable) %>%
  summarize(mean_count = mean(Count), sem_count = sd(Count) / sqrt(n()))

boxplot_plot <- ggplot(result_df_long, aes(x = Variable, y = Count)) +
  geom_boxplot(width = 0.4, fill = "lightblue", color = "blue", outlier.shape = NA) +
  labs(x = "Variable", y = "Mean Count") +
  geom_jitter(shape=16, position=position_jitter(0.2))
  theme_minimal()

plot(boxplot_plot)


# Remove HT outlier

no_ht <- result_df[-13, ]
no_ht_long <- pivot_longer(no_ht, cols = c(Conti,Bina, Nofun), names_to = "Variable", values_to = "Count")

boxplot_plot <- ggplot(no_ht_long, aes(x = Variable, y = Count)) +
  geom_boxplot(width = 0.4, fill = "lightblue", color = "blue", outlier.shape = NA) +
  labs(x = "Variable", y = "Mean Count") +
  geom_jitter(shape=16, position=position_jitter(0.2))
  theme_minimal()

pdf(file="231023_CARMA_finemapping_alltraits_counts_boxplot.pdf")
plot(boxplot_plot)
dev.off()

# scp uqesinni@bunya.rcc.uq.edu.au:/scratch/user/uqesinni/TRIAGEi_analysis/231023_CARMA_finemapping_alltraits_CDFplot.pdf .

############################################

# Number of high PIP variants

library(reshape2)
library(ggplot2)

# Melt the data frame to create a PIP column and an Annotation column
melted_data <- melt(result_df, id.vars = "Trait", measure.vars = c("Bina", "Conti", "Nofun"))

# Create the CDF plot
cdf_plot <- ggplot(data = melted_data, aes(x = value, color = variable)) +
  stat_ecdf() +
  xlab("PIP Score") +
  ylab("Cumulative Probability") +
  ggtitle("CDF of PIP Scores by Annotation") +
  scale_color_manual(values = c("Bina" = "red", "Conti" = "blue", "Nofun" = "green"))

# Show the plot
print(cdf_plot)

#######
#Remove outlier

# Melt the data frame to create a PIP column and an Annotation column
no_ht_melt <- melt(no_ht, id.vars = "Trait", measure.vars = c("Bina", "Conti", "Nofun"))

# Create the CDF plot
cdf_plot <- ggplot(data = no_ht_melt, aes(x = value, color = variable)) +
  stat_ecdf() +
  xlab("PIP Score") +
  ylab("Cumulative Probability") +
  ggtitle("CDF of PIP Scores by Annotation") +
  scale_color_manual(values = c("Bina" = "red", "Conti" = "blue", "Nofun" = "green"))

# Show the plot
pdf(file="231023_CARMA_finemapping_alltraits_CDFplot.pdf")
plot(cdf_plot)
dev.off()


# Number of high PIP variants called
sum(no_ht$Bina)
# [1] 1487
sum(no_ht$Conti)
# [1] 1553
sum(no_ht$Nofun)
# [1] 1139




###############################################################################################
###############################################################################################

# 05/02/24 

# New CDF proportions for all PIP scores


# 1) Get PIP variant info for 28 traits in three groups

# scp -r uqesinni@bunya.rcc.uq.edu.au:/QRISdata/Q3856/01_TRIAGE_intergenic/06-External_Collaborators/finemapping /scratch/user/uqesinni/TRIAGEi_analysis/

library(data.table)
library(dplyr)
library(stringr)

# testing:
# a <- fread('/scratch/user/uqesinni/TRIAGEi_analysis/finemapping/CARMA_WHRadjBMI_agg.txt')

folder_path <- "/scratch/user/uqesinni/TRIAGEi_analysis/finemapping/"
file_list <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE) #28 traits
combined_results <- data.frame()

for (file in file_list) {
    dt <- fread(file)
    dt[, V1 := paste0(rsid, "_", tools::file_path_sans_ext(basename(file)))]  # Create a new column 'V1' by combining 'rsid' and the basename
    dt <- dt[, . (V1, type, PIP)]
    combined_results <- bind_rows(combined_results, dt)  
}

head(combined_results)
dim(combined_results)
# [1] 149529309         3

write.table(combined_results, '240205_CARMA_finemapping_results_PIPscores.txt', quote = F, col.names = T, row.names = F, sep = '\t') 


# CDF plot for PIP scores

library(reshape2)
library(ggplot2)
library(data.table)

combined_results <- fread('240205_CARMA_finemapping_results_PIPscores.txt')

# Create the CDF plot
cdf_plot <- ggplot(data = combined_results, aes(x = PIP, color = type)) +
  stat_ecdf() +
  xlab("PIP Score") +
  ylab("Cumulative Probability") +
  ggtitle("CDF of PIP Scores by Annotation") +
  scale_color_manual(values = c("bina" = "red", "conti" = "blue", "nofun" = "green"))

# Show the plot
print(cdf_plot)

########################

#Subset

subs <- combined_results[combined_results$PIP > 0.1, ]
subs <- subs[order(subs$PIP), ]

cdf_plot <- ggplot(data = subs, aes(x = PIP, color = type)) +
  stat_ecdf() +
  xlab("PIP Score") +
  ylab("Cumulative Probability") +
  ggtitle("CDF of PIP Scores by Annotation") +
  scale_y_continuous(limits = c(0.5, 1)) 

# Show the plot
print(cdf_plot)

pdf(file="240205_CARMA_finemapping_alltraits_allPIPscores_CDFplot.pdf")
plot(cdf_plot)
dev.off()

# scp -r uqesinni@bunya.rcc.uq.edu.au:/scratch/user/uqesinni/TRIAGEi_analysis/240205_CARMA_finemapping_alltraits_allPIPscores_CDFplot.pdf .

# Statistics

# Extract PIP values for 'nofun' and 'bina'
nofun_group <- combined_results$PIP[combined_results$type == 'nofun']
bina_group <- combined_results$PIP[combined_results$type == 'bina']
conti_group <- combined_results$PIP[combined_results$type == 'conti']

# one-way Kolmogorov-Smirnov (KS) test

ks_result <- ks.test(nofun_group, bina_group)
#   Asymptotic two-sample Kolmogorov-Smirnov test

# data:  nofun_group and bina_group
# D = 0.10135, p-value < 2.2e-16
# alternative hypothesis: two-sided

ks_result <- ks.test(nofun_group, conti_group)
#   Asymptotic two-sample Kolmogorov-Smirnov test

# data:  nofun_group and conti_group
# D = 0.12389, p-value < 2.2e-16
# alternative hypothesis: two-sided










































###################################################################

# Testing space:


library(data.table)
library(stats) 
library(ggplot2)


data <- fread('/scratch/user/uqesinni/TRIAGEi_analysis/finemapping/CARMA_AMena_agg.txt.gz')

# a <- fread('/QRISdata/Q3856/01_TRIAGE_intergenic/06-External_Collaborators/finemapping/CARMA_Asthma_agg.txt.gz')

dim(data)
# [1] 5610171      10

summary(data$PIP)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00e+00 7.10e-06 1.34e-05 2.10e-04 2.36e-05 1.00e+00


################
# subs <- data[1:100, ]
# spip <- subs$PIP

# cdf_proportions <- lapply(subs$PIP, function(x) {
#   sum(subs$PIP <= x) / length(subs$PIP)
# })

# xx <- data.frame(PIP = subs$PIP, cdf_proportion = unlist(cdf_proportions))

# plot(xx, xlab="Data", ylab="Cumulative Probability",  
#      main="Empirical Cumulative Distribution Function")

###################

cdf_props <- lapply(data$PIP, function(x) {
  sum(data$PIP <= x) / length(data$PIP)
})

outs <- data.frame(PIP = data$PIP, cdf_proportion = unlist(cdf_props))


plot(xx, xlab="Data", ylab="Cumulative Probability",  
     main="Empirical Cumulative Distribution Function")


#############


nofun <- subset(data, (data$type == "nofun"))
bina <- subset(data, (data$type == "bina"))
conti <- subset(data, (data$type == "conti"))

cdf_props <- lapply(conti$PIP, function(x) {
  sum(conti$PIP <= x) / length(conti$PIP)
})

outs <- data.frame(PIP = conti$PIP, cdf_proportion = unlist(cdf_props))






subs <- subset(data, select=c(PIP, type))

prop.table(table(subs$type))
#      bina     conti     nofun 
# 0.3333333 0.3333333 0.3333333

table(data[data$PIP >0.9, "type"])
# type
#  bina conti nofun 
#    25    27    21 

PIP <- data$PIP
PIP <- sort(PIP)

ecdf_vals <- ecdf(PIP)

plot(ecdf_vals, xlab="Data", ylab="Cumulative Probability",  
     main="Empirical Cumulative Distribution Function")



plot(subs, xlab="Data", ylab="Cumulative Probability")



# a <- seq(min(PIP), max(PIP), length = 100)

cdf_function <- function(x) { 
    # Computes prob. for a single value
    mean(PIP <= x)
}

cdf_values <- sapply(PIP, cdf_function)

plot(PIP, cdf_values)






ecdf_function <- ecdf(pip_scores)

cdf_values <- ecdf_function(pip_values)








PIP <- data$PIP
#  summary(PIP)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00e+00 7.10e-06 1.34e-05 2.10e-04 2.36e-05 1.00e+00 

ecdf_func <- ecdf(PIP)
# summary(ecdf_func)
# Empirical CDF:	  3560399 unique values with summary
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.0000072 0.0000134 0.0002113 0.0000237 1.0000000



PIP_props <- cumsum(PIP) / sum(PIP)

ecdf_func <- ecdf(PIP_props)

plot(ecdf_func, xlab="Data", ylab="Cumulative Probability",  
     main="Empirical Cumulative Distribution Function")











