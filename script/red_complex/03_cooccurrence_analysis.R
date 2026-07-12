library(qiime2R)
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)
library(readxl)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
table_biom <- read_qza("/home/yemai/Project_Xu/red_complex/species-table-rel-new.qza")
table_biom_data <- table_biom$data
rm(table_biom)
colSums(table_biom_data)
metadata <- read_tsv("/home/yemai/Project_Xu/red_complex/metadata_oral_qc_dedup.tsv")
metadata$continent_group <- replace_na(metadata$continent_group, "Others")
metadata <- metadata[metadata$sample_id %in% colnames(table_biom_data),]
colnames(metadata)
table_biom_data <- table_biom_data[,metadata$sample_id]
original_rownames <- rownames(table_biom_data)
new_rownames <- sub(".*;s__", "", original_rownames)
empty_indices <- which(new_rownames == "" | is.na(new_rownames))
if (length(empty_indices) > 0) {
  new_rownames[empty_indices] <- paste0("OTU", empty_indices)
}
rownames(table_biom_data) <- new_rownames
target_species <- c("Porphyromonas_A gingivalis", "Tannerella forsythia", "Treponema_B_986141 denticola")
biom_filter <- table_biom_data[target_species, ]
biom_filter_binary <- ifelse(biom_filter > 0, 1, 0)
head(sort(table(metadata$bioproject_id), decreasing = TRUE), 10)

#PRJNA801882(5672)
biom_filter_binary_P_observed_PRJNA801882 <- biom_filter_binary[,metadata$bioproject_id%in%"PRJNA801882"]
rowSums(biom_filter_binary_P_observed_PRJNA801882)
P_observed_PRJNA801882 <- mean(colSums(biom_filter_binary_P_observed_PRJNA801882) == 3)                
species_probs <- rowMeans(biom_filter_binary_P_observed_PRJNA801882)
print(round(species_probs, 4))
P_expected_independent <- prod(species_probs)
n_samples <- ncol(biom_filter_binary_P_observed_PRJNA801882)  
new_table <- matrix(0, nrow = 3, ncol = n_samples)
rownames(new_table) <- rownames(biom_filter_binary_P_observed_PRJNA801882)
colnames(new_table) <- colnames(biom_filter_binary_P_observed_PRJNA801882)
set.seed(123)
for (i in rownames(new_table)) {
  original_row <- biom_filter_binary_P_observed_PRJNA801882[i, ]
  new_table[i, ] <- sample(original_row)
}

rowMeans(new_table)  
P1 <- mean(colSums(new_table) == 3)
n_permutations <- 1000

perm_results_PRJNA801882 <- data.frame(
  Times = paste0("P", 1:n_permutations),  
  Probs = numeric(n_permutations),        
  stringsAsFactors = FALSE
)
set.seed(123)
pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)

for (i in 1:n_permutations) {
  permuted_matrix <- t(apply(biom_filter_binary_P_observed_PRJNA801882, 1, sample))
  perm_results_PRJNA801882$Probs[i] <- mean(colSums(permuted_matrix) == 3)
  setTxtProgressBar(pb, i)
}
z_score_PRJNA801882 <- (P_observed_PRJNA801882 - mean(perm_results_PRJNA801882$Probs)) / sd(perm_results_PRJNA801882$Probs)
p_value_two_sided <- 2 * pnorm(-abs(z_score_PRJNA801882))
library(ggplot2)
p <- ggplot(perm_results_PRJNA801882,aes(x = Probs)) +
  geom_histogram(bins = 30,fill = "steelblue",color = "white",alpha = 0.7) +
  geom_vline(xintercept = P_observed_PRJNA801882,color = "red",linewidth = 1.5,linetype = "solid") +
  geom_vline(xintercept = mean(perm_results_PRJNA801882$Probs),color = "black",linewidth = 1,linetype = "dashed") +
  annotate("text",
           x = P_observed_PRJNA801882-0.02,
           y = Inf,
           label = paste("Z-score =", round(z_score_PRJNA801882, 3),"\n","p = ",format(p_value_two_sided, scientific = TRUE, digits = 3)),
           vjust = 2,
           size = 3,
           fontface = "bold") +
  labs(title = "PRJNA801882(5672)_Lung_Cancer_oral_wash\nNorth America",x = "Three-Species Coexistence Probability", y = "Frequency") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 10)
  )
p
ggsave("frequency_PRJNA801882.pdf",p,width = 4,height = 4)



#PRJEB38062(178)
biom_filter_binary_P_observed_PRJEB38062 <- biom_filter_binary[,metadata$bioproject_id%in%"PRJEB38062"]
rowSums(biom_filter_binary_P_observed_PRJEB38062)
P_observed_PRJEB38062 <- mean(colSums(biom_filter_binary_P_observed_PRJEB38062) == 3)                
species_probs <- rowMeans(biom_filter_binary_P_observed_PRJEB38062)
print(round(species_probs, 4))
P_expected_independent <- prod(species_probs)
n_samples <- ncol(biom_filter_binary_P_observed_PRJEB38062)  
new_table <- matrix(0, nrow = 3, ncol = n_samples)
rownames(new_table) <- rownames(biom_filter_binary_P_observed_PRJEB38062)
colnames(new_table) <- colnames(biom_filter_binary_P_observed_PRJEB38062)
set.seed(123)
for (i in rownames(new_table)) {
  original_row <- biom_filter_binary_P_observed_PRJEB38062[i, ]
  new_table[i, ] <- sample(original_row)
}

rowMeans(new_table)  
P1 <- mean(colSums(new_table) == 3)
n_permutations <- 1000

perm_results_PRJEB38062 <- data.frame(
  Times = paste0("P", 1:n_permutations),  
  Probs = numeric(n_permutations),       
  stringsAsFactors = FALSE
)
set.seed(123)
cat("正在进行", n_permutations, "次随机置换...\n")
pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)

for (i in 1:n_permutations) {
  permuted_matrix <- t(apply(biom_filter_binary_P_observed_PRJEB38062, 1, sample))
  perm_results_PRJEB38062$Probs[i] <- mean(colSums(permuted_matrix) == 3)
  setTxtProgressBar(pb, i)
}
z_score_PRJEB38062 <- (P_observed_PRJEB38062 - mean(perm_results_PRJEB38062$Probs)) / sd(perm_results_PRJEB38062$Probs)
p_value_two_sided <- 2 * pnorm(-abs(z_score_PRJEB38062))
library(ggplot2)
p <- ggplot(perm_results_PRJEB38062,aes(x = Probs)) +
  geom_histogram(bins = 30,fill = "steelblue",color = "white",alpha = 0.7) +
  geom_vline(xintercept = P_observed_PRJEB38062,color = "red",linewidth = 1.5,linetype = "solid") +
  geom_vline(xintercept = mean(perm_results_PRJEB38062$Probs),color = "black",linewidth = 1,linetype = "dashed") +
  annotate("text",
           x = P_observed_PRJEB38062-0.01,
           y = Inf,
           label = paste("Z-score =", round(z_score_PRJEB38062, 3),"\n","p = ",format(p_value_two_sided, scientific = TRUE, digits = 3)),
           vjust = 2,
           size = 3,
           fontface = "bold") +
  labs(title = "PRJEB38062(178)_Oral_Mucositis_Saliva\nAsia",x = "Three-Species Coexistence Probability", y = "Frequency") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 10)
  )
p
ggsave("frequency_PRJEB38062.pdf",p,width = 4,height = 4)
table(metadata$Healthy_Status[metadata$BioProject%in%"PRJEB38062"])


#PRJEB19122(22) 
biom_filter_binary_P_observed_PRJEB19122 <- biom_filter_binary[,metadata$bioproject_id%in%"PRJEB19122"]
rowSums(biom_filter_binary_P_observed_PRJEB19122)
P_observed_PRJEB19122 <- mean(colSums(biom_filter_binary_P_observed_PRJEB19122) == 3)                
species_probs <- rowMeans(biom_filter_binary_P_observed_PRJEB19122)
print(round(species_probs, 4))
P_expected_independent <- prod(species_probs)
n_samples <- ncol(biom_filter_binary_P_observed_PRJEB19122)  
new_table <- matrix(0, nrow = 3, ncol = n_samples)
rownames(new_table) <- rownames(biom_filter_binary_P_observed_PRJEB19122)
colnames(new_table) <- colnames(biom_filter_binary_P_observed_PRJEB19122)
set.seed(123)
for (i in rownames(new_table)) {
  original_row <- biom_filter_binary_P_observed_PRJEB19122[i, ]
  new_table[i, ] <- sample(original_row)
}

rowMeans(new_table)  
P1 <- mean(colSums(new_table) == 3)
n_permutations <- 1000

perm_results_PRJEB19122 <- data.frame(
  Times = paste0("P", 1:n_permutations),  
  Probs = numeric(n_permutations),       
  stringsAsFactors = FALSE
)
set.seed(123)
cat("正在进行", n_permutations, "次随机置换...\n")
pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)

for (i in 1:n_permutations) {
  permuted_matrix <- t(apply(biom_filter_binary_P_observed_PRJEB19122, 1, sample))
  perm_results_PRJEB19122$Probs[i] <- mean(colSums(permuted_matrix) == 3)
  setTxtProgressBar(pb, i)
}
z_score_PRJEB19122 <- (P_observed_PRJEB19122 - mean(perm_results_PRJEB19122$Probs)) / sd(perm_results_PRJEB19122$Probs)
p_value_two_sided <- 2 * pnorm(-abs(z_score_PRJEB19122))
library(ggplot2)
p <- ggplot(perm_results_PRJEB19122,aes(x = Probs)) +
  geom_histogram(bins = 30,fill = "steelblue",color = "white",alpha = 0.7) +
  geom_vline(xintercept = P_observed_PRJEB19122,color = "red",linewidth = 1.5,linetype = "solid") +
  geom_vline(xintercept = mean(perm_results_PRJEB19122$Probs),color = "black",linewidth = 1,linetype = "dashed") +
  annotate("text",
           x = P_observed_PRJEB19122-0.05,
           y = Inf,
           label = paste("Z-score =", round(z_score_PRJEB19122, 3),"\n","p = ",format(p_value_two_sided, scientific = TRUE, digits = 3)),
           vjust = 2,
           size = 3,
           fontface = "bold") +
  labs(title = "PRJEB19122(22)_Periodontal_Disease_Tooth\nNorthern America",x = "Three-Species Coexistence Probability", y = "Frequency") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 10)
  )
p
ggsave("frequency_PRJEB19122.pdf",p,width = 4,height = 4)




#PRJNA866676(49)
biom_filter_binary_P_observed_PRJNA866676 <- biom_filter_binary[,metadata$bioproject_id%in%"PRJNA866676"]
rowSums(biom_filter_binary_P_observed_PRJNA866676)
P_observed_PRJNA866676 <- mean(colSums(biom_filter_binary_P_observed_PRJNA866676) == 3)                
species_probs <- rowMeans(biom_filter_binary_P_observed_PRJNA866676)
print(round(species_probs, 4))
P_expected_independent <- prod(species_probs)
n_samples <- ncol(biom_filter_binary_P_observed_PRJNA866676) 
new_table <- matrix(0, nrow = 3, ncol = n_samples)
rownames(new_table) <- rownames(biom_filter_binary_P_observed_PRJNA866676)
colnames(new_table) <- colnames(biom_filter_binary_P_observed_PRJNA866676)
set.seed(123)
for (i in rownames(new_table)) {
  original_row <- biom_filter_binary_P_observed_PRJNA866676[i, ]
  new_table[i, ] <- sample(original_row)
}


rowMeans(new_table)  
P1 <- mean(colSums(new_table) == 3)
n_permutations <- 1000

perm_results_PRJNA866676 <- data.frame(
  Times = paste0("P", 1:n_permutations),  
  Probs = numeric(n_permutations),        
  stringsAsFactors = FALSE
)
set.seed(123)


pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)


for (i in 1:n_permutations) {
  permuted_matrix <- t(apply(biom_filter_binary_P_observed_PRJNA866676, 1, sample))
  perm_results_PRJNA866676$Probs[i] <- mean(colSums(permuted_matrix) == 3)
  setTxtProgressBar(pb, i)
}
z_score_PRJNA866676 <- (P_observed_PRJNA866676 - mean(perm_results_PRJNA866676$Probs)) / sd(perm_results_PRJNA866676$Probs)
p_value_two_sided <- 2 * pnorm(-abs(z_score_PRJNA866676))
library(ggplot2)
p <- ggplot(perm_results_PRJNA866676,aes(x = Probs)) +
  geom_histogram(bins = 20,fill = "steelblue",color = "white",alpha = 0.7) +
  geom_vline(xintercept = P_observed_PRJNA866676,color = "red",linewidth = 1.5,linetype = "solid") +
  geom_vline(xintercept = mean(perm_results_PRJNA866676$Probs),color = "black",linewidth = 1,linetype = "dashed") +
  annotate("text",
           x = P_observed_PRJNA866676-0.02,
           y = Inf,
           label = paste("Z-score =", round(z_score_PRJNA866676, 3),"\n","p = ",format(p_value_two_sided, scientific = TRUE, digits = 3)),
           vjust = 2,
           #hjust = ifelse(P_observed > mean(perm_results$Probs), 1.2, -0.2),
           #color = "red",
           size = 3,
           fontface = "bold") +
  labs(title = "PRJNA866676(49)_Oral_OSCC\nAsia",x = "Three-Species Coexistence Probability", y = "Frequency") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 10)
  )
p
ggsave("frequency_PRJNA866676.pdf",p,width = 4,height = 4)



#PRJEB70783
biom_filter_binary_P_observed_PRJEB70783 <- biom_filter_binary[,metadata$bioproject_id%in%"PRJEB70783"]
rowSums(biom_filter_binary_P_observed_PRJEB70783)
P_observed_PRJEB70783 <- mean(colSums(biom_filter_binary_P_observed_PRJEB70783) == 3)                
species_probs <- rowMeans(biom_filter_binary_P_observed_PRJEB70783)
print(round(species_probs, 4))
P_expected_independent <- prod(species_probs)
n_samples <- ncol(biom_filter_binary_P_observed_PRJEB70783) 
new_table <- matrix(0, nrow = 3, ncol = n_samples)
rownames(new_table) <- rownames(biom_filter_binary_P_observed_PRJEB70783)
colnames(new_table) <- colnames(biom_filter_binary_P_observed_PRJEB70783)
set.seed(123)
for (i in rownames(new_table)) {
  original_row <- biom_filter_binary_P_observed_PRJEB70783[i, ]
  new_table[i, ] <- sample(original_row)
}

rowMeans(new_table)  
P1 <- mean(colSums(new_table) == 3)
n_permutations <- 1000

perm_results_PRJEB70783 <- data.frame(
  Times = paste0("P", 1:n_permutations),  
  Probs = numeric(n_permutations),     
  stringsAsFactors = FALSE
)
set.seed(123)

cat("正在进行", n_permutations, "次随机置换...\n")
pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)


for (i in 1:n_permutations) {
  permuted_matrix <- t(apply(biom_filter_binary_P_observed_PRJEB70783, 1, sample))
  perm_results_PRJEB70783$Probs[i] <- mean(colSums(permuted_matrix) == 3)
  setTxtProgressBar(pb, i)
}
z_score_PRJEB70783 <- (P_observed_PRJEB70783 - mean(perm_results_PRJEB70783$Probs)) / sd(perm_results_PRJEB70783$Probs)
p_value_two_sided <- 2 * pnorm(-abs(z_score_PRJEB70783))
library(ggplot2)
p <- ggplot(perm_results_PRJEB70783,aes(x = Probs)) +
  geom_histogram(bins = 20,fill = "steelblue",color = "white",alpha = 0.7) +
  geom_vline(xintercept = P_observed_PRJEB70783,color = "red",linewidth = 1.5,linetype = "solid") +
  geom_vline(xintercept = mean(perm_results_PRJEB70783$Probs),color = "black",linewidth = 1,linetype = "dashed") +
  annotate("text",
           x = P_observed_PRJEB70783-0.02,
           y = Inf,
           label = paste("Z-score =", round(z_score_PRJEB70783, 3),"\n","p = ",format(p_value_two_sided, scientific = TRUE, digits = 3)),
           vjust = 2,
           #hjust = ifelse(P_observed > mean(perm_results$Probs), 1.2, -0.2),
           #color = "red",
           size = 3,
           fontface = "bold") +
  labs(title = "PRJEB70783(220)_Oral_disease_Saliva\nNorthern America",x = "Three-Species Coexistence Probability", y = "Frequency") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 10)
  )
p
ggsave("frequency_PRJEB70783.pdf",p,width = 4,height = 4)
