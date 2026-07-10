library(qiime2R)
library(readr)
library(dplyr)
library(ggplot2)
table_biom <- read_qza("/home/yemai/Project_Xu/red_complex/species-counts.qza")
metadata <- read_tsv("/home/yemai/Project_Xu/red_complex/metadata_oral_qc_dedup.tsv")
table(metadata$continent_group)
table_biom_data <- table_biom$data
dim(table_biom_data)
table_biom_data[1:4,1:4]
metadata <- metadata[metadata$sample_id %in% colnames(table_biom_data),]
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
species_count <- colSums(biom_filter > 0)
none_count <- sum(species_count == 0)      
single_count <- sum(species_count == 1)    
double_count <- sum(species_count == 2)    
triple_count <- sum(species_count == 3)    
results <- data.frame(
  Category = c("Triple (All three)", "Double (Any two)", "Single (Only one)", "None"),
  Count = c(triple_count, double_count, single_count, none_count),
  Percentage = round(c(triple_count, double_count, single_count, none_count) / 
                       ncol(biom_filter) * 100, 2)
)
results <- results %>%
  mutate(number = case_when(
    Category == "Single (Only one)" ~ "Only one",
    Category == "Double (Any two)" ~ "Any two",
    Category == "Triple (All three)" ~ "All three",
    Category == "None" ~ "None",
  ))
results$number <- factor(results$number,levels=c("All three","Any two","Only one","None"))
p <- ggplot(results, aes(x = number, y = Count, fill = number)) +
  geom_bar(stat = "identity", width = 0.7) + theme_classic() +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
            vjust = -0.5, size = 4, fontface = "bold") +
  ylim(0, max(results$Count) * 1.1) +
  labs(title = "Red Complex in Oral Samples",
       x = "",
       y = "Number of Samples") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "none")
p
ggsave("statistics/three_bacteria.pdf",p,width = 6,height = 5)

presence <- biom_filter > 0
true_counts <- colSums(presence)
three_true <- sum(true_counts == 3)   
two_true <- sum(true_counts == 2)     
one_true <- sum(true_counts == 1)     
zero_true <- sum(true_counts == 0) 

two_true_samples <- colSums(presence) == 2
biom_two_true <- biom_filter[, two_true_samples]
presence_two_true <- presence[, two_true_samples]
pairwise_stats <- data.frame(
  Pair = c("Pg_Tf", "Pg_Td", "Tf_Td"),
  Count = c(
    sum(presence_two_true[1, ] & presence_two_true[2, ]),
    sum(presence_two_true[1, ] & presence_two_true[3, ]),
    sum(presence_two_true[2, ] & presence_two_true[3, ])
  ),
  Percentage = NA  
)
pairwise_stats$Percentage <- round(pairwise_stats$Count / ncol(biom_filter) * 100, 1)
p1 <- ggplot(pairwise_stats, aes(x = Pair, y = Count, fill = Pair)) +
  geom_bar(stat = "identity", width = 0.7) + theme_classic() +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
            vjust = -0.5, size = 4, fontface = "bold") +
  ylim(0, max(pairwise_stats$Count) * 1.2) +
  labs(title = "",
       x = "",
       y = "Number of Samples") +
  theme(plot.title = element_text(hjust = 0, size = 16, face = "bold"),
        legend.position = "none")
p1
ggsave("statistics/two_bacteria.pdf",p1,width = 6,height = 3.5)

one_true_samples <- colSums(presence) == 1
biom_one_true <- biom_filter[, one_true_samples]
presence_one_true <- presence[, one_true_samples]
single_stats <- data.frame(
  Bacteria = c("P.gingivalis", "T.forsythia", "T.denticola"),
  Count = c(
    sum(presence[1, ] & !presence[2, ] & !presence[3, ]),
    sum(!presence[1, ] & presence[2, ] & !presence[3, ]),
    sum(!presence[1, ] & !presence[2, ] & presence[3, ])
  ),
  Percentage = NA
)

single_stats$Percentage <- round(single_stats$Count / ncol(biom_filter) * 100, 1)
p2 <- ggplot(single_stats, aes(x = Bacteria, y = Count, fill = Bacteria)) +
  geom_bar(stat = "identity", width = 0.7) + theme_classic() +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
            vjust = -0.5, size = 4, fontface = "bold") +
  ylim(0, max(single_stats$Count) * 1.2) +
  labs(title = "",
       x = "",
       y = "Number of Samples") +
  theme(plot.title = element_text(hjust = -0.5, size = 16, face = "bold"),
        legend.position = "none")
p2
ggsave("statistics/one_bacteria.pdf",p2,width = 6,height = 3.5)
