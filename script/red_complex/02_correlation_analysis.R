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
biom_filter_PG_TF <- t(biom_filter)[,c(1,2)]
biom_filter_PG_TD <- t(biom_filter)[,c(1,3)]
biom_filter_TF_TD <- t(biom_filter)[,c(2,3)]
range(biom_filter_PG_TF_clean$Porphyromonas_A.gingivalis)
colnames(metadata)
#PG_TF
biom_filter_PG_TF_clean <- biom_filter_PG_TF[rowSums(biom_filter_PG_TF != 0) > 0, ]
biom_filter_PG_TF_clean <- data.frame(biom_filter_PG_TF_clean) %>% tibble::rownames_to_column(var = "sample_id") %>% left_join(metadata %>% select(sample_id,continent_group,bioproject_id),by = "sample_id")
colnames(biom_filter_PG_TF_clean)
biom_filter_PG_TF_clean <- biom_filter_PG_TF_clean %>%
  mutate(State_Region = factor(continent_group,
                               levels = c(sort(setdiff(unique(continent_group), "Others")), "Others")))
table(biom_filter_PG_TF_clean$continent_group)
p_PG_TF <- ggplot(biom_filter_PG_TF_clean, aes(x = Porphyromonas_A.gingivalis, y = Tannerella.forsythia,color = State_Region)) +
  geom_point()+ theme_bw() +
  geom_smooth(method = "lm", color = "gray") +
  stat_cor(
    method = "spearman",           
    label.x.npc = "left",          
    label.y.npc = "top",           
    cor.coef.name = "spearman_R",         
    size = 5,                      
    color = "black"                
  ) +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) +
  scale_color_brewer(palette = "Set3", name = "State Region")
p_PG_TF
ggsave("p_PG_TF_spearman.pdf",p_PG_TF,width = 5.5,height = 4)
p_PG_TF1 <- ggplot(biom_filter_PG_TF_clean, aes(x = Porphyromonas_A.gingivalis, y = Tannerella.forsythia,color = State_Region)) +
  geom_point()+ theme_bw() +
  geom_smooth(method = "lm", color = "gray") +
  stat_cor(
    method = "spearman",           
    label.x.npc = "left",          
    label.y.npc = "top",           
    cor.coef.name = "spearman_R",         
    size = 3,                      
    color = "black"                
  ) +
  facet_wrap(~ State_Region) +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),legend.position = "none") +
  scale_color_brewer(palette = "Set3", name = "State Region")
p_PG_TF1
ggsave("p_PG_TF_spearman_facet_wrap.pdf",p_PG_TF1,width = 8,height = 8)
#biom_filter_PG_TD
colnames(metadata)
biom_filter_PG_TD_clean <- biom_filter_PG_TD[rowSums(biom_filter_PG_TD != 0) > 0, ]
biom_filter_PG_TD_clean <- data.frame(biom_filter_PG_TD_clean) %>% tibble::rownames_to_column(var = "sample_id") %>% left_join(metadata %>% select(sample_id,continent_group,bioproject_id),by = "sample_id")
colnames(biom_filter_PG_TD_clean)
biom_filter_PG_TD_clean <- biom_filter_PG_TD_clean %>%
  mutate(State_Region = factor(continent_group,
                               levels = c(sort(setdiff(unique(continent_group), "Others")), "Others")))
p_PG_TD <- ggplot(biom_filter_PG_TD_clean, aes(x = Porphyromonas_A.gingivalis, y = Treponema_B_986141.denticola,color = State_Region)) +
  geom_point()+ theme_bw() +
  geom_smooth(method = "lm", color = "gray") +
  stat_cor(
    method = "spearman",           
    label.x.npc = "left",          
    label.y.npc = "top",           
    cor.coef.name = "spearman_R",         
    size = 5,                      
    color = "black"                
  ) +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) +
  scale_color_brewer(palette = "Set3", name = "State Region")
p_PG_TD
ggsave("p_PG_TD_spearman.pdf",p_PG_TD,width = 5.5,height = 4)
p_PG_TD1 <- ggplot(biom_filter_PG_TD_clean, aes(x = Porphyromonas_A.gingivalis, y = Treponema_B_986141.denticola,color = State_Region)) +
  geom_point()+ theme_bw() +
  geom_smooth(method = "lm", color = "gray") +
  stat_cor(
    method = "spearman",           
    label.x.npc = "left",          
    label.y.npc = "top",           
    cor.coef.name = "spearman_R",         
    size = 3,                      
    color = "black"                
  ) +
  facet_wrap(~ State_Region) +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),legend.position = "none") +
  scale_color_brewer(palette = "Set3", name = "State Region")
p_PG_TD1
ggsave("p_PG_TD_spearman_facet_wrap.pdf",p_PG_TD1,width = 8,height = 8)
#biom_filter_TF_TD
biom_filter_TF_TD_clean <- biom_filter_TF_TD[rowSums(biom_filter_TF_TD != 0) > 0, ]
biom_filter_TF_TD_clean <- data.frame(biom_filter_TF_TD_clean) %>% tibble::rownames_to_column(var = "sample_id") %>% left_join(metadata %>% select(sample_id,continent_group,bioproject_id),by = "sample_id")
biom_filter_TF_TD_clean <- biom_filter_TF_TD_clean %>%
  mutate(State_Region = factor(continent_group,
                               levels = c(sort(setdiff(unique(continent_group), "Others")), "Others")))
p_TF_TD <- ggplot(biom_filter_TF_TD_clean, aes(x = Tannerella.forsythia, y = Treponema_B_986141.denticola,color = State_Region)) +
  geom_point()+ theme_bw() +
  geom_smooth(method = "lm", color = "gray") +
  stat_cor(
    method = "spearman",           
    label.x.npc = "left",         
    label.y.npc = "top",           
    cor.coef.name = "spearman_R",         
    size = 5,                      
    color = "black"             
  ) +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) +
  scale_color_brewer(palette = "Set3", name = "State Region")
p_TF_TD
ggsave("p_TF_TD_spearman.pdf",p_TF_TD,width = 5.5,height = 4)
p_TF_TD1 <- ggplot(biom_filter_TF_TD_clean, aes(x = Tannerella.forsythia, y = Treponema_B_986141.denticola,color = State_Region)) +
  geom_point()+ theme_bw() +
  geom_smooth(method = "lm", color = "gray") +
  stat_cor(
    method = "spearman",           
    label.x.npc = "left",          
    label.y.npc = "top",           
    cor.coef.name = "spearman_R",         
    size = 3,                      
    color = "black"                
  ) +
  facet_wrap(~ State_Region) +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),legend.position = "none") +
  scale_color_brewer(palette = "Set3", name = "State Region")
p_TF_TD1
ggsave("p_TF_TD_spearman_facet_wrap.pdf",p_TF_TD1,width = 8,height = 8)
