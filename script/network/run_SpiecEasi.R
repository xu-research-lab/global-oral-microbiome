library(SpiecEasi)
library(biomformat)
library(phyloseq)
args = commandArgs(trailingOnly=TRUE)
study_id <- args[1]
meta <- read.csv("../../data/network/metadata_hc_saliva.csv", sep=",", row.names = 2)
input_biom <- biomformat::read_biom(paste0("/beegfs/dongbiao/oral/", study_id, "/results/extracted-feature-table/feature-table.biom"))
amgut1_filt <- as.matrix(biom_data(input_biom))

inter_sid <- intersect(colnames(amgut1_filt), rownames(meta))
amgut1_filt <- amgut1_filt[, inter_sid]
amgut1_filt <- amgut1_filt[apply(amgut1_filt, 1, function(x) sum(x > 0) >= 5), ]

sparcc_amgut <- sparcc(t(amgut1_filt), iter = 50, inner_iter = 10, th = 0.01)
sparcc_cor <- as.data.frame(sparcc_amgut$Cor)
rownames(sparcc_cor) <- rownames(amgut1_filt)
colnames(sparcc_cor) <- rownames(amgut1_filt)

write.csv(sparcc_cor, paste0("../../data/network/sparcc_hc_saliva/", study_id, ".csv"))

