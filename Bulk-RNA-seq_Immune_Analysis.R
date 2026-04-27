if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GSVA", "msigdbr"), ask = FALSE)
install.packages(c("tidyestimate", "tidyverse", "ggplot2", "ggpubr", "pheatmap"))

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(tidyestimate)
library(GSVA)
library(msigdbr)
library(pheatmap)


bulk_counts <- read.table("Exp_GCT.txt", header = TRUE, row.names = 1, 
                          sep = "\t", check.names = FALSE) 


sample_info <- read.table("Info.txt", header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE)


colnames(bulk_counts) <- gsub("^X", "", colnames(bulk_counts))


common_samples <- intersect(colnames(bulk_counts), sample_info$SampleID)
bulk_counts <- bulk_counts[, common_samples]
sample_info <- sample_info[match(common_samples, sample_info$SampleID), ]

keep <- rowSums(bulk_counts >= 10) >= 3
bulk_counts_filt <- bulk_counts[keep, ]
cpm <- t(t(bulk_counts_filt) / colSums(bulk_counts_filt)) * 1e6
log2cpm <- log2(cpm + 1)

estimate_results <- log2cpm %>%
  filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = FALSE) %>%
  estimate_score(is_affymetrix = FALSE)

scores <- as.data.frame(estimate_results)
colnames(scores)[colnames(scores) == "sample"] <- "Sample"
plot_df <- merge(scores, sample_info, by.x = "Sample", by.y = "SampleID")

BiocManager::install("GSVA", ask = FALSE)
library(GSVA)
library(msigdbr)

immune_gene_sets <- msigdbr(species = "Homo sapiens", collection = "C8")
immune_list <- split(immune_gene_sets$gene_symbol, immune_gene_sets$gs_name)
immune_list <- immune_list[sapply(immune_list, length) > 10]


mat <- as.matrix(log2cpm)
param <- ssgseaParam(mat, immune_list)
ssgsea_res <- gsva(param)

ssgsea_df <- as.data.frame(t(ssgsea_res))
ssgsea_df$Sample <- rownames(ssgsea_df)

all_sets <- colnames(ssgsea_df)
grep("CD8|CD4|B_CELL|MACROPHAGE|NK|TREG|MEMORY|ACTIVATED", all_sets, value = TRUE, ignore.case = TRUE)

selected_cells <- c(
  "DURANTE_ADULT_OLFACTORY_NEUROEPITHELIUM_CD8_T_CELLS",
  "DURANTE_ADULT_OLFACTORY_NEUROEPITHELIUM_CD4_T_CELLS",
  "DURANTE_ADULT_OLFACTORY_NEUROEPITHELIUM_B_CELLS",
  "DURANTE_ADULT_OLFACTORY_NEUROEPITHELIUM_MACROPHAGES",
  "DURANTE_ADULT_OLFACTORY_NEUROEPITHELIUM_NK_CELLS",
  "HE_LIM_SUN_FETAL_LUNG_C4_TREG_CELL"
)

selected_cells <- intersect(selected_cells, colnames(ssgsea_df))

heat_data <- ssgsea_df[, selected_cells, drop = FALSE]
rownames(heat_data) <- ssgsea_df$Sample

annotation_col <- data.frame(Histology = sample_info$Histology,
                             row.names = sample_info$SampleID)

pretty_names <- c(
  "CD8 T cells",
  "CD4 T cells",
  "B cells",
  "Macrophages",
  "NK cells",
  "Treg"
)
colnames(heat_data) <- pretty_names

pheatmap::pheatmap(
  t(scale(heat_data)),
  annotation_col = annotation_col,
  main = "ssGSEA Immune Cell Enrichment by Histology",
  filename = "ssGSEA_heatmap.pdf",
  width = 10,
  height = 6
)


library(ggplot2)
library(ggpubr)
library(rstatix)


p_site <- ggplot(plot_df, aes(x = Lesion.Site, y = immune, fill = Lesion.Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  theme_minimal() +
  labs(title = "ESTIMATE Immune Score by Lesion Site", y = "Immune Score") +
  stat_compare_means(method = "kruskal.test", label.y = max(plot_df$immune) + 200)
ggsave("ESTIMATE_by_site.pdf", p_site, width = 8, height = 5)

library(ggplot2)
library(ggpubr)
library(rstatix)

kw_test <- plot_df %>% kruskal_test(immune ~ Histology)
p_value <- kw_test$p

p <- ggplot(plot_df, aes(x = Histology, y = immune, fill = Histology)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8, shape = 21, color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = "ESTIMATE Immune Score by Histological Subtype",
    subtitle = paste0("Kruskal-Wallis p = ", format(p_value, digits = 3, scientific = TRUE)),
    y = "Immune Score",
    x = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave("ESTIMATE_Immune_by_Histology.pdf", p, width = 8, height = 5)

