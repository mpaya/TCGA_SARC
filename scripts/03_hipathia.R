#!/opt/R/4.0.5/bin/Rscript

root_dir <- "~/TCGA_SARC"
setwd(root_dir)
source("scripts/my_functions.R")

library(hipathia)
library(org.Hs.eg.db)
library(tidyverse)
library(Rtsne)
library(pheatmap)
library(survival)
library(survminer)
library(RColorBrewer)

setwd("results")

## load data
pheno <- read_tsv("recount3_mesoderm_metadata.tsv", col_types = cols())
limma.TMM <- get(load("recount3_mesoderm_TMMvoom.RData"))
limma.TMM <- limma.TMM$E


## run hipathia
path_list <- read_tsv( "../data/physiological_paths.tsv", col_types = cols(), col_names = FALSE)
pathways <- load_pathways(species = "hsa", path_list[[2]])

cnt <- limma.TMM %>% normalize_data(.)
results <- hipathia(cnt, pathways)
path_vals <- get_paths_data(results, matrix = T)
path_vals <- normalize_paths(path_vals, pathways)

write.table(path_vals, file = "recount3_mesoderm_TMM_hipathia.tsv", sep = "\t", 
            col.names = T, row.names = T)
save(results, file = "rec3_TMM_hip_results.Rdata")


##################
### some plots ###
##################
setwd("01_figures")

pheno <- read_tsv("../recount3_mesoderm_metadata.tsv", col_types = cols())
path_vals <- read.table("../recount3_mesoderm_TMM_hipathia.tsv", sep = "\t", 
                        header = TRUE, row.names = 1)

# distribution of activity values
png("rec3_TMM_hipathia_box.png", width = 550, height = 500)
boxplot(path_vals[,1:50], ylab = "norm counts", main = "Distribution of TMM activity values")
dev.off()

### Apply z-scaling to data and draw heatmap and tSNE
z_data <- as.data.frame(scal(path_vals))
z_data[is.na(z_data)] <- 0
z_data <- z_data[!apply(z_data, 1, function(x) any(is.infinite(x))),]

## save large intermediate data
save(z_data, file = "data/TMM-hipathia-z_data.Rdata")


# heatmap new classification
pheno2 <- pheno %>% filter(!is.na(short_histo)) %>% 
  mutate(gtex_tis = ifelse(project == "GTEx", tissue, NA), 
         sarc = ifelse(project == "SARC", short_histo, NA))
annotation_col = data.frame(
  Subtype = factor(pheno2$sarc),
  Tissue = factor(pheno2$gtex_tis),
  Type = factor(pheno2$project)
)
rownames(annotation_col) = pheno2$fixed_names

getColors <- colorRampPalette(brewer.pal(12, "Paired"))
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
type_colors <- getColors(length(levels(annotation_col$Type)))
names(type_colors) <- levels(annotation_col$Type)
tis_colors <- getColors(length(levels(annotation_col$Tissue)))
names(tis_colors) <- levels(annotation_col$Tissue)
subtype_colors <- colorBlindBlack8[1:length(levels(annotation_col$Subtype))]
names(subtype_colors) <- levels(annotation_col$Subtype)
annot_colors <- list(Type = type_colors, Tissue = tis_colors, Subtype = subtype_colors)

# heatmap of all circuits
pheatmap(z_data[,pheno2$fixed_names], annotation_col = annotation_col, annotation_colors =annot_colors,
         show_rownames = F, show_colnames = F, clustering_method = "complete", 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "correlation", 
         breaks = seq(-2, 2, length.out = 101), width = 14, height = 10, 
         filename = "rec3_TMM_hipathia_heatmap_TCGA.png")

gc()


# t-SNE
pheno2 <- pheno %>% mutate(gtex_tis = ifelse(project == "GTEx", tissue, NA), 
                           sarc = ifelse(project == "SARC", short_histo, NA))

tis_colors <- colorRampPalette(brewer.pal(8, "Pastel2"))(length(unique(na.omit(pheno2$gtex_tis))))
names(tis_colors) <- unique(na.omit(pheno2$gtex_tis))
subtype_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(na.omit(pheno2$sarc))))
names(subtype_colors) <- unique(na.omit(pheno2$sarc))
annot_colors <- c(tis_colors, subtype_colors)

set.seed(142)
tSNE_fit <- z_data %>% t %>% Rtsne(check_duplicates = F)
tSNE_df <- tSNE_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1 = "V1",
         tSNE2 = "V2") %>%
  cbind(pheno)

tSNE_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = short_histo,
             shape = short_histo)) +
  geom_point() + 
  scale_shape_manual(values = seq(1, 18)) + 
  scale_color_manual(values = annot_colors[levels(factor(tSNE_df$short_histo))]) + 
  theme_bw() + 
  ggtitle("tSNE plot of mesenchymal activity values")
ggsave("rec3_TMM_hipathia_tSNE_TCGA.png", width = 11, height = 9)

## tSNE for sarcoma samples
pheno2 <- filter(pheno2, project == "SARC")
set.seed(142)
tSNE_fitTCGA <- z_data[,pheno2$fixed_names] %>% t %>% Rtsne(check_duplicates = F)
tSNE_dfTCGA <- tSNE_fitTCGA$Y %>% 
  as.data.frame() %>%
  rename(tSNE1 = "V1",
         tSNE2 = "V2") %>%
  cbind(pheno2)

tSNE_dfTCGA %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = short_histo,
             shape = short_histo)) +
  geom_point() + 
  scale_shape_manual(values = seq(1, 18)) + 
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(20)) + 
  theme_bw() + 
  ggtitle("tSNE plot of TCGA sarcoma activity values")
ggsave("rec3_TMM_hipathia_tSNE_TCGA_SARC.png", width = 6, height = 4)

