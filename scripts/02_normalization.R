#!/opt/R/4.0.5/bin/Rscript

root_dir <- "~/TCGA_SARC"
setwd(root_dir)
source("scripts/my_functions.R")

library(org.Hs.eg.db)
library(edgeR)
library(tidyverse)
library(Rtsne)
library(pheatmap)
library(RColorBrewer)

setwd("results")

#####################
### normalization ###
#####################
## read data
pheno <- read_tsv("recount3_mesoderm_metadata.tsv", col_types = cols())
expr <- read.table("recount3_mesoderm_expr.tsv", sep = "\t", row.names = 1, header = T)
expr <- translate_ensemble(expr)

# first approach, direct normalization with no grouping factors
expr.TMM <- DGEList(counts = expr)
expr.TMM <- calcNormFactors(expr.TMM,  method = "TMM")

# get log counts with voom
limma.TMM <- voom(expr.TMM)
save(limma.TMM, file = "recount3_mesoderm_TMMvoom.RData")


# second, normalization with gene filtering, for DEGs
expr.TMM <- DGEList(counts = expr, group = factor(pheno$project))
keep <- filterByExpr(expr.TMM)
expr.TMM <- expr.TMM[keep,,keep.lib.sizes=FALSE]
expr.TMM <- calcNormFactors(expr.TMM,  method = "TMM")

# prepare design and run voom
tissue <- factor(pheno$short_histo)
source <- factor(make.names(pheno$source))
design <- model.matrix(~ 0 + tissue + source)
colnames(design)[1:length(levels(tissue))] <- levels(tissue)
limma.TMM <- voom(expr.TMM, design = design)
save(limma.TMM, file = "recount3_mesoderm_TMMvoomFilt.RData")



##################
### some plots ###
##################
if (!dir.exists("./01_figures/data")) {dir.create("./01_figures/data", recursive = T)}
setwd("01_figures")

pheno <- read_tsv("../recount3_mesoderm_metadata.tsv", col_types = cols())
limma.TMM <- get(load("../recount3_mesoderm_TMMvoom.RData"))
limma.TMM <- limma.TMM$E

# distribution of normalized counts
png("rec3_TMM_box_counts.png", width = 550, height = 500)
boxplot(limma.TMM[,1:50], ylab = "norm counts", main = "Distribution of TMM norm counts")
dev.off()

### Apply z-scaling to data and draw heatmap and tSNE
z_data <- as.data.frame(scal(limma.TMM))
z_data[is.na(z_data)] <- 0

save(z_data, file = "data/TMM-z_data.Rdata")

# heatmap new classification
pheno2 <- pheno %>% mutate(gtex_tis = ifelse(project == "GTEx", tissue, NA), 
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

# heatmap of all genes, samples in TCGA SARC paper
pheatmap(z_data[,pheno2$fixed_names], annotation_col = annotation_col, annotation_colors = annot_colors,
         show_rownames = F, show_colnames = F, clustering_method = "complete", 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
         breaks = seq(-1.5, 1.5, length.out = 101), width = 14, height = 10, 
         filename = "rec3_TMM_heatmap_TCGA.png")

gc()


# t-SNE new classification
pheno2 <- pheno %>%   mutate(gtex_tis = ifelse(project == "GTEx", tissue, NA), 
                             sarc = ifelse(project == "SARC", short_histo, NA))

tis_colors <- colorRampPalette(brewer.pal(8, "Pastel2"))(length(unique(pheno2$gtex_tis)))
names(tis_colors) <- unique(pheno2$gtex_tis)
subtype_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(pheno2$sarc)))
names(subtype_colors) <- unique(pheno2$sarc)
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
  ggtitle("tSNE plot of mesenchymal expression data")
ggsave("rec3_TMM_tSNE_TCGA.png", width = 11, height = 9)

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
  ggtitle("tSNE plot of TCGA sarcoma expression values")
ggsave("rec3_TMM_tSNE_TCGA_SARC.png", width = 6, height = 4)


## save large intermediate data
save(tSNE_fit, file = "data/TMM-tSNE_fit-z_data_Rtsne.Rdata")


