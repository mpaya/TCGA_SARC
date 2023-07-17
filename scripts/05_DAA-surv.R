#!/opt/R/4.0.5/bin/Rscript

root_dir <- "~/TCGA_SARC"
setwd(root_dir)
source("scripts/my_functions.R")

library(org.Hs.eg.db)
library(edgeR)
library(hipathia)
library(tidyverse)
library(survival)
library(survminer)
library(enrichR)



#############################
### survival all circuits ###
#############################
daa_dir <-  "./results/03_DAA"
if (!dir.exists(daa_dir)) {dir.create(daa_dir)}
setwd(daa_dir)

path_vals <- read.table("../recount3_mesoderm_TMM_hipathia.tsv", sep = "\t", 
                        header = TRUE, row.names = 1)
all_clin <- read_tsv("../recount3_mesoderm_metadata.tsv", col_types = cols()) %>% 
  filter(project == "SARC", type_code == 1)
all_clin$SurvObj <- with(all_clin, Surv(new_death, death_event))
path_data <- path_vals[,all_clin$fixed_names] %>% t %>% 
  as.data.frame %>% rename_with(~ gsub(" ", "", gsub("-", "", .x)))

## survival analysis
covariates <- colnames(path_data)
univ_formulas <- sapply(covariates, function(x) as.formula(paste('all_clin$SurvObj ~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = path_data)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value <- signif(x$wald["pvalue"], digits = 3)
                         wald.test <- signif(x$wald["test"], digits = 3)
                         beta <- signif(x$coef[1], digits = 3);  # coefficient beta
                         HR <- signif(x$coef[2], digits = 3);  # exp(beta)
                         HR.ci.lower <- signif(x$conf.int[,"lower .95"], 3)
                         HR.ci.upper <- signif(x$conf.int[,"upper .95"], 3)
                         HR2 <- paste0(HR, " (", HR.ci.lower, "-", HR.ci.upper, ")")
                         conc <- signif(x$concordance["C"], digits = 3)
                         res <- c(beta, HR2, HR, HR.ci.lower, HR.ci.upper, wald.test, conc, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "HR", "HR_CI_lower",  
                                       "HR_CI_upper","wald.test", "concordance", "p.value")
                         return(res)
                       })
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res <- add_column(res, circuit = rownames(path_vals), var = apply(path_data, 2, var), .before = 1)
res$p.adj <- signif(p.adjust(res$p.value, method = "fdr"), 3)

# proportional hazards assumption
res$p.zph <- unlist(lapply(univ_models, function(x) ifelse(!is.na(x[["coefficients"]]), cox.zph(x)$table[1,3], NA)))

# write survival results
write_tsv(res, "rec3_hipathia_surv.tsv")


### Survival plots of significant pathways
path_list <- read_tsv( "../../data/physiological_paths.tsv", 
                       col_types = cols(), col_names = FALSE)
pathways <- load_pathways(species = "hsa", path_list[[2]])

z_data <- as.data.frame(t(scal(t(path_data))))
event_data <- t(apply(z_data, 1, function(x) 
  ifelse(x > 0.5, "High", ifelse(x < -0.5, "Low", NA)))) %>% as.data.frame
event_data <- cbind(all_clin, event_data)

signif <- rownames(filter(res, p.adj < 0.05))  # 13 significant

if ( !dir.exists("survival_plots") ){dir.create("survival_plots")}
for (s in signif) {
  pathname <- get_path_names(pathways, res[s,"circuit"])
  out <- paste0("survival_plots/circ_coxsignif-", res[s,]$circuit, ".png")
  fit <- survfit(Surv(new_death, death_event) ~ eval(parse(text = s)), data = event_data)
  p <- ggsurvplot(fit, data = event_data, conf.int = TRUE, pval = TRUE, 
                  surv.median.line = "hv", legend.labs = c("High activ", "Low activ"), 
                  palette = c("#E7B800", "#2E9FDF")) + ggtitle(pathname)
  p$plot <- p$plot + theme(plot.title.position = "plot")
  ggsave(out, p$plot, width = 5.4, height = 4)
}



###########
### DAA ###
###########

setwd(root_dir)

## load data
pheno <- read_tsv("results/recount3_mesoderm_metadata.tsv", col_types = cols())
limma.TMM <- get(load("results/recount3_mesoderm_TMMvoomFilt.RData"))
path_vals <- read.table("results/recount3_mesoderm_TMM_hipathia.tsv", sep = "\t", 
                        header = TRUE, row.names = 1)

path_list <- read_tsv( "data/physiological_paths.tsv", col_types = cols(), col_names = FALSE)
pathways <- load_pathways(species = "hsa", path_list[[2]])

tis_comp <- read_tsv("data/tcga-gtex_contrasts.tsv", col_types = cols())
con_list <- tis_comp %>% group_by(short_histo) %>% 
  summarise(con = paste0("(", paste(gtex, collapse = " + "), ")/", n())) %>% 
  ungroup %>% unite(contrast, short_histo, con, sep = " - ", remove = F)

# # change directory
setwd("./results/03_DAA")

###################
### Preparation ###
###################
# prepare design
tissue <- factor(pheno$short_histo)
source <- factor(make.names(pheno$source))
design <- model.matrix(~ 0 + tissue + source)
colnames(design)[1:length(levels(tissue))] <- levels(tissue)

## prepare contrasts
contrast.matrix <- makeContrasts(contrasts = con_list$contrast, levels = design) 

## fit
fit_vals <- lmFit(path_vals, design)
fit_vals <- contrasts.fit(fit_vals, contrast.matrix) 
fit_vals <- eBayes(fit_vals)

## get results
tts_res <- mat.or.vec(nc = 8, nr = 0)
for(i in seq(1, ncol(contrast.matrix))){
  dac_df <- topTable(fit_vals, coef = i, number = nrow(path_vals))
  dac_df <- dac_df %>% rownames_to_column("circuit") %>%
    add_column(subtype = con_list$short_histo[i], .before = 1)
  tts_res <- rbind(tts_res, dac_df)
}
write_tsv(tts_res, "rec3_subtypes_DAA.tsv")


## look at common circuits
com_dacl <- rbind(tts_res %>% filter(adj.P.Val < 0.05 & t > 0) %>% group_by(circuit) %>% 
                    mutate(n = length(circuit), activity = mean(AveExpr), logFC = mean(logFC), 
                           dac_t = mean(t), dac_p = mean(adj.P.Val)) %>% ungroup %>% filter(n == max(n)) %>% 
                    dplyr::select(circuit, logFC, dac_t, dac_p, activity) %>% unique,
                  tts_res %>% filter(adj.P.Val < 0.05 & t < 0) %>% group_by(circuit) %>% 
                    mutate(n = length(circuit), activity = mean(AveExpr), logFC = mean(logFC), 
                           dac_t = mean(t), dac_p = mean(adj.P.Val)) %>% ungroup %>% filter(n == max(n)) %>% 
                    dplyr::select(circuit, logFC, dac_t, dac_p, activity) %>% unique) %>% 
  add_column(names = get_path_names(pathways, .$circuit), .after = 1)
write_tsv(com_dacl, "rec3_subtypes_common_DAC.tsv")

daa_surv <- read_tsv("rec3_hipathia_surv.tsv", col_types = cols())
dac_surv <- left_join(com_dacl, daa_surv, by = "circuit")
write_tsv(dac_surv, "rec3_subtypes_common_DAC_surv.tsv")

# aggregate DAC results by # hits
num_dac <- rbind(tts_res %>% filter(adj.P.Val < 0.05 & t > 0) %>% group_by(circuit) %>% 
                   mutate(n = length(circuit), activity = mean(AveExpr), logFC = mean(logFC), 
                          dac_t = mean(t), dac_p = mean(adj.P.Val)) %>% ungroup  %>% 
                   dplyr::select(circuit, n, logFC, dac_t, dac_p, activity) %>% unique %>% 
                   arrange(desc(n), dac_p),
                 tts_res %>% filter(adj.P.Val < 0.05 & t < 0) %>% group_by(circuit) %>% 
                   mutate(n = length(circuit), activity = mean(AveExpr), logFC = mean(logFC), 
                          dac_t = mean(t), dac_p = mean(adj.P.Val)) %>% ungroup  %>% 
                   dplyr::select(circuit, n, logFC, dac_t, dac_p, activity) %>% unique %>% 
                   arrange(desc(n), dac_p)) %>% 
  add_column(names = get_path_names(pathways, .$circuit), .after = 1)
write_tsv(num_dac, "rec3_subtypes_common_DAC_num.tsv")


######################################
### hipathia functional enrichment ###
######################################

## enrichR

phys_genes <- read_tsv("../../data/circuitGenes.tsv", col_types = "cccc")
tts_res <- read_tsv("rec3_subtypes_DAA.tsv", col_types = cols())

dac_genes <- tts_res %>% filter(adj.P.Val < 0.05) %>%
  left_join(filter(phys_genes, !is.na(eff_gene)), by = "circuit")
  
# dbs <- listEnrichrDbs()
dbs <- "GO_Biological_Process_2021"
moredbs <- c(
  "KEGG_2021_Human", "GO_Molecular_Function_2021", "GO_Biological_Process_2021", 
  "GWAS_Catalog_2019", "ClinVar_2019", "DisGeNET", "Human_Phenotype_Ontology"
)

for (sarc_sub in unique(tts_res$subtype)) {
  daa_sub <- filter(dac_genes, subtype == sarc_sub)
  enriched <- enrichr(unique(daa_sub$eff_gene), dbs)
  write_tsv(enriched$GO_Biological_Process_2021, paste0("func_annot/enrichR-table-",sarc_sub,".tsv"))
  maxchar <- max(nchar(enriched$GO_Biological_Process_2021$Term[1:25]))
  w <- maxchar * 8 + 380
  p <- plotEnrich(enriched[[1]], showTerms = 25, numChar = maxchar, y = "Count", 
                  title = paste("GO enrichment of", sarc_sub, "DACs")) + 
    theme(axis.text = element_text(size = 15), plot.title = element_text(size = 20))
  ggsave(paste0("func_annot/enrichR-DAC-",sarc_sub,".png"), p, scale = 4, 
         width = w, height = 700, units = "px")
  
  ## dot plot
  dataset <- enriched
  a <- 1
  aux<-dataset[[a]][dataset[[a]]$Adjusted.P.value<0.05,]
  genes<-sapply(strsplit(aux$Genes,split=";"),length)
  aux$Gene_number<-genes
  aux<-head(aux,25)
  aux<-aux[order(aux$Combined.Score),]
  aux$log10_Pvalue<- (-log10(aux$Adjusted.P.value))
  ytitle<-"GO term"
  aux$Term<-factor(aux$Term,levels=aux$Term[order(aux$Combined.Score)])
  p <- ggplot(aux, aes(x=Combined.Score, y= Term, size = Gene_number, fill = log10_Pvalue)) +
    geom_point(shape = 21) +
    ggtitle(paste("GO enrichment of", sarc_sub, "DACs")) +
    labs(x = "Combined score", y = ytitle) +
    scale_fill_continuous(low = "khaki1", high = "firebrick2")
  ggsave(paste0("func_annot/enrichR-DAC_eff-",sarc_sub,".png"), p, scale = 2.4, 
         width = w, height = 700, units = "px")
  
}


# enrichR in selected circuits by the number of subtypes where they are DAC
# dbs <- listEnrichrDbs()
dbs <- "GO_Biological_Process_2021"

num_dac <- read_tsv("rec3_subtypes_common_DAC_num.tsv", show_col_types = F)
phys_genes <- read_tsv("../../data/circuitGenes.tsv", col_types = "cccc")


sel_dac <- filter(num_dac, n > 4)$circuit  # 108 paths
dac_genes <- filter(phys_genes, circuit %in% sel_dac, !is.na(eff_gene))$eff_gene
enriched <- enrichr(unique(dac_genes), dbs)
write_tsv(enriched$GO_Biological_Process_2021, paste0("func_annot/enrichR-table-DAC_Ngt4.tsv"))
maxchar <- max(nchar(enriched$GO_Biological_Process_2021$Term[1:25]))
w <- maxchar * 8 + 350
## Maria's plot
dataset <- enriched
a <- 1
aux<-dataset[[a]][dataset[[a]]$Adjusted.P.value<0.05,]
genes<-sapply(strsplit(aux$Genes,split=";"),length)
aux$Gene_number<-genes
aux<-head(aux,25)
aux<-aux[order(aux$Combined.Score),]
aux$log10_Pvalue<- (-log10(aux$Adjusted.P.value))
ytitle<-"GO term"
aux$Term<-factor(aux$Term,levels=aux$Term[order(aux$Combined.Score)])
p <- ggplot(aux, aes(x=Combined.Score, y= Term, size = Gene_number, fill = log10_Pvalue)) +
  geom_point(shape = 21) +
  ggtitle(paste("GO enrichment of path effectors")) +
  labs(x = "Combined score", y = ytitle) +
  scale_fill_continuous(low = "khaki1", high = "firebrick2")
ggsave(paste0("func_annot/enrichR-DAC_eff-DAC_Ngt4.png"), p, scale = 2.4, 
       width = w, height = 700, units = "px")



########################
### heatmap of  DACs ###
########################

library(pheatmap)
library(RColorBrewer)
library(viridis)

# heatmap aggregating pathway information, percent of pathway covered

sarc_pheno <- read_tsv("../recount3_mesoderm_metadata.tsv", col_types = cols()) %>% 
  filter(project == "SARC")

tts_res <- read_tsv("rec3_subtypes_DAA.tsv", col_types = cols())
phys_genes <- read_tsv("../../data/circuitGenes.tsv", col_types = "cccc")
all_paths <- unique(select(phys_genes, pathway, circuit, circuit_name))

tts_paths <- left_join(tts_res, all_paths, by = "circuit") %>% 
  mutate(is_sig = ifelse(adj.P.Val < 0.05, 1, 0)) %>% 
  group_by(subtype, pathway) %>% 
  mutate(mean_fc = mean(logFC), n_path = n(), n_sig = sum(is_sig), 
         perc_signif = signif(n_sig/n_path*100, 3)) %>% 
  ungroup %>% filter(!is.na(pathway))

path_data <- tts_paths %>% 
  select(subtype, pathway, perc_signif) %>% unique %>% 
  pivot_wider(names_from = subtype, values_from = perc_signif) %>% 
  column_to_rownames("pathway")
named_data <- path_data
rownames(named_data) <- path_list$X1[match(rownames(named_data), path_list$X2)]

pheatmap(named_data, clustering_method = "complete", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "correlation",  
         color = inferno(100), fontsize_row = 10, fontsize_col = 15, 
         width = 8, height = 12, 
         filename = "rec3_subtype_DAC_heatmap_percPathways.png")



#################
### survival2 ###
#################
# load nodes and relationships to effectors
phys_genes <- read_tsv("../../data/circuitGenes.tsv", col_types = cols())
node_rel <- read_tsv("../../data/NodeRelation2effector.tsv", col_types = cols())
phys_genes <- select(node_rel, circuit, node, node_name, rel_eff, gene) %>%
  left_join(phys_genes, by = c("circuit", "gene")) %>% 
  relocate(pathway, circuit, circuit_name, .before = 1)

## join survival-significant results
tts_daa <- read_tsv("rec3_subtypes_DAA.tsv", col_types = cols()) %>%
  mutate(de_dir = ifelse(adj.P.Val > 0.05, "no_de", ifelse(t > 0, "up", "down"))) %>%
  group_by(circuit) %>% mutate(up = table(de_dir)["up"], down = table(de_dir)["down"])
daa_surv <- read_tsv("rec3_hipathia_surv.tsv", col_types = cols())

daa_res <- left_join(tts_daa, daa_surv, by = "circuit") %>% filter(p.adj < 0.05)
colnames(daa_res)[3:ncol(daa_res)] <- paste0("daa_", colnames(daa_res)[3:ncol(daa_res)])

## load DEA results
tts_dea <- read_tsv("../02_DEA/rec3_subtypes_DEA.tsv", col_types = cols()) %>%
  mutate(de_dir = ifelse(adj.P.Val > 0.05, "no_de", ifelse(t > 0, "up", "down"))) %>%
  group_by(gene) %>% mutate(up = table(de_dir)["up"], down = table(de_dir)["down"])
dea_surv <- read_tsv("../02_DEA/rec3_TMM_surv.tsv", col_types = cols())

dea_res <- left_join(tts_dea, dea_surv, by = "gene") %>% filter(de_dir != "no_de" | p.value < 0.05)
colnames(dea_res)[3:ncol(dea_res)] <- paste0("dea_", colnames(dea_res)[3:ncol(dea_res)])

## join and write results
full_res <- left_join(daa_res, phys_genes, by = "circuit") %>%
  left_join(dea_res, by = intersect(colnames(dea_res), colnames(.))) %>% 
  filter(!is.na(eff_gene) | !is.na(dea_t))

write_tsv(full_res, "rec3_subtypes_DAA-DEA_surv_nodes.tsv")


# combine survival data on DACs and their genes
com_paths <- c(tts_daa %>% filter(adj.P.Val < 0.05 & t > 0) %>% group_by(circuit) %>% 
                 summarise(n = length(circuit)) %>% filter(n == max(n)) %>% pull(circuit),
               tts_daa %>% filter(adj.P.Val < 0.05 & t < 0) %>% group_by(circuit) %>% 
                 summarise(n = length(circuit)) %>% filter(n == max(n)) %>% pull(circuit))


daa_res <- left_join(tts_daa, daa_surv, by = "circuit") %>% filter(circuit %in% com_paths)
colnames(daa_res)[3:ncol(daa_res)] <- paste0("daa_", colnames(daa_res)[3:ncol(daa_res)])
dea_res <- left_join(tts_dea, dea_surv, by = "gene") %>% filter(de_dir != "no_de" | p.value < 0.05)
colnames(dea_res)[3:ncol(dea_res)] <- paste0("dea_", colnames(dea_res)[3:ncol(dea_res)])

full_res <- left_join(daa_res, phys_genes, by = "circuit") %>%
  left_join(dea_res, by = intersect(colnames(dea_res), colnames(.))) %>% 
  filter(!is.na(eff_gene) | !is.na(dea_t))

write_tsv(full_res, "rec3_subtypes_comDAC-DEA_surv_nodes.tsv")


