#!/opt/R/4.0.5/bin/Rscript

root_dir <- "~/TCGA_SARC"
setwd(root_dir)
source("scripts/my_functions.R")

library(edgeR)
library(tidyverse)
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)
library(viridis)


# change directory
dea_dir <-  "./results/02_DEA"
if (!dir.exists(dea_dir)) {dir.create(dea_dir)}
setwd(dea_dir)

##########################
### survival all genes ###
##########################
limma.TMM <- get(load("../recount3_mesoderm_TMMvoom.RData"))
all_clin <- read_tsv("../recount3_mesoderm_metadata.tsv", col_types = cols()) %>% 
  filter(project == "SARC")
all_clin$SurvObj <- with(all_clin, Surv(new_death, death_event))
expr_data <- limma.TMM$E[,all_clin$fixed_names] %>% t %>% 
  as.data.frame %>% rename_with(~ make.names(.x))

## survival analysis
covariates <- colnames(expr_data)
univ_formulas <- sapply(covariates, function(x) as.formula(paste('all_clin$SurvObj ~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = expr_data)})

# Extract data 
univ_res_deg <- lapply(univ_models,
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
res <- as.data.frame(t(as.data.frame(univ_res_deg, check.names = FALSE)))
res <- add_column(res, gene = rownames(limma.TMM$E), var = apply(expr_data, 2, var), .before = 1)
res$p.adj <- signif(p.adjust(res$p.value, method = "fdr"), 3)

# proportional hazards assumption
res$p.zph <- unlist(lapply(univ_models, function(x) ifelse(!is.na(x[["coefficients"]]), cox.zph(x)$table[1,3], NA)))

# write survival results
write_tsv(res, "rec3_TMM_surv.tsv")


#################
### DEA/TFTEA ###
#################
annot <- read_tsv("../../data/tf-target_interactions_Wdorot.tsv", col_types = cols())

## fit
fit <- lmFit(limma.TMM, design)
fit <- contrasts.fit(fit, contrast.matrix) 
fit <- eBayes(fit)

## get results
tts_deg <- mat.or.vec(nc = 8, nr = 0)
tf_res <- mat.or.vec(nr = 0, nc = 10)
for(i in seq(1, ncol(contrast.matrix))){
  # get DEA results
  deg_df <- topTable(fit, coef = i, number = nrow(limma.TMM$E), sort.by = "t")
  deg_df <- deg_df %>% rownames_to_column("gene") %>%
    add_column(subtype = con_list$short_histo[i], .before = 1)
  tts_deg <- rbind(tts_deg, deg_df)
  # run TFTEA
  stat <- deg_df[,'t']
  names (stat) <- deg_df$gene
  res <- uvGsa (rankstat = stat, annotation = annot)
  res <- res %>% rownames_to_column("TF") %>%
    add_column(subtype = con_list$short_histo[i], .before = 1)
  tf_res <- rbind(tf_res, res)
}

# write results
write_tsv(tts_deg, "rec3_subtypes_DEA.tsv")
write_tsv(tf_res, "rec3_subtypes_TFTEA.tsv")


## aggregate DEG results and add survival data
com_deg <- rbind(tts_deg %>% filter(adj.P.Val < 0.05 & t > 0) %>% group_by(gene) %>% 
                   mutate(n = length(gene), expression = mean(AveExpr), deg_t = mean(t), 
                          deg_p = mean(adj.P.Val)) %>% ungroup %>% filter(n == max(n)) %>% 
                   dplyr::select(gene, deg_t, deg_p, expression) %>% unique,
                 tts_deg %>% filter(adj.P.Val < 0.05 & t < 0) %>% group_by(gene) %>% 
                   mutate(n = length(gene), expression = mean(AveExpr), deg_t = mean(t), 
                          deg_p = mean(adj.P.Val)) %>% ungroup %>% filter(n == max(n)) %>% 
                   dplyr::select(gene, deg_t, deg_p, expression) %>% unique) %>% 
  left_join(as.data.frame(org.Hs.egSYMBOL) %>% mutate(gene_id = as.numeric(gene_id)), 
            by = c("gene" = "gene_id")) %>% relocate(symbol, .before = 2)
write_tsv(com_deg, "rec3_subtypes_common_DEG.tsv")

dea_surv <- read_tsv("rec3_TMM_surv.tsv", col_types = cols())
deg_surv <- left_join(com_deg, dea_surv, by = "gene")
write_tsv(deg_surv, "rec3_subtypes_common_DEG_surv.tsv")

# aggregate TF results
com_tf <- rbind(tf_res %>% filter(adj.X < 0.05 & LOR.X > 0) %>% group_by(TF) %>% 
                  mutate(n = length(TF), LOR = mean(LOR.X), adj_p = mean(adj.X)) %>% ungroup %>% 
                  select(TF, n, LOR, adj_p) %>% unique %>% arrange(desc(n), adj_p),
                tf_res %>% filter(adj.X < 0.05 & LOR.X < 0) %>% group_by(TF) %>% 
                  mutate(n = length(TF), LOR = mean(LOR.X), adj_p = mean(adj.X)) %>% ungroup %>% 
                  select(TF, n, LOR, adj_p) %>% unique %>% arrange(desc(n), adj_p))
write_tsv(com_tf, "rec3_subtypes_common_TF.tsv")


# join counted TFs with annotations and DEA/surv results
tsog <- read_tsv("../../data/human_ONG-TSG_dbs.tsv", col_types = cols())
cosmic <- read_csv("../../data/COSMIC_v96_31MAY2022-cancer_gene_census.csv", col_types = cols())

tts_dea <- read_tsv("rec3_subtypes_DEA.tsv", col_types = cols()) %>%
  mutate(de_dir = ifelse(adj.P.Val > 0.05, "no_de", ifelse(t > 0, "up", "down"))) %>%
  group_by(gene) %>% mutate(up = table(de_dir)["up"], down = table(de_dir)["down"]) %>% 
  ungroup %>% .[,c(2,4,10, 11)] %>% unique
dea_res <- left_join(tts_dea, dea_surv[,c(1,4,5,9,10,11,12)], by = "gene")
colnames(dea_res)[3:ncol(dea_res)] <- paste0("dea_", colnames(dea_res)[3:ncol(dea_res)])

tf_annot <- inner_join(select(annot, TF, TF_id), com_tf, by = "TF") %>% 
  unique %>% arrange(desc(n)) %>% 
  left_join(select(tsog, TF = DB_symbol, TSG_ONG = regulation), by = "TF") %>% 
  left_join(select(cosmic, TF = `Gene Symbol`, COSMIC = `Role in Cancer`), by = "TF") %>% 
  mutate(ts = ifelse(str_detect(COSMIC, "TSG"), 1, 0), 
         og = ifelse(str_detect(COSMIC, "oncogene"), 1, 0), 
         COSMIC = ifelse(ts == 1 & og == 1, "both", ifelse(ts == 1, "TSG", ifelse(og == 1, "ONG", NA)))) %>% 
  left_join(dea_res, by = c("TF_id" = "gene")) %>% mutate_if(is.numeric, signif) %>% 
  select(-ts, -og)

write_tsv(tf_annot, "rec3_subtypes_common_TF_annot.tsv")



###################
### heatmap TFs ###
###################

tf_res <- read_tsv("rec3_subtypes_TFTEA.tsv", col_types = cols())
tf_data <- tf_res %>% select(subtype, TF, LOR.X) %>% 
  pivot_wider(names_from = subtype, values_from = LOR.X) %>% 
  column_to_rownames("TF")

tf_data_na <- tf_res %>% mutate(LOR.X = ifelse(adj.X < 0.05, LOR.X, 0)) %>% 
  select(subtype, TF, LOR.X) %>% 
  pivot_wider(names_from = subtype, values_from = LOR.X) %>% 
  column_to_rownames("TF")
tf_data_na <- tf_data_na[rowSums(tf_data_na) != 0,]

pheatmap(tf_data_na, clustering_method = "complete", 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "correlation", 
         breaks = seq(-0.7, 0.7, length.out = 101), color = turbo(100), 
         width = 10, height = 15, fontsize_row = 7, filename = "rec3_TF_heatmap_all-viridis.png")

