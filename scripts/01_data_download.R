#!/opt/R/4.0.5/bin/Rscript

library("recount3")
library(TCGAbiolinks)
library(TCGAutils)
library(tidyverse)

source("TCGAbiolinks_functions.R")
root_dir <- "~/TCGA_SARC"
setwd(root_dir)
if (!dir.exists("results")){dir.create("results")}

################
### get data ###
################
# load list of tissue-sarcoma comparisons
tis_list <- read_tsv("data/tcga-gtex_relation.tsv", col_types = cols())


# get GTEx data
gtex_tis <- toupper(unique(tis_list$smts)) %>% str_replace(" ", "_")
human_projects <- available_projects()

gtex_pheno <- data.frame()
gtex_expr <- data.frame()
idx <- 1
for (tis in gtex_tis) {
  proj_info <- subset(human_projects, project == tis)
  rse_gene_GTEx <- create_rse(proj_info)
  gtex_pheno <- rbind(gtex_pheno, as.data.frame(colData(rse_gene_GTEx)))
  if (idx == 1) {gtex_expr <- rbind(gtex_expr, as.data.frame(assay(rse_gene_GTEx)))}
  if (idx > 1) {gtex_expr <- cbind(gtex_expr, as.data.frame(assay(rse_gene_GTEx)))}
  idx <- idx + 1
}

## write full table
gtex_pheno <- gtex_pheno %>%
  inner_join(unique(dplyr::select(tis_list, smts, smtsd, gtex)), 
             by = c("gtex.smts" = "smts", "gtex.smtsd" = "smtsd")) %>% 
  select_if(~!is.list(.)) %>% 
  .[,apply(., 2, function(x) !isTRUE(unique(is.na(x))))]
write_tsv(gtex_pheno, "results/recount3_GTEx_metadata.tsv")

## selected pheno data
gtex_pheno <- gtex_pheno %>% 
  dplyr::select(sample = external_id, rail_id, tissue = gtex, 
                tis_desc = gtex, source = gtex.smcenter, gtex = gtex.smtsd) %>% 
  add_column(project = "GTEx", fixed_names = make.names(.$sample), .after = 1)

## get expression for selected samples
gtex_expr <- gtex_expr[,gtex_pheno$sample]


# get TCGA data
human_projects <- available_projects()
proj_info <- subset(human_projects, project == "SARC")
rse_gene_SARC <- create_rse(proj_info)

## write full table
clinical <- as.data.frame(colData(rse_gene_SARC)) %>% 
  mutate(sample = rownames(.)) %>%
  select_if(~!is.list(.)) %>% 
  .[,apply(., 2, function(x) !isTRUE(unique(is.na(x))))]
write_tsv(clinical, "results/recount3_SARC_clinical.tsv")

## prepare data for survival analysis
# get values of days to death
ind_keep <- grep('days_to_death',colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_follow', colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum (is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

# and put everything together
all_clin <- data.frame(death_collapsed, fl_collapsed) %>% 
  mutate_all(function(x) as.numeric(as.character(x)))
colnames(all_clin) <- c('death_days', 'followUp_days')

# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(all_clin$death_days)){
  all_clin$new_death[i] <- ifelse(is.na(all_clin$death_days[i]),
                                  all_clin$followUp_days[i], all_clin$death_days[i])
}
all_clin$death_event <- ifelse(clinical$tcga.gdc_cases.diagnoses.vital_status == 'alive', 0,1)

## join data with selected columns
sarc_pheno <- clinical %>% 
  inner_join(unique(dplyr::select(tis_list, tis_desc, sarcoma, sarc_abbrev, sarc_desc)), 
             by = c("tcga.xml_primary_pathology_histological_type" = "sarcoma")) %>% 
  dplyr::select(sample, rail_id, project = study, tissue = sarc_abbrev, tis_desc, 
                sarcoma = tcga.xml_primary_pathology_histological_type, 
                submitter_id = tcga.gdc_cases.submitter_id, 
                type_code = tcga.cgc_sample_sample_type_code, 
                type_name = tcga.cgc_sample_sample_type, 
                gender = tcga.gdc_cases.demographic.gender, 
                source = tcga.cgc_sample_tissue_source_site, 
                year_of_birth = tcga.gdc_cases.demographic.year_of_birth, 
                race = tcga.gdc_cases.demographic.race, 
                has_radiation = tcga.cgc_radiation_therapy_id, 
                has_chemo = tcga.cgc_drug_therapy_id, 
                drug_name = tcga.cgc_drug_therapy_drug_name, 
                therapy_type = tcga.cgc_drug_therapy_pharmaceutical_therapy_type, 
                tumor_status_pre = tcga.cgc_case_tumor_status, 
                tumor_status_fu = tcga.cgc_follow_up_tumor_status) %>% 
  mutate(has_radiation = ifelse(!is.na(has_radiation), 1, 0), 
         has_chemo = ifelse(!is.na(has_chemo), 1, 0), 
         tumor_status = ifelse(!is.na(tumor_status_fu), tumor_status_fu, tumor_status_pre)) %>% 
  add_column(fixed_names = make.names(.$sample), .after = 1) %>% 
  cbind(all_clin)

## correct and add data from TCGA reclassification
library("readxl")
file <- "data/2017 TCGA SARC study TableS1.xlsx"
rev_pheno <- read_excel(file, sheet = 2, skip = 1) %>%
  separate("TCGA barcode", into = c("submitter_id",NA,"type_code"), 
           sep = c(-3, -2), remove = F) %>% 
  mutate(type_code = as.numeric(type_code)) %>% 
  rename_all(function(x)gsub(" ", "_", x)) %>% select(-gender, -tumor_status)

sarc_pheno <- left_join(sarc_pheno, rev_pheno, by = c("submitter_id", "type_code")) %>% 
  mutate(new_histo = ifelse(short_histo == "DDLPS", "DDLS", ifelse(
    short_histo %in% c("STLMS", "ULMS"), "LMS", short_histo)))

dif_pheno <- filter(sarc_pheno, !is.na(short_histo), tissue != new_histo)
fixed_sarc <- dif_pheno$new_histo
names(fixed_sarc) <- dif_pheno$submitter_id
my_desc <- dif_pheno$tis_desc
names(my_desc) <- dif_pheno$tissue

# at this point, metadata may be filtered
sarc_pheno <- sarc_pheno %>% filter(!is.na(short_histo)) %>% mutate(tissue = ifelse(
  submitter_id %in% names(fixed_sarc), fixed_sarc[submitter_id], tissue),) %>% 
  mutate(tis_desc = ifelse(submitter_id %in% dif_pheno$submitter_id, my_desc[.$tissue], tis_desc)) %>% 
  mutate(sarcoma = ifelse(submitter_id %in% dif_pheno$submitter_id, my_desc[.$tissue], sarcoma))

sarc_pheno2 <- as.data.frame(sarc_pheno)
sarc_pheno2 <- mutate(sarc_pheno2, drug_name = tolower(drug_name))

write_tsv(sarc_pheno2, "results/recount3_sarc_pheno.tsv")


# get expression for sarcoma samples
sarc_expr <- assay(rse_gene_SARC)[,sarc_pheno$sample]


# join sarc and gtex data
pheno <- bind_rows(sarc_pheno, gtex_pheno) %>% 
  mutate(short_histo = ifelse(project == "GTEx", tissue, short_histo), 
         histology = ifelse(project == "GTEx", tissue, histology))
expr <- cbind(sarc_expr, gtex_expr)
write_tsv(pheno, "results/recount3_mesoderm_metadata.tsv")
write.table(expr, file = "results/recount3_mesoderm_expr.tsv", sep = "\t", 
            col.names = NA, row.names = T)  # expression with Ensembl IDs

