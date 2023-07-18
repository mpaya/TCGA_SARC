# Functional analysis of Soft Tissue Sarcoma signaling pathways

Soft Tissue Sarcomas (STS) are a set of rare cancers that originate in multiple parts of the body in tissues of mesenchymal origin. 
The primary form of therapy is surgery, although neo/adjuvant therapy may help improve patient outcome. 
The availability of molecular data has allowed the discovery of specific drivers for some sarcoma subtypes and thus improve actions for targeted therapy. 
In this work, we conducted pathway analysis with the mechanistic models implemented in HiPathia to further explore the biological mechanisms driving STS in conjunction with patient survival to propose potential therapies. 
For that purpose, we analyzed bulk RNA-Seq from the TCGA_SARC project for cancer samples and GTEx for non-diseased sarcomagenic tissues. 
Raw reads from both projects were downloaded from the Recount3 resource. 
We performed three main types of analyses:
* Differential activation of signaling pathways
* Survival analysis on signaling pathways
* Transcription factor-target enrichment analysis (TFTEA)

## Data files

All the files required to reproduce the analysis are provided in the data folder.

* Relation of tissues to download from recount3.
  * File name: tcga-gtex_relation.tsv
* Updated annotations for TCGA_SARC samples.
  * File name: 2017 TCGA SARC study TableS1.xlsx
  * Source: https://doi.org/10.1016/j.cell.2017.10.014
* Set of sarcoma-normal tissue contrasts.
  * File name: tcga-gtex_contrasts.tsv
* Relation of physiological KEGG pathways for analysis with HiPathia.
  * File name: physiological_paths.tsv
* Annotation of genes involved in each KEGG pathway.
  * File name: circuitGenes.tsv
* Activating/inhibiting relatinships of genes to effector genes on each circuit.
  * File name: NodeRelation2effector.tsv
* Transcription factor-target regulons.
  * File name: tf-target_interactions_Wdorot.tsv
  * Source 1: https://www.nature.com/articles/srep39709
  * Source 2: https://doi.org/10.1158/0008-5472.CAN-17-1679
* COSMIC annotations.
  * File name: COSMIC_v96_31MAY2022-cancer_gene_census.csv
  * Source: https://cancer.sanger.ac.uk/cell_lines
* Human oncogenes and tumor suppressor genes.
  * File name: human_ONG-TSG_dbs.tsv
  * Source 1: https://ongene.bioinfo-minzhao.org/
  * Source 2: https://bioinfo.uth.edu/TSGene/
 
## Authors and contributors

- Miriam Payá-Milans <miriam.paya@juntadeandalucia.es>
- Marina Esteban-Medina <marina.esteban@juntadeandalucia.es>
- Carlos Loucera <carlos.loucera@juntadeandalucia.es>
- Joaquin Dopazo <joaquin.dopazo@juntadeandalucia.es>
- Maria Peña-Chilet <mariapch84@gmail.com>
