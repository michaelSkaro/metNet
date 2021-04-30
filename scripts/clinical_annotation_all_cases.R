# Towards annotation of all TCGA metastatic loci
library(TCGAbiolinks)
library(xml2)
library(tidyverse)

projects <- c("TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD",
              "TCGA-DLBC","TCGA-ESCA","TCGA-GBM","TCGA-HNSC","TCGA-KICH","TCGA-KIRC",
              "TCGA-KIRP","TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC",
              "TCGA-MESO","TCGA-OV","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ",
              "TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THCA","TCGA-THYM",
              "TCGA-UCEC","TCGA-UCS","TCGA-UVM")

setwd("/mnt/storage/mskaro1/Clinical_annotation")

proj <- projects[2]


ACC <- function(proj){
    dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), header = TRUE,na.strings=c("","NA"))
    dat1 <- dat %>%
      dplyr::select(fileID,ct_scan_findings,distant_metastasis_anatomic_site,metastatic_neoplasm_initial_diagnosis_anatomic_site,
                  `metastatic_neoplasm_initial_diagnosis_anatomic_site[1]`,`metastatic_neoplasm_initial_diagnosis_anatomic_site[2]`,
                  `metastatic_neoplasm_initial_diagnosis_anatomic_site[3]`,new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_occurrence_anatomic_site_text,
                  number_of_lymphnodes_positive_by_he,other_malignancy_anatomic_site)
    #write.csv(dat1,str_glue("met_anno/{proj}_met_anno.txt"))
    met_anno <- data.table::fread("met_anno/TCGA-ACC_met_anno.txt", na.strings=c("","NA"))
    dat <- left_join(dat, met_anno, by="fileID")
    write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Clinical_annotation_metastatic_locations_{proj}.csv"))
  }
BLCA <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), header = TRUE,na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(fileID,metastatic_site,
                  `metastatic_site[1]`,`metastatic_site[2]`,`metastatic_site[3]`,
                  new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_occurrence_anatomic_site_text, 
                  other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,other_metastatic_site) 
  write.csv(dat1,str_glue("met_anno/{proj}_met_anno.txt"))
  met_anno <- data.table::fread(str_glue("met_anno/{proj}_met_anno.txt")) %>%
    dplyr::select(-index)
  
  dat <- left_join(dat,met_anno, by ="fileID")
  write.csv(df, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Clinical_annotation_metastatic_locations_{proj}.csv"))
  
}
BRCA <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"),header = TRUE,na.strings=c("","NA"))
  dat1 <- dat %>%
    dplyr::select(
      fileID,bcr_patient_barcode,metastatic_site_at_diagnosis,metastatic_site_at_diagnosis_other,
      `metastatic_site_at_diagnosis[1]`,
      `metastatic_site_at_diagnosis[2]`,
      `metastatic_site_at_diagnosis[3]`,
      `metastatic_site_at_diagnosis[4]`,
      new_neoplasm_event_occurrence_anatomic_site,
      new_tumor_event_after_initial_treatment,
      other_malignancy_anatomic_site
    ) %>%
    unite("Loci", metastatic_site_at_diagnosis:other_malignancy_anatomic_site, sep= ",", 
          remove = FALSE)
  met_anno <- dat1 %>% dplyr::select(fileID,Loci)
  #write.csv(met_anno,str_glue("met_anno/{proj}_met_anno.txt"))
  met_anno <- data.table::fread("met_anno/TCGA-BRCA_met_anno.txt") %>%
    dplyr::select(-index)
  
  df <- left_join(dat,met_anno, by ="fileID")
  
}
CESC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), header = TRUE,na.strings=c("","NA"))
  dat1 <-  dat %>% dplyr::select(
    fileID,
    bcr_patient_barcode,
    diagnostic_ct_result_outcome,
    `diagnostic_ct_result_outcome[1]`,
    `diagnostic_ct_result_outcome[2]`,
    `diagnostic_ct_result_outcome[3]`,
    `diagnostic_ct_result_outcome[4]`,
    `diagnostic_ct_result_outcome[5]`,
    `diagnostic_ct_result_outcome[6]`,
    diagnostic_mri_result_outcome,
    `diagnostic_mri_result_outcome[1]`,
    `diagnostic_mri_result_outcome[2]`,
    `diagnostic_mri_result_outcome[3]`,
    `diagnostic_mri_result_outcome[4]`,
    `diagnostic_mri_result_outcome[5]`,
    `diagnostic_mri_result_outcome[6]`,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    new_neoplasm_event_type,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text
  )
  #write.csv(dat1,str_glue("met_anno/{proj}_met_anno.txt"))
  met_anno <- data.table::fread(str_glue("met_anno/{proj}_met_anno.txt", na.strings=c("","NA"))) %>%
    dplyr::select(fileID, Metastasis)
  write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Clinical_annotation_metastatic_locations_{proj}.csv"))
  

}
CHOL <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), header = TRUE, na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text
  )
  write.csv(dat1,str_glue("met_anno/{proj}_met_anno.txt"))
  df <- data.table::fread(str_glue("met_anno/{proj}_met_anno.txt", na.strings=c("","NA"))) %>%
    dplyr::select(fileID, Metastasis)
  dat <- left_join(dat,df, by ="fileID")
  write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Clinical_annotation_metastatic_locations_{proj}.csv"))
  
}
COAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(
    fileID,
    new_neoplasm_event_type,
    non_nodal_tumor_deposits,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text
  )
  write.csv(dat1,str_glue("met_anno/{proj}_met_anno.txt"))
  df <- data.table::fread(str_glue("met_anno/{proj}_met_anno.txt", na.strings=c("","NA"))) %>%
    dplyr::select(fileID, Metastasis)
  dat <- left_join(dat,df, by ="fileID")
  write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Clinical_annotation_metastatic_locations_{proj}.csv"))
  
}
DLBC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(
    fileID,
    bone_marrow_involvement,
    lymph_node_involvement_site,
    `lymph_node_involvement_site[1]`,
    `lymph_node_involvement_site[2]`,
    `lymph_node_involvement_site[3]`,
    `lymph_node_involvement_site[4]`,
    `lymph_node_involvement_site[5]`,
    `lymph_node_involvement_site[6]`,
    `lymph_node_involvement_site[7]`,
    `lymph_node_involvement_site[8]`,
    `lymph_node_involvement_site[9]`,
    `lymph_node_involvement_site[10]`,
    tumor_tissue_site,
    `tumor_tissue_site[1]`,
    `tumor_tissue_site[2]`,
    `tumor_tissue_site[3]`,
    `tumor_tissue_site[5]`,
    `tumor_tissue_site[6]`,
    `tumor_tissue_site[6]`
  )
  #write.csv(dat1,str_glue("met_anno/{proj}_met_anno.txt"))
  met_anno <- data.table::fread(str_glue("met_anno/{proj}_met_anno.txt", na.strings=c("","NA"))) %>%
    dplyr::select(fileID, Metastasis)
  dat <- left_join(dat,met_anno, by ="fileID")
  write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Clinical_annotation_metastatic_locations_{proj}.csv"))
  
}
ESCA <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
GBM <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
HNSC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
KICH <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
KIRC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
KIRP <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
LAML <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
LGG <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
LIHC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
LUAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
LUSC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
MESO <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
OV <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
PAAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
PCPG <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
PRAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
READ <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
SARC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
SKCM <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
STAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
TGCT <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
THCA <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
THYM <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
UCEC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
UCS <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
UVM <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  
}
