#CUP analysis
setwd("/mnt/storage/mskaro1/CUP/GDCdata")

library(TCGAbiolinks)

# Download and prepare 450k methylation beta values for analysis

#projects <- TCGAbiolinks:::getGDCprojects()$project_id
#projects <- projects[grepl('^TCGA',projects,perl=T)]
projects <- c("TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD","TCGA-READ",
              "TCGA-DLBC","TCGA-ESCA","TCGA-GBM","TCGA-HNSC","TCGA-KICH","TCGA-KIRC",
              "TCGA-KIRP","TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC",
              "TCGA-MESO","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-SARC","TCGA-SKCM",
              "TCGA-STAD","TCGA-TGCT","TCGA-THCA","TCGA-THYM","TCGA-UCEC","TCGA-UVM")
match.file.cases.all <- NULL
for(proj in projects){
  print(proj)
  query <- GDCquery(project = proj,
                    data.category = "DNA Methylation",
                    platform = "Illumina Human Methylation 450")
  match.file.cases <- getResults(query,cols=c("cases","file_name"))
  match.file.cases$project <- proj
  match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
  tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
           error = function(e) GDCdownload(query, method = "client"))
  betas <- GDCprepare(query)
  save(betas, file = str_glue("{proj}_beta.RData"))
}
# This will create a map between idat file name, cases (barcode) and project
readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")
# code to move all files to local folder
for(file in dir(".",pattern = ".idat", recursive = T)){
  TCGAbiolinks::move(file,basename(file))
}
