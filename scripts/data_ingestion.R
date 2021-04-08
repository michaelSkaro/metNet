# Data ingestion in R
#   1. Install and load all of the necessary packages for the analysis
#   2. Specify the location of the data in the directory to the idat files
#   3. Process all of the raw data for the project and select only the IDs that
#      exist in the selected training data      
#   4. Feed these data as a csv into the ML model
#   5. Biological analysis? up to andrea and team members

# Methylation ingestion
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylationArrayAnalysis")
# that takes a long time to install


# load packages required for analysis if they aren't 
# included in the methylation array analysis load

#Biocmanager::install("knitr")
#Biocmanager::install("limma")
#Biocmanager::install("minfi")
#Biocmanager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#Biocmanager::install("IlluminaHumanMethylation450kmanifest")
#Biocmanager::install("RColorBrewer")
#Biocmanager::install("missMethyl")
#Biocmanager::install("minfiData")
#Biocmanager::install("Gviz")
#Biocmanager::install("DMRcate")
#install.packages("stringr")

library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(methylationArrayAnalysis)
# This is where data will be input from the front end
setwd("/mnt/storage/mskaro1/methylation_illumina450k/data/")

# read in the annotation that matches the illu ina 450k
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# list the data in the working directory

data_files <- list.files(recursive = TRUE)
# Have the user input some kind of sample annotation that just has to be
# in a datatable with one instance per row and any number of descriptive columns they want

# will make later for reference but this is what will be passed from user
targets <- read.metharray.sheet("/mnt/storage/mskaro1/methylation_illumina450k/data/", pattern = "csv$", ignore.case = TRUE,
                                recursive = TRUE, verbose = TRUE)
                                
# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)


# calcuate the pvalues for the row cgIDs in the data set. We may need to reduce the 
# sample set and get rid of some of the outlier samples
detP <- detectionP(rgSet)
# is this necessary?

mSetSq <- preprocessIllumina(rgSet, normalize = "controls") 
# calculate the normalized beat values

bVals <- as.data.frame(getBeta(mSetSq))

write.csv(bVals, "/mnt/storage/mskaro1/methylation_illumina450k/data/Normalized_betas.csv")
