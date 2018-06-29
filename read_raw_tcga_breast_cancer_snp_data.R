# Draft file was "preProcessRawSnpFiles.R".

#' This file will take the raw SNP data, which is spread across disparate files in a non-usable format and it will load it into a matrix of rsID x TCGA id, which we can use for downstream analysis.
#' Note, this TCGA SNP data was obtained from GDC, and permission must be obtained to download these protected data (see below).
#' The data can then be downloaded using the "gdc-client" tool (provided on the GDC website). A "token" file is required to download the protected data. A "manifest" files is required to download the correct data. The "manifest" file used to download these data is provided with these instructions (see index.html file)
#' There are 2 files required here to annotate the SNP array data, i.e. map the Affymetrix probe IDs to RS IDs. These files can be obtained as part of the metadata provided with birdsuite, i.e. the file birdsuite_metadata_1.5.5.tgz which can be downloaded from https://www.broadinstitute.org/birdsuite/birdsuite-downloads. We have also provided the two relevant files (GenomeWideSNP_6.rs_snp_map & GenomeWideSNP_6_alleles.csv).
#' This script assumes the folders containing the SNP files (downloaded from GDC using gdc-client) have been stored in a directory /mnt/data_scratch/prediXcanProj/data/genotypes/rawGenotypes/


#' NB: To reproduce this analysis, make sure you use a directory structure similar to what is used in this file, or change the directory structure.

#' Set the root directory where the data will be stored. NB: this directory needs to be set / created based on your own system!! The file "theRootDir.R" is sourced by most of these scripts and should include the desired root directory used by your entire analysis. Default value of "theRootDir" variable is "/mnt/data_scratch/prediXcanProj/"
source("theRootDir.R")


#' Create the data/ and genotype/ directories if they don't already exist.
#' Create the directory to store the figures.
dir.create(paste(theRootDir, "data/genotypes/", sep = ""), showWarnings = FALSE, recursive=TRUE)

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

#' Load the metadata for the Affy SNP 6.0 array. 
theMap <- read.delim(theRootDir %&% "GenomeWideSNP_6.rs_snp_map", as.is=T)
rsIds <- theMap[,1]
names(rsIds) <- theMap[,2]

#' This file contains the actual genotypes (G, C, T, As) for each RSid
theAlleles <- read.delim(theRootDir %&% "GenomeWideSNP_6_alleles.csv", as.is=T)

#' Write out a dataframe that has affy snp ids
alleleData <- cbind(theAlleles, rsIds[theAlleles[,1]])
alleleData_fin <- alleleData[-which(is.na(alleleData[,4])), c(4,2,3)]
save(alleleData_fin, file=theRootDir %&% "data/genotypes/alleleData_fin.RData")


#' Load the mappings of TCGA file names to TCGA IDs. There were retrieved using the following query (with assistance from GDC help desk): https://gdc-api.nci.nih.gov/legacy/files?format=tsv&size=10000&fields=file_name,file_id,cases.samples.portions.analytes.aliquots.submitter_id&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.platform%22%2C%22value%22%3A%5B%22Affymetrix+SNP+Array+6.0%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.primary_site%22%2C%22value%22%3A%5B%22Breast%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Genotypes%22%5D%7D%7D%5D%7D
fNameToTcga <- read.delim(theRootDir %&% "tcgaFileNameToTcgaIdMapping.txt", as.is=T)
tcgaIds <- fNameToTcga[,3]
names(tcgaIds) <- fNameToTcga[,2]

#' Iterate through the individual files and load into a big matrix. I have downloaded genotype data called from both tumor and normal samples.
setwd(theRootDir %&% "data/genotypes/rawGenotypes/")
folderNames <- dir()
genotypesList <- list()
fileNamesList <- list()
callsList <- list()
confidenceList <- list()
for(i in 1:length(folderNames))
{
  filesInThisFolder <- dir(folderNames[i])
  dataFileName <- filesInThisFolder[grep(".data.txt", filesInThisFolder)]
  fileNamesList[[i]] <- dataFileName
  theData <- read.delim(paste(folderNames[i], "/", dataFileName, sep=""), as.is=T)[-1,]
  callsList[[i]] <- as.numeric(theData[,2])
  names(callsList[[i]]) <- theData[,1]
  confidenceList[[i]] <- as.numeric(theData[,3])
  names(confidenceList[[i]]) <- theData[,1]
  print(i)
}

#' Convert these lists to matrices...
callsMatrix <- do.call(cbind, callsList)
colnames(callsMatrix) <- unlist(fileNamesList)
confidenceMatrix <- do.call(cbind, confidenceList)
colnames(confidenceMatrix) <- unlist(fileNamesList)

#' Map the rsIds onto the rownames in the matrix...
rsIdsOrd <- rsIds[rownames(callsMatrix)]
rownames(callsMatrix) <- rsIds[rownames(callsMatrix)]
rownames(confidenceMatrix) <- rsIds[rownames(confidenceMatrix)]

#' Map the "file names" that are currently the column names to TCGA IDs.
colnames(callsMatrix) <- tcgaIds[colnames(callsMatrix)]
colnames(confidenceMatrix) <- tcgaIds[colnames(confidenceMatrix)]


#' For some reason there's a duplicated samples, but the genotypes don't match. Hmmm. I will need to make sure everything is okay here.
which(tcgaIds[colnames(callsMatrix)] == "TCGA-EW-A424-10A-01D-A242-01")
# CROZE_p_TCGA_254_256_255_N_GenomeWideSNP_6_A03_1301360.birdseed.data.txt 
#                                                                      309 
# CROZE_p_TCGA_254_256_255_N_GenomeWideSNP_6_A02_1301374.birdseed.data.txt 
#                                                                     1693 
sum(callsMatrix[,309] != callsMatrix[,1693])
# [1] 424042


#' Save matrices separately beacuse these are huge.
save(callsMatrix, file=theRootDir %&% "data/genotypes/geneotypeMatrices/callsMatrix.RData")
save(confidenceMatrix, file=theRootDir %&% "data/genotypes/geneotypeMatrices/confidenceMatrix.RData")

print(sessionInfo())