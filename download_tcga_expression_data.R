#' This script will download all of the TCGA data that we have used in this project. This data is obtained from the firebrowse.org repository.


#' To run from the Bionimbus PDC, first run the following commands from the command line to set up the HTTP proxy:
#' export http_proxy=http://cloud-proxy:3128
#' export https_proxy=http://cloud-proxy:3128

#' Set the root directory where the data will be stored. NB: this directory needs to be set / created based on your own system!! The file "theRootDir.R" is sourced by most of these scripts and should include the desired root directory used by your entire analysis.
source("theRootDir.R")

#' create the "dataIn/" directory if it doesn't already exist.
dir.create(paste(theRootDir, "dataIn/", sep=""), showWarnings = FALSE)


#' The data is organized on firebrowse by cancer type. They use these disease abbreviations to access the different folders containing the various data.
diseaseAbbrvs <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "FPPP", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "STES", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")


#' Download all the TCGA RNA-seq data.
missingAbrvsRnaSeq <- c(10, 31) # there is no RNA-seq data for "FPPP" or "STAD"
rnaSeqDiseaseAbbrvs <- diseaseAbbrvs[-missingAbrvsRnaSeq]
rnaSeqFilesDir <- paste(theRootDir, "dataIn/rnaSeq/", sep="")
dir.create(rnaSeqFilesDir, showWarnings = FALSE) # make this directory if it doesn't exist.
for(i in 1:length(rnaSeqDiseaseAbbrvs))
{
  fname <- paste("gdac.broadinstitute.org_", rnaSeqDiseaseAbbrvs[i], ".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2015082100.0.0.tar.gz", sep="")
  download.file(paste("http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/", rnaSeqDiseaseAbbrvs[i], "/20150821/gdac.broadinstitute.org_", rnaSeqDiseaseAbbrvs[i], ".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2015082100.0.0.tar.gz", sep=""), paste(rnaSeqFilesDir, fname, sep=""))
}

# Unzip the downloaded ".tar.gz" RNA-seq data! NB, this command has been tested in Linux. It may not work in Windows. If it does not work, please extract these files manually using software such as 7zip.
thegzFiles <-  paste(rnaSeqFilesDir, dir(rnaSeqFilesDir), sep="")
sapply(thegzFiles, untar, exdir=rnaSeqFilesDir)



print(sessionInfo())