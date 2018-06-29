#' # This code uses RUV based methods (described in PMID:25150836) to batch correct the TCGA gene expression data.

#' Set the root directory where the data will be stored. NB: this directory needs to be set / created based on your own system!! The file "theRootDir.R" is sourced by most of these scripts and should include the desired root directory used by your entire analysis. Default value of "theRootDir" variable is "/mnt/data_scratch/prediXcanProj/"
source("theRootDir.R")

#' Create the directory to store the figures.
dir.create(paste(theRootDir, "figures/", sep = ""), showWarnings = FALSE)

#' Load the pRRophetic library
library("pRRophetic")

#' A function for getting p-value from linear regression fit in R. Credit source: http://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

#' Read the directories in which the files are contained.
#' These data were downloaded for firebrowse.org.
#' We provide a script to automate the downloading of these data (see download_tcga_data.R).
theRnaSeqDir <- paste(theRootDir, "dataIn/rnaSeq/", sep="") # the directory containing the RNA-seq data.
theDirs <- dir(theRnaSeqDir)
theDirs <- theDirs[-grep(".tar.gz", theDirs, fixed=T)] # ignore the .tar.gz files.

#' Some of the TCGA data are redundant, i.e. the same samples are contained in different datasets, we need to remove these duplicated samples.
cancerTypeNames <- sapply(sapply(strsplit(theDirs, ".", fixed=T), function(a)return(strsplit(a[[3]], "_"))), function(b)return(b[2])) #
removeTypes <- c("COADREAD", "GBMLGG", "KIPAN", "STES") ## NB these must not be inlcuded as they are totally redundant, i.e. these samples are identical to those contained in other folders.
theDirsFilt <- theDirs[!cancerTypeNames %in% removeTypes]
cancerTypeNames <- cancerTypeNames[!cancerTypeNames %in% removeTypes]

#' Load ALL of the data. N.B. This requires a very large amound of memory. It has been tested on a machine with 128Gb or RAM.
#' Note, this code assumes you have obtained the same RNA-seq data as downloaded by the "download_tcga_data.R" script
tpmMatList <- list()
for(i in 1:length(theDirsFilt))
{
  theFile <- dir(paste(theRnaSeqDir, theDirsFilt[i], sep=""))[grep("MANIFEST", dir(paste(theRnaSeqDir, theDirsFilt[i], sep="")), invert=T)]
  
  tpmDatMat <- read.delim(paste(theRnaSeqDir, theDirsFilt[i], "/", theFile, sep=""), as.is=T)

  tpmDatMat_tpm <- apply(tpmDatMat[-1,which(tpmDatMat[1,] == "scaled_estimate")], 2, as.numeric)
  tpmDatMat_tpm <- tpmDatMat[-1,which(tpmDatMat[1,] == "scaled_estimate")]
  tpmDatMat_tpm <- apply(tpmDatMat_tpm, 2, as.numeric)

  geneNames <- do.call(cbind, strsplit(tpmDatMat[, "Hybridization.REF"], "|", fixed=TRUE))[1,][-1]
  rownames(tpmDatMat_tpm) <- geneNames
  colnames(tpmDatMat_tpm) <- substr(colnames(tpmDatMat_tpm), 1, 28)

  tpmDatMat_tpm_logged <- log((tpmDatMat_tpm*1000000)+1) # transform the data

  tpmMatList[[i]] <- tpmDatMat_tpm_logged
  
}
rnames <- lapply(tpmMatList, rownames)

#' Get the cancer types.
names(tpmMatList) <- cancerTypeNames
numSampls <- sapply(tpmMatList, ncol)
cancerTypesVec <- character()
for(i in 1:length(cancerTypeNames)){cancerTypesVec <- c(cancerTypesVec, rep(cancerTypeNames[i], numSampls[i]))}

allExprData <- do.call(cbind, tpmMatList)

#' Save the gene expression matrix here, we will use this again.
save(allExprData, file=paste(theRootDir, "dataIn/allExprData.RData", sep="")) # allExprData


#' Now that the data are loaded, we wish to calculate the prinicple components that will be used to "remove unwanted variation" (RUV).
#' First create a matrix of the expression data that is standardized by cancer type, as we do not wish to remove this variabilty.
allCancerTypes <- unique(cancerTypesVec)
standardizeByCancerType <- allExprData
for(i in 1:length(allCancerTypes))
{
  for(j in 1:nrow(allExprData))
  {
    vec <- standardizeByCancerType[j, cancerTypesVec %in% allCancerTypes[i]]
    standardizeByCancerType[j, cancerTypesVec %in% allCancerTypes[i]] <- ((vec-mean(vec))/sd(vec))    
  }
}
save(standardizeByCancerType, file=paste(theRootDir, "dataIn/standardizeByCancerType.RData", sep=""))

#' Find a set of 250 genes that are expressed in all samples and exhibit the lowest variabilty.
zeroExprSums <- apply(allExprData, 1, function(r)sum(r == 0)) # get the number of samples in which each gene isn't expressed.
consistentlyExpressedGenes <- which(zeroExprSums == 0) # the genes that are expressed in every sample.
varsExprssed <- apply(allExprData[consistentlyExpressedGenes, ], 1, var)
veryLowVarExpressed <- names(sort(varsExprssed)[1:250]) # expressed genes with the lowest variabiltiy
medExprGene <- apply(allExprData, 1, median)

#' Calculate the principal componets of these genes.
rowNoNas <- which(apply(standardizeByCancerType, 1, function(row)return(sum(is.na(row)))) == 0)
noNasLowVar <- intersect(rownames(standardizeByCancerType[rowNoNas, ]), veryLowVarExpressed)
ruvPcs_standardized <- prcomp(t(standardizeByCancerType[noNasLowVar[-1], ]))

#' As a sanity check, are these RUV components actually correlated with batch ID, if they are, that is a very good thing, because it shows that these principal components have identified the real batches in a completely unbiased way.
batchIds <- sapply(strsplit(colnames(allExprData), ".", fixed=T), function(l)l[6])
summary(lm(ruvPcs_standardized$x[,1]~factor(batchIds)+factor(cancerTypesVec)))
theRsquareds <- numeric()
thePvals <- numeric()
for(i in 1:100)
{
  theMod <- lm(ruvPcs_standardized$x[,i]~factor(batchIds))
  theRsquareds[i] <- summary(theMod)$r.squared
  thePvals[i] <- lmp(theMod)
}
print(theRsquareds[1:10])
print(thePvals[1:10])

#' Supplementary Figure: the proportion of variabilty in "batch" caputred by each of the RUV principal components. This clearly levels off after 10, suggesting that 10 PCs, the number proposed by the original authors, is appropriate in this case.
pdf(paste(theRootDir, "figures/theRsquareds_ruv_components_against_batch.pdf", sep=""), width=4, height=4)
plot(theRsquareds, pch=20, ylab="R squared", xlab="Principal Component")
abline(v=10, col="red")
dev.off()

#' finally, for each gene, regress out the 10 RUV PCs that we have calculated on the matrix of standardized expression data, then save this updated matrx.
tenRuvNewStandardApproach <- allExprData
for(i in 1:nrow(allExprData))
{
  tenRuvNewStandardApproach[i,] <- residuals(glm(allExprData[i,]~ruvPcs_standardized$x[, 1:10], family="quasipoisson"))
}
save(tenRuvNewStandardApproach, cancerTypesVec, file=paste(theRootDir, "dataIn/tenRuvNewStandardApproach.RData", sep=""))



print(sessionInfo())