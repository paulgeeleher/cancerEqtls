#' This script runs the eQTL analysis in R, using a conventional model and an interaction model, contolling for copy number variation.

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

#' Load the filtered subsetted files needed for the GWAS analysis
load(file=theRootDir %&% "Results/rDatas/gwasFilteredInputData_WithQNexprAndCpe.RData") # tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin, geneLocation_filt_qnDat, theCpeProps_matched, covariatesVec_filt, expressionMatNorm_filt, expressionMat_filt, genotypesMat_filt, genotypePcs, snpsLocation_filt, geneLocation_filt, 

#' Load the methylation data, which was obtained from FireBrowse (http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BRCA/20160128/gdac.broadinstitute.org_BRCA.Methylation_Preprocess.Level_3.2016012800.0.0.tar.gz)
read.delim("/mnt/data_scratch/prediXcanProj/dataIn/methylation/gdac.broadinstitute.org_BRCA.Methylation_Preprocess.Level_3.2016012800.0.0/BRCA.meth.by_mean.data.txt")
methyData <- read.delim(theRootDir %&% "dataIn/methylation/gdac.broadinstitute.org_BRCA.Methylation_Preprocess.Level_3.2016012800.0.0/BRCA.meth.by_mean.data.txt", as.is=T)
rownames(methyData) <- methyData[,1]
methyDataMat <- data.matrix(methyData[-1, -1])

is01 <- which(substring(colnames(methyDataMat), 14, 15) == "01")
methyDataMat01 <- methyDataMat[, is01]
sum(duplicated(substring(colnames(methyDataMat01), 9, 12))) # there are no duplicated sample IDs once we remove the non 01 (i.e. non primary site) samples
colnames(methyDataMat01) <- substring(colnames(methyDataMat01), 9, 12)

#' Only keep genes that have a Methy estimate in at least 200 samples (some have many missing values)
dim(methyDataMat01)
dim(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin)
commonGNames <- intersect(rownames(methyDataMat01), rownames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin))
numNas <- apply(methyDataMat01[commonGNames, ], 1, function(x)return(sum(is.na(x))))
names(numNas) <- rownames(methyDataMat01[commonGNames, ])
commonGeneNames <- names(numNas)[numNas < 200]
commonSamples <- intersect(colnames(methyDataMat01), colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin))
length(commonGeneNames)
length(commonSamples)

# Filter the datasets
methyMatComm <- methyDataMat01[commonGeneNames, commonSamples]
exprMatComm <- tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[commonGeneNames, commonSamples]
genotypeMatComm <- genotypesMat_filt[, commonSamples]
geneLocation_filt_qnDat_Comm <- geneLocation_filt_qnDat[commonGeneNames, ] 


#' Create a list mapping SNPs to genes, with a 1 megabase window around the gene.
snpsLocation_filt[,4] <- snpsLocation_filt[,3] - 500000 # lets do a megabase window around the gene.
snpsLocation_filt[,4][snpsLocation_filt[,4] < 0] <- 0
snpsLocation_filt[,5] <- snpsLocation_filt[,3] + 500000
colnames(snpsLocation_filt)[c(4,5)] <- c("minusMB", "plusMB")
rownames(snpsLocation_filt) <- snpsLocation_filt[,1]

#' We will match snps to genes using the genomic ranges package in R, which will allow us to do this very quickly.
library(GenomicRanges) # this library should be able to do this fairly efficiently I think.
geneRanges <- GRanges(seqnames=Rle(geneLocation_filt_qnDat_Comm[,2]), ranges=IRanges(geneLocation_filt_qnDat_Comm[,3], geneLocation_filt_qnDat_Comm[,4]))
names(geneRanges) <- geneLocation_filt_qnDat_Comm[,1]
snpLocations <- GRanges(seqnames=Rle(snpsLocation_filt[,2]), ranges=IRanges(snpsLocation_filt[,4], snpsLocation_filt[,5]))
names(snpLocations) <- snpsLocation_filt[,1]

countOut <- findOverlaps(geneRanges,snpLocations, type="any", ignore.strand=T)
countOut_mat <- as.matrix(countOut)
snpSplit <- split(countOut_mat[,1], countOut_mat[,2])

snpsToGenesMat <- countOut_mat
snpsToGenesMat[,1] <- geneLocation_filt_qnDat_Comm[,1][countOut_mat[,1]]
snpsToGenesMat[,2] <- snpsLocation_filt[,1][countOut_mat[,2]]
snpsToGenesList <- split(snpsToGenesMat[,1], snpsToGenesMat[,2])
length(snpsToGenesList)
snpsToGenesList_methyl <- snpsToGenesList

#' Load PEER factors from the normalized gene expression data (tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin). (create in run_cis_eqtl_analysis_in_r.R)
load(file=theRootDir %&% "Results/rDatas/thePeers.RData") # theFactors_NoInteractionAnalysis, theFactors_correctedAnaysis, theFactors_correctedAnaysis_filtered, 
rownames(theFactors_correctedAnaysis) <- colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin)
rownames(theFactors_NoInteractionAnalysis) <- colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin)
theFactors_correctedAnaysis_Comm <- theFactors_correctedAnaysis[commonSamples, ]
theFactors_NoInteractionAnalysis_Comm <- theFactors_NoInteractionAnalysis[commonSamples, ]


#' 2 analyses here:
#' Conventional with Methy covariate.
# #' Interaction with Methy covariate.
coefsList_conv_methys <- list(length(snpsToGenesList))
# coefsList_int_methys <- list(length(snpsToGenesList))
coefsList_intInv_methys <- list(length(snpsToGenesList))

## YOU need a different covariate for each gene (each response variable), I dont think that can be coded in the formula below.
## Doesnt really make sense that this would be possible, i.e. youd need multiple matix inveresions, 
## so one might as well code these separately (or do it in C!!? (probably a modest performance gain, maybe not worht it))
# This is gonna have to be a nested loop. This will be slow.
thisInvProp <- (1-theCpeProps_matched[commonSamples])
for(i in 1:length(snpsToGenesList)) # [1] 615568
{
  coefsList_conv_methys[[i]] <- list()
#   coefsList_int_methys[[i]] <- list()
  coefsList_intInv_methys[[i]] <- list()
  for(j in 1:length(snpsToGenesList[[i]]))
  {
      coefsList_conv_methys[[i]][[j]] <- coef(summary(lm(exprMatComm[snpsToGenesList[[i]][j],]~genotypeMatComm[names(snpsToGenesList)[i], ] + genotypePcs$x[commonSamples,1:3] + theFactors_NoInteractionAnalysis_Comm + methyMatComm[snpsToGenesList[[i]][j], ])))
#       coefsList_int_methys[[i]][[j]] <- coef(summary(lm(exprMatComm[snpsToGenesList[[i]][j],]~(genotypeMatComm[names(snpsToGenesList)[i], ]*theCpeProps_matched[commonSamples]) + genotypePcs$x[commonSamples,1:3] + theFactors_correctedAnaysis_Comm + methyMatComm[snpsToGenesList[[i]][j], ]))) # This will get estimates for the "normal" component.
      coefsList_intInv_methys[[i]][[j]] <- coef(summary(lm(exprMatComm[snpsToGenesList[[i]][j],]~(genotypeMatComm[names(snpsToGenesList)[i], ]*thisInvProp) + genotypePcs$x[commonSamples,1:3] + theFactors_correctedAnaysis_Comm + methyMatComm[snpsToGenesList[[i]][j], ]))) # by inverting the proportions like so, we will calculate an effect for the "cancer" component, not the normal component.
  }
  names(coefsList_conv_methys[[i]]) <- snpsToGenesList[[i]]
#   names(coefsList_int_methys[[i]]) <- snpsToGenesList[[i]]
  names(coefsList_intInv_methys[[i]]) <- snpsToGenesList[[i]]
  print(i)
}
names(coefsList_conv_methys) <- names(snpsToGenesList)
# names(coefsList_int_methys) <- names(snpsToGenesList)
names(coefsList_intInv_methys) <- names(snpsToGenesList)

# dir.create(theRootDir %&% "data/forMatrixEQTL/output/", recursive=TRUE, showWarnings = FALSE)
# Save the data, also, snpsToGenesList_methyl, which I need for processing this output.
save(snpsToGenesList_methyl, coefsList_conv_methys, coefsList_intInv_methys, file=theRootDir %&% "data/forMatrixEQTL/output/coefLists_methys.RData")
# save(coefsList_int_methys, file=theRootDir %&% "data/forMatrixEQTL/output/coefsList_int_methys.RData")

 

print(sessionInfo())