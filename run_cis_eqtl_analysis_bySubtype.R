#' This script will run conventaional and interaction models, stratified by PAM50 defined subtypes. We have may have enough samples for  Basal-like, HER2-enriched, Luminal A, Luminal B.

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

#'  Load the PAM50 Data obtained from the Supplementary Information of this paper: https://www.nature.com/articles/nature11412.
pam50Subtypes <- read.csv(theRootDir %&% "tcgaBrcaPam50SubtypesFrom_pmid23000897.csv", as.is=T, header=T, skip=1)
table(pam50Subtypes[, "PAM50.mRNA"]) # The numbers of samples for the various PAM50 subtypes.
pam50Info <- pam50Subtypes[, "PAM50.mRNA"]
names(pam50Info) <- pam50Subtypes[, "Complete.TCGA.ID"]


#' There's also PAM50 information provided here: https://breast-cancer-research.biomedcentral.com/articles/10.1186/s13058-016-0724-2 (additional file 2, sheet 2 of this excel file). I found this because of this thread on the bioconductor support form: https://support.bioconductor.org/p/107669/
#' Except in this article, there is information for many more samples.....
pam50SubtypesFull <- read.csv(theRootDir %&% "tcgaPam50_fromNetanelyEtAl.csv", as.is=T, header=F, skip=1)
pam50SubtypesFull01 <- pam50SubtypesFull[substring(pam50SubtypesFull[,1], 14, 15) == "01", ]
table(pam50SubtypesFull01[, 7])
#  Basal   Her2   LumA   LumB Normal 
#    183     78    534    203     37
pam50SubtypesFullNorm <- pam50SubtypesFull01[, 7]
names(pam50SubtypesFullNorm) <- substring(pam50SubtypesFull01[,1], 1, 12)

# Use the same string as provided by the TCGA so I can compare these and see if the PAM50 calls are concordant
# TCGA and Netanely et al PAM 50 calls are mostly (88%) concordant.
pam50SubtypesFullNorm[pam50SubtypesFullNorm == "LumA"] <- "Luminal A"
pam50SubtypesFullNorm[pam50SubtypesFullNorm == "LumB"] <- "Luminal B"
pam50SubtypesFullNorm[pam50SubtypesFullNorm == "Basal"] <- "Basal-like"
pam50SubtypesFullNorm[pam50SubtypesFullNorm == "Her2"] <- "HER2-enriched"
commonSamples <- names(pam50SubtypesFullNorm)[names(pam50SubtypesFullNorm) %in% names(na.omit(pam50Info))]
sum(pam50SubtypesFullNorm[commonSamples] == pam50Info[commonSamples]) # [1] 453
sum(!pam50SubtypesFullNorm[commonSamples] == pam50Info[commonSamples]) # [1] 63



#' Load the filtered subsetted files needed for the GWAS analysis
load(file=theRootDir %&% "Results/rDatas/gwasFilteredInputData_WithQNexprAndCpe.RData") # tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin, geneLocation_filt_qnDat, theCpeProps_matched, covariatesVec_filt, expressionMatNorm_filt, expressionMat_filt, genotypesMat_filt, genotypePcs, snpsLocation_filt, geneLocation_filt, 

#' Create a list mapping SNPs to genes, with a 1 megabase window around the gene.
snpsLocation_filt[,4] <- snpsLocation_filt[,3] - 500000 # lets do a megabase window around the gene.
snpsLocation_filt[,4][snpsLocation_filt[,4] < 0] <- 0
snpsLocation_filt[,5] <- snpsLocation_filt[,3] + 500000
colnames(snpsLocation_filt)[c(4,5)] <- c("minusMB", "plusMB")
rownames(snpsLocation_filt) <- snpsLocation_filt[,1]

#' We will match snps to genes using the genomic ranges package in R, which will allow us to do this very quickly.
library(GenomicRanges) # this library should be able to do this fairly efficiently I think.
geneRanges <- GRanges(seqnames=Rle(geneLocation_filt_qnDat[,2]), ranges=IRanges(geneLocation_filt_qnDat[,3], geneLocation_filt_qnDat[,4]))
names(geneRanges) <- geneLocation_filt_qnDat[,1]
snpLocations <- GRanges(seqnames=Rle(snpsLocation_filt[,2]), ranges=IRanges(snpsLocation_filt[,4], snpsLocation_filt[,5]))
names(snpLocations) <- snpsLocation_filt[,1]

countOut <- findOverlaps(geneRanges,snpLocations, type="any", ignore.strand=T)
countOut_mat <- as.matrix(countOut)
snpSplit <- split(countOut_mat[,1], countOut_mat[,2])

snpsToGenesMat <- countOut_mat
snpsToGenesMat[,1] <- geneLocation_filt_qnDat[,1][countOut_mat[,1]]
snpsToGenesMat[,2] <- snpsLocation_filt[,1][countOut_mat[,2]]
snpsToGenesList <- split(snpsToGenesMat[,1], snpsToGenesMat[,2])

length(snpsToGenesList)
# save(snpsToGenesList, file="/mnt/data_scratch/prediXcanProj/Results/rDatas/snpsToGenesList.RData")


#' Lets create expression and genotype matrices for the 4 various subtypes, A, B, basal, HER2.
aSamples <- intersect(substring(names(pam50SubtypesFullNorm)[which(pam50SubtypesFullNorm == "Luminal A")], 9, 12), colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin)) # samples for which there's both expression and PAM50 annotations
bSamples <- intersect(substring(names(pam50SubtypesFullNorm)[which(pam50SubtypesFullNorm == "Luminal B")], 9, 12), colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin))
basalSamples <- intersect(substring(names(pam50SubtypesFullNorm)[which(pam50SubtypesFullNorm == "Basal-like")], 9, 12), colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin))
her2Samples <- intersect(substring(names(pam50SubtypesFullNorm)[which(pam50SubtypesFullNorm == "HER2-enriched")], 9, 12), colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin))

# the numbers of samples left in each subtype (for the new PAM50 calls)
length(aSamples) # [1] 451
length(bSamples) # [1] 175
length(basalSamples) # [1] 159
length(her2Samples) # [1] 66

varsAsText <- c("aSamples", "bSamples", "basalSamples", "her2Samples")

#' Create PEER factors for each of the 4 subtypes. NOTE: Because of the smaller sample sizes I am reducing the number of PEER factors to 20
library("peer")
for(i in 1:4)
{
  peerExprMat <- t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[, get(varsAsText[i])])
  model = PEER()
  PEER_setPhenoMean(model,as.matrix(peerExprMat))
  mod <- model.matrix(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[1]], get(varsAsText[i]), drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[1], get(varsAsText[i])]*theCpeProps_matched[get(varsAsText[i])])+genotypePcs$x[get(varsAsText[i]),1:3])
  PEER_setCovariates(model, mod[,3:6]) # include the purity estimate and the PCs as co-variates.
  PEER_setNk(model,20)
  PEER_update(model)
#   factors_correctedAnalysis_luminalA = PEER_getX(model)
  assign(paste("factors_correctedAnalysis_", varsAsText[i], sep=""), PEER_getX(model))
#   theFactors_correctedAnaysis_luminalA <- factors_correctedAnalysis_luminalA[,5:39] #' drop the real covariates (Pcs, purity), which are included in this object that peer returns.
  assign(paste("theFactors_correctedAnaysis_", varsAsText[i], sep=""), get(paste("factors_correctedAnalysis_", varsAsText[i], sep=""))[, 5:24])
  print(i)
}

# YOU ARE HERE: NEED TO MAKE SURE THIS RAN CORRECTLY AND RUN THIS FOR THE CONVENTIONAL APPROACH BELOW:::


#' Create PEER factors for each of the 4 subtypes, for the conventional approach.
for(i in 1:4)
{
  peerExprMat <- t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[, get(varsAsText[i])])
  model = PEER()
  PEER_setPhenoMean(model,as.matrix(peerExprMat))
  mod <- model.matrix(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[1]], get(varsAsText[i]), drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[1], get(varsAsText[i])]*theCpeProps_matched[get(varsAsText[i])])+genotypePcs$x[get(varsAsText[i]),1:3])
  PEER_setCovariates(model, mod[,4:6]) # use only the PCs from the genotype matrix as the co-variates.
  PEER_setNk(model,20)
  PEER_update(model)
  assign(paste("factors_conventionalAnalysis_", varsAsText[i], sep=""), PEER_getX(model))
  assign(paste("theFactors_conventionalAnaysis_", varsAsText[i], sep=""), get(paste("factors_conventionalAnalysis_", varsAsText[i], sep=""))[, 4:23])
  print(i)
}

# save these PEER factors for posterity
save(factors_conventionalAnalysis_aSamples, factors_correctedAnalysis_aSamples, factors_conventionalAnalysis_bSamples, factors_correctedAnalysis_bSamples, factors_conventionalAnalysis_basalSamples, factors_correctedAnalysis_basalSamples, factors_conventionalAnalysis_her2Samples, factors_correctedAnalysis_her2Samples, file=theRootDir %&% "Results/rDatas/subtypePeers.RData")
# load(file=theRootDir %&% "Results/rDatas/subtypePeers.RData")

#' 8 analyses here:
#' Conventional for lum A, B, basal, her2
#' Interaction for lum A, B, basal, her2
#' Run the association analysis. (should take about a day)
coefsList_conv_LuminalA <- list(length(snpsToGenesList))
coefsList_conv_LuminalB <- list(length(snpsToGenesList))
coefsList_conv_Basal <- list(length(snpsToGenesList))
coefsList_conv_Her2 <- list(length(snpsToGenesList))

coefsList_intInv_LuminalA <- list(length(snpsToGenesList)) # these will contain the cancer eQTL effect sizes.
coefsList_intInv_LuminalB <- list(length(snpsToGenesList))
coefsList_intInv_Basal <- list(length(snpsToGenesList))
coefsList_intInv_Her2 <- list(length(snpsToGenesList))
propInvA <- (1 - theCpeProps_matched[aSamples])
propInvB <- (1 - theCpeProps_matched[bSamples])
propInvBasal <- (1 - theCpeProps_matched[basalSamples])
propInvHer2 <- (1 - theCpeProps_matched[her2Samples])

# Note, when I was running this, I only had enough memory to run it 4 at a time.
for(i in 1:length(snpsToGenesList)) # [1] 615568
{
#   # conventional model
  coefsList_conv_LuminalA[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], aSamples,drop=FALSE])~genotypesMat_filt[names(snpsToGenesList)[i],aSamples ]+genotypePcs$x[aSamples,1:3]+factors_conventionalAnalysis_aSamples)))
  coefsList_conv_LuminalB[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], bSamples,drop=FALSE])~genotypesMat_filt[names(snpsToGenesList)[i], bSamples]+genotypePcs$x[bSamples,1:3]+factors_conventionalAnalysis_bSamples)))
  coefsList_conv_Basal[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], basalSamples,drop=FALSE])~genotypesMat_filt[names(snpsToGenesList)[i],basalSamples]+genotypePcs$x[basalSamples,1:3]+factors_conventionalAnalysis_basalSamples)))
  coefsList_conv_Her2[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], her2Samples,drop=FALSE])~genotypesMat_filt[names(snpsToGenesList)[i],her2Samples]+genotypePcs$x[her2Samples,1:3]+factors_conventionalAnalysis_her2Samples)))
 

  # inverse interaction model (to recover "cancer" effects)
  coefsList_intInv_LuminalA[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], aSamples,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i], aSamples]*propInvA)+genotypePcs$x[aSamples,1:3]+factors_correctedAnalysis_aSamples)))
  coefsList_intInv_LuminalB[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], bSamples,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i], bSamples]*propInvB)+genotypePcs$x[bSamples,1:3]+factors_correctedAnalysis_bSamples)))
  coefsList_intInv_Basal[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], basalSamples,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i],basalSamples]*propInvBasal)+genotypePcs$x[basalSamples,1:3]+factors_correctedAnalysis_basalSamples)))
  coefsList_intInv_Her2[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], her2Samples,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i],her2Samples]*propInvHer2)+genotypePcs$x[her2Samples,1:3]+factors_correctedAnalysis_her2Samples)))
  
  print(i)
}
 
#' Add the names to these.
names(coefsList_conv_LuminalA) <- names(snpsToGenesList)
names(coefsList_conv_LuminalB) <- names(snpsToGenesList)
names(coefsList_conv_Basal) <- names(snpsToGenesList)
names(coefsList_conv_Her2) <- names(snpsToGenesList)

names(coefsList_intInv_LuminalA) <- names(snpsToGenesList)
names(coefsList_intInv_LuminalB) <- names(snpsToGenesList)
names(coefsList_intInv_Basal) <- names(snpsToGenesList)
names(coefsList_intInv_Her2) <- names(snpsToGenesList)

dir.create(theRootDir %&% "data/forMatrixEQTL/output/", recursive=TRUE, showWarnings = FALSE)

save(coefsList_conv_LuminalA, coefsList_conv_LuminalB, coefsList_conv_Basal, coefsList_conv_Her2, file=theRootDir %&% "data/forMatrixEQTL/output/subtype_conv_CoefLists.RData") # Contains the first 8 coefLists.
save(coefsList_intInv_LuminalA, coefsList_intInv_LuminalB, coefsList_intInv_Basal, coefsList_intInv_Her2, file=theRootDir %&% "data/forMatrixEQTL/output/subtypeIntInvCoefLists.RData")

print(sessionInfo())