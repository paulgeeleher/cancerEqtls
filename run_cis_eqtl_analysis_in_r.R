#' This script runs the TCGA eQTL analysis in R. The script does 5 analyses:
#' Conventional With/wihtout PEER.
#' Interaction with/without PEER.
#' Interaction with PEER, but minus factors that are correlated with theCpeProps_matched.


# Original script: rEqtlAnalysis_FINAL.R

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

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

#' Calculate PEER factors from the normalized gene expression data (tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin).
#' NB: PEER factors are calculated in a way that is aware of the co-variates, hence I could just include the interaction term and proportions as co-variates in the analysis.
#' NB: Calculate PEERs separately for these and for the conventional analysis, as the co-variates change.
library("peer")
peerExprMat <- t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin)
model = PEER()
PEER_setPhenoMean(model,as.matrix(peerExprMat))
mod <- model.matrix(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[1]], ,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[1], ]*theCpeProps_matched)+genotypePcs$x[,1:3])
PEER_setCovariates(model, mod[,3:6]) # include the purity estimate and the PCs as co-variates.
PEER_setNk(model,35)
PEER_update(model)
factors_correctedAnalysis = PEER_getX(model)
theFactors_correctedAnaysis <- factors_correctedAnalysis[,5:39] #' drop the real covariates (Pcs, purity), which are included in this object that peer returns.

#' How correlated are the PEER factors and tumor purity? I want to test the effect of removing PEER factors that are correlated with tumor purity.
pValCorProp <- numeric()
estCor <- numeric()
for(i in 1:ncol(theFactors_correctedAnaysis))
{
  pValCorProp[i] <- cor.test(theFactors_correctedAnaysis[,i], theCpeProps_matched)$p.value
  estCor[i] <- cor.test(theFactors_correctedAnaysis[,i], theCpeProps_matched)$estimate
  print(i)
}
cbind(pValCorProp, estCor)
theFactors_correctedAnaysis_filtered <- theFactors_correctedAnaysis[, -which(pValCorProp < 0.05)] # 8 peer factors are removed because of correlation with the tumor purity, which I do not wish to remove.


#' create PEER factors for the uncorrected analysis.
modelPeerNoInt = PEER()
PEER_setPhenoMean(modelPeerNoInt,as.matrix(peerExprMat))
PEER_setCovariates(modelPeerNoInt, mod[,4:6]) # use only the PCs from the genotype matrix as the co-variates.
PEER_setNk(modelPeerNoInt,35)
PEER_update(modelPeerNoInt)
factors_NoInteractionAnalysis = PEER_getX(modelPeerNoInt)
theFactors_NoInteractionAnalysis <- factors_NoInteractionAnalysis[, 4:38]


#' Create PEER factors for an analysis where CPE estimates are arranged at random.
set.seed(12345)
randomCpe <- sample(theCpeProps_matched)
modelRandom = PEER()
PEER_setPhenoMean(modelRandom,as.matrix(peerExprMat))
modRand <- model.matrix(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[1]], ,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[1], ]*randomCpe)+genotypePcs$x[,1:3])
PEER_setCovariates(modelRandom, modRand[,3:6]) # include the purity estimate and the PCs as co-variates.
PEER_setNk(modelRandom,35)
PEER_update(modelRandom)
factors_randomCpeAnalysis = PEER_getX(modelRandom)
theFactors_randomCpeAnaysis <- factors_randomCpeAnalysis[,5:39] #' Drop the real covariates (Pcs, purity), which are included in this object that peer returns.

# save these (used in subsequent analysis)
dir.create(theRootDir %&% "Results/rDatas/", recursive=TRUE, showWarnings = FALSE)
save(theFactors_NoInteractionAnalysis, theFactors_correctedAnaysis, theFactors_correctedAnaysis_filtered, file=theRootDir %&% "Results/rDatas/thePeers.RData") # 



#' 5 analyses here:
#' Conventional With/wihtout PEER.
#' Interaction with/without PEER.
#' Interaction with PEER, but minus factors that are correlated with theCpeProps_matched

#' Calculate the PEER factors for the corrected analysis, i.e. include the interaction terms as co-variates
#' Run the association analysis.
#' There will be several analysis, one with no PEER factors and an interaction term, one that includes PEER factors and no interaction term, and one that includes PEER factors and the interaction term.
#' The GTEx consortium recommends 35 PEER factors for samples sizes >250, hence we have used 35 PEER factors.
coefsList_conv_noPeer <- list(length(snpsToGenesList))
coefsList_conv_Peer <- list(length(snpsToGenesList))
coefsList_int_noPeer <- list(length(snpsToGenesList))
coefsList_int_Peer <- list(length(snpsToGenesList))
coefsList_int_selectPeers <- list(length(snpsToGenesList))
# Should complete in 40hrs
for(i in 1:length(snpsToGenesList)) # [1] 615568
{
  coefsList_conv_noPeer[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], ,drop=FALSE])~genotypesMat_filt[names(snpsToGenesList)[i], ]+genotypePcs$x[,1:3])))
  coefsList_conv_Peer[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], ,drop=FALSE])~genotypesMat_filt[names(snpsToGenesList)[i], ]+genotypePcs$x[,1:3]+theFactors_NoInteractionAnalysis)))
  coefsList_int_noPeer[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], ,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i], ]*theCpeProps_matched)+genotypePcs$x[,1:3])))
  coefsList_int_Peer[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], ,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i], ]*theCpeProps_matched)+genotypePcs$x[,1:3]+theFactors_correctedAnaysis)))
  coefsList_int_selectPeers[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], ,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i], ]*theCpeProps_matched)+genotypePcs$x[,1:3]+theFactors_correctedAnaysis_filtered)))
# print(i)
}
names(coefsList_conv_noPeer) <- names(snpsToGenesList)
names(coefsList_conv_Peer) <- names(snpsToGenesList)
names(coefsList_int_noPeer) <- names(snpsToGenesList)
names(coefsList_int_Peer) <- names(snpsToGenesList)
names(coefsList_int_selectPeers) <- names(snpsToGenesList)

dir.create(theRootDir %&% "data/forMatrixEQTL/output/", recursive=TRUE, showWarnings = FALSE)
save(coefsList_conv_noPeer, coefsList_conv_Peer, coefsList_int_noPeer, coefsList_int_Peer, coefsList_int_selectPeers, file=theRootDir %&% "data/forMatrixEQTL/output/allTheCoefLists.RData")


# Also run these with the interaction terms reversed, so as to calculate the genotype P-value in the "cancer" cells rather than the normal cell.
coefsList_intInv_noPeer <- list(length(snpsToGenesList))
coefsList_intInv_Peer <- list(length(snpsToGenesList))
coefsList_intInv_selectPeers <- list(length(snpsToGenesList))
theCpeProps_matched_INV <- 1 - theCpeProps_matched
for(i in 1:length(snpsToGenesList)) # [1] 615568
{
  coefsList_intInv_noPeer[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], ,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i], ]*theCpeProps_matched_INV)+genotypePcs$x[,1:3])))
  coefsList_intInv_Peer[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], ,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i], ]*theCpeProps_matched_INV)+genotypePcs$x[,1:3]+theFactors_correctedAnaysis)))
  coefsList_intInv_selectPeers[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], ,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i], ]*theCpeProps_matched_INV)+genotypePcs$x[,1:3]+theFactors_correctedAnaysis_filtered)))
#   print(i)
}
names(coefsList_intInv_noPeer) <- names(snpsToGenesList)
names(coefsList_intInv_Peer) <- names(snpsToGenesList)
names(coefsList_intInv_selectPeers) <- names(snpsToGenesList)
save(coefsList_intInv_noPeer, coefsList_intInv_Peer, coefsList_intInv_selectPeers, file=theRootDir %&% "data/forMatrixEQTL/output/allInvCoefLists.RData")


coefsList_int_selectPeers <- coefsList_intInv_selectPeers


#' Finally run this for the randomly arranged CPE estimates. (Should take approx 8 hours)
coefList_intInv_Peer_RandomCpe <- list(length(snpsToGenesList))
for(i in 1:length(snpsToGenesList)) # [1] 615568
{
  coefList_intInv_Peer_RandomCpe[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], ,drop=FALSE])~(genotypesMat_filt[names(snpsToGenesList)[i], ]*randomCpe)+genotypePcs$x[,1:3]+theFactors_randomCpeAnaysis)))
#   print(i)
}

save(coefList_intInv_Peer_RandomCpe, randomCpe, file=theRootDir %&% "data/forMatrixEQTL/output/coefList_intInv_Peer_RandomCpe.RData")


print(sessionInfo())