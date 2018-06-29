## Run the analysis, but only on the 10% of samples with the most "cancer cell" content.

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

# Load the PEER factors.
load(file=theRootDir %&% "Results/rDatas/thePeers.RData") # theFactors_NoInteractionAnalysis, theFactors_correctedAnaysis, theFactors_correctedAnaysis_filtered, 

#' What are the 10% of samples with the highest cancer cell content?
top10pctInfo <- sort(theCpeProps_matched, decreasing=T)[1:as.integer((length(theCpeProps_matched) / 100)*10)]
length(top10pctInfo)
min(top10pctInfo)
median(top10pctInfo)
top10pct <- names(top10pctInfo)

#' Conventional and interaction anaylsis, on 10pct of samples with most cancer content.
coefsList_conv_Peer_10pct <- list(length(snpsToGenesList))
for(i in 1:length(snpsToGenesList)) # [1] 615568; Should complete in 6hrs?
{
  coefsList_conv_Peer_10pct[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList[[i]], top10pct,drop=FALSE])~genotypesMat_filt[names(snpsToGenesList)[i], top10pct]+genotypePcs$x[top10pct,1:3]+theFactors_NoInteractionAnalysis[top10pct, ])))
  print(i)
}
names(coefsList_conv_Peer_10pct) <- names(snpsToGenesList)

dir.create(theRootDir %&% "data/forMatrixEQTL/output/", recursive=TRUE, showWarnings = FALSE)
save(coefsList_conv_Peer_10pct, file=theRootDir %&% "data/forMatrixEQTL/output/coefList10pct.RData")











