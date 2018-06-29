##' Here, I will re-run the eQTL analysis for the previous Breast Cancer GWAS *significant* imputed Snps that were created on the michigan imputation server.

library("ggplot2")
library("reshape2")
library("RColorBrewer")

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

### NNBBB, Include the source of the imputed data script....
#' Load the imputed SNP data. 
load(file=theRootDir %&% "data/brcaGwasSigSnpsImputed.RData") # charGenotypesMat, numGenotypesMat

#' Load the original breast cancer data.
load(file=theRootDir %&% "Results/rDatas/gwasFilteredInputData_WithQNexprAndCpe.RData") # tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin, geneLocation_filt_qnDat, theCpeProps_matched, covariatesVec_filt, expressionMatNorm_filt, expressionMat_filt, genotypesMat_filt, genotypePcs, snpsLocation_filt, geneLocation_filt, 

#' Load the PEER factor data.
load(file=theRootDir %&% "Results/rDatas/thePeers.RData") # theFactors_NoInteractionAnalysis, theFactors_correctedAnaysis, theFactors_correctedAnaysis_filtered, 


## Map these new SNPs to Genes!
# First I have to find the HG18 location of these SNPs (I can get this from biomaRt: http://may2009.archive.ensembl.org/biomart/martview/)
# I am using HG18 to be consistent with the annotation provided by affymetrix for the SNP 6.0 array that were used in the original analysis (of the unimputed snps).
# write.table(colnames(charGenotypesMat), quote=F, row.names=F) # This line will output the data.
hg18BrcaGwasHits <- read.delim(theRootDir %&% "brcaGwasSigSnpsHg18Locs_fromBiomart.txt", as.is=T) # the resulting tsv file (from biomart) is saved here, load it. The last row contains duplicated snp on a contig: drop it
hg18BrcaGwasHits <- hg18BrcaGwasHits[!duplicated(hg18BrcaGwasHits[,1]), ]

#' Create a list mapping SNPs to genes, with a 1 megabase window around the gene.
hg18BrcaGwasHits[,4] <- hg18BrcaGwasHits[,3] - 500000 # lets do a megabase window around the gene.
hg18BrcaGwasHits[,4][hg18BrcaGwasHits[,4] < 0] <- 0
hg18BrcaGwasHits[,5] <- hg18BrcaGwasHits[,3] + 500000
colnames(hg18BrcaGwasHits)[c(4,5)] <- c("minusMB", "plusMB")
rownames(hg18BrcaGwasHits) <- hg18BrcaGwasHits[,1]

#' We will match snps to genes using the genomic ranges package in R, which will allow us to do this very quickly.
library(GenomicRanges) # this library should be able to do this fairly efficiently I think.
geneRanges <- GRanges(seqnames=Rle(geneLocation_filt_qnDat[,2]), ranges=IRanges(geneLocation_filt_qnDat[,3], geneLocation_filt_qnDat[,4]))
names(geneRanges) <- geneLocation_filt_qnDat[,1]
snpLocations <- GRanges(seqnames=Rle(hg18BrcaGwasHits[,2]), ranges=IRanges(hg18BrcaGwasHits[,4], hg18BrcaGwasHits[,5]))
names(snpLocations) <- hg18BrcaGwasHits[,1]

countOut <- findOverlaps(geneRanges,snpLocations, type="any", ignore.strand=T)
countOut_mat <- as.matrix(countOut)
snpSplit <- split(countOut_mat[,1], countOut_mat[,2])

snpsToGenesMat <- countOut_mat
snpsToGenesMat[,1] <- geneLocation_filt_qnDat[,1][countOut_mat[,1]]
snpsToGenesMat[,2] <- hg18BrcaGwasHits[,1][countOut_mat[,2]]
snpsToGenesList_Imputed_GWAS_Snps <- split(snpsToGenesMat[,1], snpsToGenesMat[,2])
length(snpsToGenesList_Imputed_GWAS_Snps)

# Correct the column names on the genotype matrix, and subset to the breast cancer.
impGwasSampleNames <- substring(colnames(numGenotypesMat), 9, 12)

# There are 32 samples for which there are no imputations.
length(theCpeProps_matched) # [1] 894
sum(impGwasSampleNames %in% names(theCpeProps_matched)) # [1] 862
numGenotypesMat_brca <- numGenotypesMat[, impGwasSampleNames %in% names(theCpeProps_matched)] # breast cancer genotypes.
colnames(numGenotypesMat_brca) <- substring(colnames(numGenotypesMat_brca), 9, 12)
charGenotypesMat_brca <- charGenotypesMat[, impGwasSampleNames %in% names(theCpeProps_matched)]
colnames(charGenotypesMat_brca) <- substring(colnames(charGenotypesMat_brca), 9, 12)
sampsWithImpGwasSnps <- colnames(numGenotypesMat_brca)



#' Calculate the MAFs for these SNPs. I may want to remove if its extremely low in certain cases, as the models for these will be insufficient.
theMaf <- numeric(nrow(numGenotypesMat_brca))
numSamps <- ncol(numGenotypesMat_brca)
for(i in 1:nrow(numGenotypesMat_brca))
{
  tab <- table(factor(numGenotypesMat_brca[i,], levels=c(0,1,2))) # Count the number of 0, 1, 2 in each row of the genotype matrix.
  f1 <- ((tab["0"] * 2) + tab["1"]) / (sum(tab)*2) # count the total number of one of the alleles and divide by the total number alleles (total number of (0, 1, 2) x 2), to get the allele frequency for that allele
  f2 <- 1 - f1
  theMaf[i] <- min(c(f1, f2))
  print(i)
}
names(theMaf) <- rownames(numGenotypesMat_brca)
theMaf_mappedToGenes <- theMaf[names(theMaf) %in% names(snpsToGenesList_Imputed_GWAS_Snps)]

#' Find the snps with a MAF of > 0.05
snpsToGenesList_Imputed_GWAS_Snps_MafFilt <- snpsToGenesList_Imputed_GWAS_Snps[names(theMaf_mappedToGenes[theMaf_mappedToGenes > 0.05])] # rs6507583

#' Run the eQTL analysis for "cancer" for the imputed SNPs in GWAS.
coefsList_intInv_Peer_GWAS_sigSnps <- list(length(snpsToGenesList_Imputed_GWAS_Snps_MafFilt))
coefsList_Bulk_Peer_GWAS_sigSnps <- list(length(snpsToGenesList_Imputed_GWAS_Snps_MafFilt))
coefsList_intNormal_Peer_GWAS_sigSnps <- list(length(snpsToGenesList_Imputed_GWAS_Snps_MafFilt))
theCpeProps_matched_INV <- 1 - theCpeProps_matched
eQtlPValue <- numeric()
for(i in 1:length(snpsToGenesList_Imputed_GWAS_Snps_MafFilt)) # [1] 615568
{
  # Cancer
  coefsList_intInv_Peer_GWAS_sigSnps[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList_Imputed_GWAS_Snps_MafFilt[[i]], sampsWithImpGwasSnps,drop=FALSE])~(numGenotypesMat_brca[names(snpsToGenesList_Imputed_GWAS_Snps_MafFilt)[i], sampsWithImpGwasSnps]*theCpeProps_matched_INV[sampsWithImpGwasSnps])+genotypePcs$x[sampsWithImpGwasSnps,1:3]+theFactors_correctedAnaysis[sampsWithImpGwasSnps, ])))
  
  # Bulk (interaction term and the term for cancer cell proportions have been removed and the PEER factors created for this data are used.)
  coefsList_Bulk_Peer_GWAS_sigSnps[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList_Imputed_GWAS_Snps_MafFilt[[i]], sampsWithImpGwasSnps,drop=FALSE])~numGenotypesMat_brca[names(snpsToGenesList_Imputed_GWAS_Snps_MafFilt)[i], sampsWithImpGwasSnps]+genotypePcs$x[sampsWithImpGwasSnps,1:3]+theFactors_correctedAnaysis[sampsWithImpGwasSnps, ])))
  
  # Normal (I have used theCpeProps_matched instead of theCpeProps_matched_INV, thus estimating the main effect for "normal", not "cancer"). Could actually be very intesting if we see anything here....
  coefsList_intNormal_Peer_GWAS_sigSnps[[i]] <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[snpsToGenesList_Imputed_GWAS_Snps_MafFilt[[i]], sampsWithImpGwasSnps,drop=FALSE])~(numGenotypesMat_brca[names(snpsToGenesList_Imputed_GWAS_Snps_MafFilt)[i], sampsWithImpGwasSnps]*theCpeProps_matched[sampsWithImpGwasSnps])+genotypePcs$x[sampsWithImpGwasSnps,1:3]+theFactors_correctedAnaysis[sampsWithImpGwasSnps, ])))
    
  print(i)
}
names(coefsList_intInv_Peer_GWAS_sigSnps) <- names(snpsToGenesList_Imputed_GWAS_Snps_MafFilt)
names(coefsList_Bulk_Peer_GWAS_sigSnps) <- names(snpsToGenesList_Imputed_GWAS_Snps_MafFilt)
names(coefsList_intNormal_Peer_GWAS_sigSnps) <- names(snpsToGenesList_Imputed_GWAS_Snps_MafFilt)



#' This function is used to prepare the results obtained from above, (i.e. the model coefficients) and convert them into usable matrices (below).
returnMatices <- function(theCoefList, rowsToGet, snpsToGenesList) # rowsToGet is the rows of the coef matrix that we want, e.g. the genotype, interaction term etc.
{
  # I'd like a matrix of snp x gene with p-values, estimtes, which I will have to extract from "coefsList".
  numTargetGenes <- sapply(snpsToGenesList, length)
  withMultipleTargets <- which(numTargetGenes != 1) # I will have to run this separately for snps with 1 or more than 1 target gene, because the function has returned objects of different classes ("listof", or "matrix" if only one response variable)
  singleTarget <- which(numTargetGenes == 1)

  # First run for snps with many possible gene targets
  # Initialize the output matrix, because this apparently speeds things up a lot...
  theMatSize <- dim(theCoefList[[withMultipleTargets[1]]][[1]])[1] # Keep 2, 3 and the interaction term (the genotype effect, the proportion effect, the interaction term)
  snpNamesOrd <- names(snpsToGenesList) # snp names
  allSnpsToGenes <- unlist(snpsToGenesList[withMultipleTargets]) # this will be used to alocate the correct amount of memory.
  
  pMat <- numeric(length(rowsToGet)*length(allSnpsToGenes)) # allocate memory up front cause this is 150x faster.
  dim(pMat) <- c(length(rowsToGet), length(allSnpsToGenes))
  betaMat  <- numeric(length(rowsToGet)*length(allSnpsToGenes)) # allocate memory up front cause this is 150x faster.
  dim(betaMat) <- c(length(rowsToGet), length(allSnpsToGenes))
  errMat  <- numeric(length(rowsToGet)*length(allSnpsToGenes)) # allocate memory up front cause this is 150x faster.
  dim(errMat) <- c(length(rowsToGet), length(allSnpsToGenes))
  
  colnamesVec <- character(length(allSnpsToGenes))
  cumCols <- 0
  
  for(i in withMultipleTargets) 
  {
    # get the p-values in a matrix
    pS <- array(unlist(theCoefList[[i]]), dim=c(theMatSize,4,length(theCoefList[[i]])))[rowsToGet, 4, ]
    cumColsNew <- cumCols + length(theCoefList[[i]])
    pMat[1:length(rowsToGet), (cumCols+1):(cumColsNew)] <- pS
    
    # Get the beta values in a matrix
    betas <- array(unlist(theCoefList[[i]]), dim=c(theMatSize,4,length(theCoefList[[i]])))[rowsToGet, 1, ]
    betaMat[1:length(rowsToGet), (cumCols+1):(cumColsNew)] <- betas
    
    # get the std Errors in a matrix
    errs <- array(unlist(theCoefList[[i]]), dim=c(theMatSize,4,length(theCoefList[[i]])))[rowsToGet, 2, ]
    errMat[1:length(rowsToGet), (cumCols+1):(cumColsNew)] <- errs
    
    colnamesVec[(cumCols+1):(cumColsNew)] <- paste(substring(names(theCoefList[[i]]), 10), "~", snpNamesOrd[i], sep="")
    cumCols <- cumColsNew
    print(i)
  }  
    
  colnames(pMat) <- colnamesVec
  colnames(betaMat) <- colnamesVec
  colnames(errMat) <- colnamesVec

  # pull out the results for snps with a single target
  allSnpsToGenes_sing <- unlist(snpsToGenesList[singleTarget]) # this will be used to alocate the correct amount of memory.
  pMat_sing <- numeric(length(rowsToGet)*length(allSnpsToGenes_sing)) # allocate memory up front cause this is 150x faster.
  dim(pMat_sing) <- c(length(rowsToGet), length(allSnpsToGenes_sing))
  betaMat_sing <- numeric(length(rowsToGet)*length(allSnpsToGenes_sing)) # allocate memory up front cause this is 150x faster.
  dim(betaMat_sing) <- c(length(rowsToGet), length(allSnpsToGenes_sing))
  errMat_sing <- numeric(length(rowsToGet)*length(allSnpsToGenes_sing)) # allocate memory up front cause this is 150x faster.
  dim(errMat_sing) <- c(length(rowsToGet), length(allSnpsToGenes_sing))
  colnamesVec_sing <- character(length(allSnpsToGenes_sing))
  colNum <- 1
  for(i in singleTarget)
  {
    colnamesVec_sing[colNum] <- paste(snpsToGenesList[[i]], "~", names(snpsToGenesList[i]), sep="")
    pMat_sing[,colNum] <- theCoefList[[i]][rowsToGet, 4]
    betaMat_sing[,colNum] <- theCoefList[[i]][rowsToGet, 1]
    errMat_sing[,colNum] <- theCoefList[[i]][rowsToGet, 2]
    colNum <- colNum + 1
#     print(i)
  }
  colnames(pMat_sing) <- colnamesVec_sing
  colnames(betaMat_sing) <- colnamesVec_sing
  colnames(errMat_sing) <- colnamesVec_sing
  pMatAll <- cbind(pMat, pMat_sing)
  betaMatAll <- cbind(betaMat, betaMat_sing)
  errMatAll <- cbind(errMat, errMat_sing)
  
  return(list(pMatAll=pMatAll, betaMatAll=betaMatAll, errMatAll=errMatAll))
}

#' Process cancer results using function above.
rowsToGet <- c(2, 3, 42)
theData_Cancer <- returnMatices(coefsList_intInv_Peer_GWAS_sigSnps, rowsToGet, snpsToGenesList_Imputed_GWAS_Snps_MafFilt)

#' Bulk
rowsToGet <- 2
theData_bulk <- returnMatices(coefsList_Bulk_Peer_GWAS_sigSnps, rowsToGet, snpsToGenesList_Imputed_GWAS_Snps_MafFilt)

#' Normal
rowsToGet <- c(2, 3, 42)
theData_Normal <- returnMatices(coefsList_intNormal_Peer_GWAS_sigSnps, rowsToGet, snpsToGenesList_Imputed_GWAS_Snps_MafFilt)

sum(p.adjust(theData_bulk$pMatAll[1,], method="BH") < 0.05) # [1] 24
sigBulkEqtls <- names(which(p.adjust(theData_bulk$pMatAll[1,], method="BH") < 0.05))
fdrs_sigBulkEqtls_gwasSigSnps <- p.adjust(theData_bulk$pMatAll[1,], method="BH")[sigBulkEqtls]
fdrs_bulk_gwasSigSnps <- p.adjust(theData_bulk$pMatAll[1,], method="BH")
fdr_cancer_gwasSigSnps <- p.adjust(theData_Cancer$pMatAll[1,], method="BH")
sum(p.adjust(theData_Cancer$pMatAll[1,], method="BH") < 0.05) # [1] 9
sum(p.adjust(theData_Normal$pMatAll[1,], method="BH") < 0.05) # [1] 1

#' Are there any new eQTLs identifed in cancer?
sigCancerEqtls <- names(which(p.adjust(theData_Cancer$pMatAll[1,], method="BH") < 0.05))
nomSigCancerEqtls <- names(which(theData_Cancer$pMatAll[1,] < 0.05))
fdrsCancer_forSigBulkEqtls_gwasSigSnps <- p.adjust(theData_Cancer$pMatAll[1,], method="BH")[sigBulkEqtls]
minSigPCancer <- max(theData_Cancer$pMatAll[1,sigCancerEqtls]) # P-value where the FDR drops under 0.05

length(unique(do.call(cbind, strsplit(colnames(theData_Cancer$pMatAll), "~"))[2,])) # [1] 81
ncol(theData_Cancer$pMatAll) # [1] 565
length(sort(table(do.call(cbind, strsplit(sigBulkEqtls, "~"))[2,]))) # [1] 16

sort(table(do.call(cbind, strsplit(sigBulkEqtls, "~"))[2,]), decreasing=T)

sigSnpsBulk <- names(table(do.call(cbind, strsplit(sigBulkEqtls, "~"))[2,]))

sort(table(do.call(cbind, strsplit(sigCancerEqtls, "~"))[2,]), decreasing=T)

sigSnpsCancer <- names(table(do.call(cbind, strsplit(sigCancerEqtls, "~"))[2,]))
noLongerSigSnps <- sigSnpsBulk[!sigSnpsBulk %in% sigSnpsCancer]
nomSigSnpsCancer <- names(table(do.call(cbind, strsplit(nomSigCancerEqtls, "~"))[2,]))
noLongerEvenNomSigSnps <- sigSnpsBulk[!sigSnpsBulk %in% nomSigSnpsCancer]


#' These are the "not even nominally significant in cancer" GWAS-sig/bulk-sig eQTLs
print(a[(a[,2] > 0.05), ])



#' MAKE A TABLE FOR THESE SNPS.
# data from reRunQtlAnalysisForImputedSnps.R
gwasSigEqtl_sigInBulk_notNominallySigInCancer <- c("EIF1AD~rs3903072", "SNX32~rs3903072", "SSBP4~rs4808801", "TOM1L1~rs6504950", "CYP51A1~rs6964587", "ARHGEF5~rs720475", "ANKLE1~rs8170", "OCEL1~rs8170", "PLVAP~rs8170")
gwasSigEqtl_sigInBulk_cancerSigEqtls <- c('CDYL2~rs13329835', 'RANBP9~rs204247', 'CTSW~rs3903072', 'LRRC25~rs4808801', 'ATG10~rs7707921', 'ATP6AP1L~rs7707921', 'RPS23~rs7707921', 'C5orf35~rs889312', 'TOX3~rs3803662')


#' These annotations files and create a map of variant ID -> RS ID and ensemble ID to gene name. (they are used below as the GTEx data is not provided with RSids or gene symbols)
ensToGeneSymTab <- read.delim(theRootDir %&% "gencode.v19.genes.v6p.patched_contigs.parsed.txt", as.is=T)
varIdToRsIdTab <- read.delim(theRootDir %&% "GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.txt", as.is=T)

#' Map ensembl ids to gene symbols
geneSyms <- ensToGeneSymTab$gene_name
names(geneSyms) <- ensToGeneSymTab$gene_id
length(geneSyms) # [1] 56238
length(unique(geneSyms)) # [1] 54301
length(unique(names(geneSyms))) # [1] 56238

#' Map GTEx VariantIds to rs ids
geneEns <- varIdToRsIdTab$RS_ID_dbSNP142_CHG37p13
names(geneEns) <- varIdToRsIdTab$VariantID
length(geneEns) # [1] 11552519
length(unique(geneEns)) # [1] 11545740
length(unique(names(geneEns))) # [1] 11552513

#' Get this variant ID -> RS ID map for the breast cancer GWAS SNPs.
gwasSnpsWithEqtl <- names(snpsToGenesList_Imputed_GWAS_Snps_MafFilt)
varIdToRsIdTab_gwasSnps <- varIdToRsIdTab[which(varIdToRsIdTab$RS_ID_dbSNP142_CHG37p13 %in% gwasSnpsWithEqtl), ]
gwasVariants_gtexVariantId <- varIdToRsIdTab_gwasSnps[, "VariantID"]

#' Load the GTEx fibroblast data. Then annotate with gene symbols and rs ids
fibroblastGtex <- read.delim(theRootDir %&% "data/Cells_Transformed_fibroblasts_Analysis.v6p.all_snpgene_pairs.txt", as.is=T)
fibroblastGtex_gwasSigVars <- fibroblastGtex[fibroblastGtex[, "variant_id"] %in% gwasVariants_gtexVariantId, ]
fibroblastGtex_gwasSigVars$geneIdSym <- geneSyms[fibroblastGtex_gwasSigVars$gene_id]
fibroblastGtex_gwasSigVars$rsId <- geneEns[fibroblastGtex_gwasSigVars$variant_id]
fibroblastGtex_gwasSigVars$geneToSnps <- paste(fibroblastGtex_gwasSigVars$geneIdSym, "~", fibroblastGtex_gwasSigVars$rsId, sep="")
rownames(fibroblastGtex_gwasSigVars) <- fibroblastGtex_gwasSigVars$geneToSnps

#' GTEx breast data
breastResultsGtex <- read.delim(theRootDir %&% "data/Breast_Mammary_Tissue_Analysis.v6p.all_snpgene_pairs.txt", as.is=T)
breastGtex_gwasSigVars <- breastResultsGtex[breastResultsGtex[, "variant_id"] %in% gwasVariants_gtexVariantId, ]
breastGtex_gwasSigVars$geneIdSym <- geneSyms[breastGtex_gwasSigVars$gene_id]
breastGtex_gwasSigVars$rsId <- geneEns[breastGtex_gwasSigVars$variant_id]
breastGtex_gwasSigVars$geneToSnps <- paste(breastGtex_gwasSigVars$geneIdSym, "~", breastGtex_gwasSigVars$rsId, sep="")
rownames(breastGtex_gwasSigVars) <- breastGtex_gwasSigVars$geneToSnps

#' GTEx lcl data
lclGtex <- read.delim(theRootDir %&% "data/Cells_EBV-transformed_lymphocytes_Analysis.v6p.all_snpgene_pairs.txt", as.is=T)
lclGtex_gwasSigVars <- lclGtex[lclGtex[, "variant_id"] %in% gwasVariants_gtexVariantId, ]
lclGtex_gwasSigVars$geneIdSym <- geneSyms[lclGtex_gwasSigVars$gene_id]
lclGtex_gwasSigVars$rsId <- geneEns[lclGtex_gwasSigVars$variant_id]
lclGtex_gwasSigVars$geneToSnps <- paste(lclGtex_gwasSigVars$geneIdSym, "~", lclGtex_gwasSigVars$rsId, sep="")
rownames(lclGtex_gwasSigVars) <- lclGtex_gwasSigVars$geneToSnps


#' Make the table
gwasSigTab <- data.frame(P_cancer=theData_Cancer$pMatAll[1,sigBulkEqtls], beta_cancer=theData_Cancer$betaMatAll[1,sigBulkEqtls], fdr_cancer=fdr_cancer_gwasSigSnps[sigBulkEqtls], P_bulk=theData_bulk$pMatAll[1,sigBulkEqtls], fdr_bulk=fdrs_sigBulkEqtls_gwasSigSnps, P_GTEx_fibroblast=fibroblastGtex_gwasSigVars[sigBulkEqtls, "pval_nominal"], Slope_GTEx_fibroblast=fibroblastGtex_gwasSigVars[sigBulkEqtls, "slope"], P_GTEx_breast=breastGtex_gwasSigVars[sigBulkEqtls, "pval_nominal"], Slope_GTEx_breast=breastGtex_gwasSigVars[sigBulkEqtls, "slope"], P_GTEx_Lcl=lclGtex_gwasSigVars[sigBulkEqtls, "pval_nominal"], Slope_GTEx_Lcl=lclGtex_gwasSigVars[sigBulkEqtls, "slope"])
gwasSigTab_byCancerP <- gwasSigTab[order(gwasSigTab[, "P_cancer"]), ]
allGtexTab <- rbind(rbind(fibroblastGtex_gwasSigVars, breastGtex_gwasSigVars), lclGtex_gwasSigVars)

#' NB: As previously, I Need to flip some of the beta values for these eQLTs as I've used a different reference allele in the imputed TCGA data as had been used in GTEx.
# Find reference and alternate allele for GTEx for each snp
findAllels <- do.call(rbind, strsplit(allGtexTab[, "variant_id"], "_"))
refAllelGtex <- findAllels[,3]
names(refAllelGtex) <- allGtexTab[, "rsId"]
altAlleleGtex <- findAllels[,4]
names(altAlleleGtex) <- allGtexTab[, "rsId"]
refAllelGtex_noDups <- refAllelGtex[!duplicated(names(refAllelGtex))]
altAllelGtex_noDups <- altAlleleGtex[!duplicated(names(altAlleleGtex))]

# Find what has been used as ref and alt in our imputed TCGA data.
snpsGwasGtex <- names(altAllelGtex_noDups)
charGenotypesMat_filt <- charGenotypesMat[snpsGwasGtex, ]
numGenotypesMat_filt <- numGenotypesMat[snpsGwasGtex, ]
refAlleleTcgaImputed <- character()
altAlleleTcgaImputed <- character()
for(i in 1:nrow(charGenotypesMat_filt))
{
  refAlleleTcgaImputed[i] <- substr(charGenotypesMat_filt[i, which(numGenotypesMat_filt[i,] == 0)[1]], 1, 1)
  altAlleleTcgaImputed[i] <- substr(charGenotypesMat_filt[i, which(numGenotypesMat_filt[i,] == 2)[1]], 1, 1)
  print(i)
}
names(refAlleleTcgaImputed) <- rownames(charGenotypesMat_filt)
names(altAlleleTcgaImputed) <- rownames(charGenotypesMat_filt)
sum(refAlleleTcgaImputed == refAllelGtex_noDups) # [1] 49
sum(altAlleleTcgaImputed == refAllelGtex_noDups) # [1] 32

snpsNeedFlipping <- names(refAlleleTcgaImputed)[!refAlleleTcgaImputed == refAllelGtex_noDups]

# Now Flip the betas for the eQTLs that contain these flip SNPs
cancerBetasAll <- theData_Cancer$betaMatAll[1,]
bulkBetasAll <- theData_bulk$betaMatAll[1,]
allEqtlSnpsOrd <- do.call(cbind, strsplit(names(cancerBetasAll), "~"))[2,]
toFlipIndex <- which(allEqtlSnpsOrd %in% snpsNeedFlipping)
cancerBetasAll_flipped <- cancerBetasAll
cancerBetasAll_flipped[toFlipIndex] <- (cancerBetasAll_flipped[toFlipIndex] * -1)
bulkBetasAll_flipped <- bulkBetasAll
bulkBetasAll_flipped[toFlipIndex] <- (bulkBetasAll_flipped[toFlipIndex] * -1)

# Add these flipped Betas to the table.
gwasSigTab <- data.frame(P_cancer=theData_Cancer$pMatAll[1,sigBulkEqtls], beta_cancer=cancerBetasAll_flipped[sigBulkEqtls], fdr_cancer=fdr_cancer_gwasSigSnps[sigBulkEqtls], P_bulk=theData_bulk$pMatAll[1,sigBulkEqtls], fdr_bulk=fdrs_sigBulkEqtls_gwasSigSnps, bulkBetasAll_flipped[sigBulkEqtls], P_GTEx_fibroblast=fibroblastGtex_gwasSigVars[sigBulkEqtls, "pval_nominal"], Slope_GTEx_fibroblast=fibroblastGtex_gwasSigVars[sigBulkEqtls, "slope"], P_GTEx_breast=breastGtex_gwasSigVars[sigBulkEqtls, "pval_nominal"], Slope_GTEx_breast=breastGtex_gwasSigVars[sigBulkEqtls, "slope"], P_GTEx_Lcl=lclGtex_gwasSigVars[sigBulkEqtls, "pval_nominal"], Slope_GTEx_Lcl=lclGtex_gwasSigVars[sigBulkEqtls, "slope"])

bulkQtlsSort <- names(sort(theData_bulk$pMatAll[1,]))
suppTab_allGWAS_snpEqtls_sortBy_bulkP <- data.frame(P_bulk=theData_bulk$pMatAll[1,bulkQtlsSort], fdr_bulk=fdrs_bulk_gwasSigSnps[bulkQtlsSort], beta_bulk=bulkBetasAll_flipped[bulkQtlsSort], P_cancer=theData_Cancer$pMatAll[1,bulkQtlsSort], beta_cancer=cancerBetasAll_flipped[bulkQtlsSort], fdr_cancer=fdr_cancer_gwasSigSnps[bulkQtlsSort], P_GTEx_fibroblast=fibroblastGtex_gwasSigVars[bulkQtlsSort, "pval_nominal"], Slope_GTEx_fibroblast=fibroblastGtex_gwasSigVars[bulkQtlsSort, "slope"], P_GTEx_breast=breastGtex_gwasSigVars[bulkQtlsSort, "pval_nominal"], Slope_GTEx_breast=breastGtex_gwasSigVars[bulkQtlsSort, "slope"], P_GTEx_Lcl=lclGtex_gwasSigVars[bulkQtlsSort, "pval_nominal"], Slope_GTEx_Lcl=lclGtex_gwasSigVars[bulkQtlsSort, "slope"])

write.csv(suppTab_allGWAS_snpEqtls_sortBy_bulkP, file=theRootDir %&% "suppTab_allGWAS_snpEqtls_sortBy_bulkP.csv")

notEvenNomSigInCancer <- gwasSigTab[gwasSigTab[, "P_cancer"] > 0.05, ]
fdrSigInCancer <- gwasSigTab[gwasSigTab[, "fdr_cancer"] < 0.05, ]

fdrSigInCancer[, c("P_bulk", "bulkBetasAll_flipped.sigBulkEqtls.", "Slope_GTEx_fibroblast", "P_GTEx_fibroblast")]
fdrSigInCancer[, c("P_bulk", "P_GTEx_fibroblast", "P_GTEx_breast", "P_GTEx_Lcl")]
fdrSigInCancer[, c("bulkBetasAll_flipped.sigBulkEqtls.", "Slope_GTEx_fibroblast", "Slope_GTEx_breast", "Slope_GTEx_Lcl")]
notEvenNomSigInCancer[, c("P_bulk", "bulkBetasAll_flipped.sigBulkEqtls.", "Slope_GTEx_fibroblast", "P_GTEx_fibroblast")]


#' Now make a heatmap of the "notEvenNomSigInCancer" associations (Figure 4 b)
notNomHeatMat <- notEvenNomSigInCancer[, c("P_cancer", "P_bulk", "P_GTEx_fibroblast", "P_GTEx_breast", "P_GTEx_Lcl")]
notNomSlopeHeatMat <- notEvenNomSigInCancer[, c("beta_cancer", "bulkBetasAll_flipped.sigBulkEqtls.", "Slope_GTEx_fibroblast", "Slope_GTEx_breast", "Slope_GTEx_Lcl")]
datLong <- melt(data.matrix(-log10(notNomHeatMat)))
colnames(notNomSlopeHeatMat) <- c("Cancer", "Bulk Tumor", "GTEx Fibroblast", "GTEx Breast", "GTEx LCL")
slopeDatLong <- melt(data.matrix(notNomSlopeHeatMat))
allDatLong <- cbind(slopeDatLong, datLong[,3])
colnames(allDatLong) <- c("eQTL", "Source", "EffectSize", "P")
allDatLong_filt <- allDatLong[allDatLong[, "P"] > -log10(0.05), ] # Only display the effect size for "significant" effects, i.e. things that are at least P < 0.05 in lcl, breast, fibroblast, tumor (Cancer will never be P < 0.05 here.).
colnames(allDatLong_filt) <- c("eQTL", "Source", "EffectSize", "P")
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")), space="Lab")
dir.create(theRootDir %&% "paper/figures/figure4/", recursive=T, showWarnings=F)
svg(file=theRootDir %&% "paper/figures/figure4/fig4b_1.svg", width=4, height=3)
ggplot(data = allDatLong_filt, aes(x=Source, y=eQTL)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_point(aes(color=EffectSize, size=P)) + scale_color_gradient2(low='#2166ac',mid="white", high="#b2182b",midpoint=0) #+ opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())
dev.off()



#' Now make a heatmap of the FDR sig associations in cancer (Supplementary Figure)
fdrSigHeatMat <- fdrSigInCancer[, c("P_cancer", "P_bulk", "P_GTEx_fibroblast", "P_GTEx_breast", "P_GTEx_Lcl")]
fdrSigSlopeHeatMat <- fdrSigInCancer[, c("beta_cancer", "bulkBetasAll_flipped.sigBulkEqtls.", "Slope_GTEx_fibroblast", "Slope_GTEx_breast", "Slope_GTEx_Lcl")]
datLong <- melt(data.matrix(-log10(fdrSigHeatMat)))
colnames(fdrSigSlopeHeatMat) <- c("Cancer", "Bulk Tumor", "GTEx Fibroblast", "GTEx Breast", "GTEx LCL")
slopeDatLong <- melt(data.matrix(fdrSigSlopeHeatMat))
allDatLong <- cbind(slopeDatLong, datLong[,3])
colnames(allDatLong) <- c("eQTL", "Source", "EffectSize", "P")
allDatLong_filt <- allDatLong[allDatLong[, "P"] > -log10(0.05), ]
colnames(allDatLong_filt) <- c("eQTL", "Source", "EffectSize", "P")
allDatLong_filt <- allDatLong_filt[!is.na(allDatLong_filt["eQTL"]), ]
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")), space="Lab")

# C5orf35 is annotated as SETD9 in the gTEX data, fix this here and make sure these are included in the figure.
# SETD9 == ENSG00000155542.7
# rs889312 == 5_56031884_C_A_b37
setd9InfoGtex_betas <- bulkBetasAll_flipped["SETD9~rs889312"]
setdEqtlRow_brca <- intersect(which(breastResultsGtex[, "gene_id"] == "ENSG00000155542.7"), which(breastResultsGtex[, "variant_id"] == "5_56031884_C_A_b37"))
setdEqtlRow_fibro <- intersect(which(fibroblastGtex[, "gene_id"] == "ENSG00000155542.7"), which(fibroblastGtex[, "variant_id"] == "5_56031884_C_A_b37"))
setdEqtlRow_lcl <- intersect(which(lclGtex[, "gene_id"] == "ENSG00000155542.7"), which(lclGtex[, "variant_id"] == "5_56031884_C_A_b37"))
breastResultsGtex[setdEqtlRow_brca, ]
fibroblastGtex[setdEqtlRow_fibro, ]
lclGtex[setdEqtlRow_lcl, ]
setd9DF <- data.frame(eQTL=rep("C5orf35~rs889312", 3), Source=c("GTEx Breast", "GTEx Fibroblast", "GTEx LCL"), EffectSize=c(breastResultsGtex[setdEqtlRow_brca, "slope"], fibroblastGtex[setdEqtlRow_fibro, "slope"], lclGtex[setdEqtlRow_lcl, "slope"]), P=c(breastResultsGtex[setdEqtlRow_brca, "pval_nominal"], fibroblastGtex[setdEqtlRow_fibro, "pval_nominal"], lclGtex[setdEqtlRow_lcl, "pval_nominal"]))
allDatLong_filt_setd9 <- rbind(allDatLong_filt, c(setd9DF))

#' This is a supplementary figure (with the cancer eQTL effects in tcga and gtex.)
svg(file=theRootDir %&% "fig4c.svg", width=4.5, height=3.5)
ggplot(data = allDatLong_filt_setd9, aes(x=Source, y=eQTL)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_point(aes(color=EffectSize, size=P)) + scale_color_gradient2(low='#2166ac',mid="white", high="#b2182b",midpoint=0) #+ opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())
dev.off()


#' These are plot for RANBP9 vs rs204247. These are supplementary figures.
svg(file=theRootDir %&% "RANBP9_rs204247_boxplot.svg", width=3, height=3)
boxplot(split(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin["RANBP9", sampsWithImpGwasSnps], charGenotypesMat_brca["rs204247", sampsWithImpGwasSnps]), vertical=T, pch=20, las=1, ylab="RANBP9 expression", xlab="rs204247", bty="l", cex.axis=.8)
dev.off()

svg(file=theRootDir %&% "RANBP9_rs204247_stripchart.svg", width=3, height=3)
plot((numGenotypesMat_brca["rs204247", sampsWithImpGwasSnps]+rnorm(length(numGenotypesMat_brca["rs204247", sampsWithImpGwasSnps]), 0, .1)), tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin["RANBP9", sampsWithImpGwasSnps], col="#00000044", pch=20, las=1, ylab="RANBP9 expression", xlab="rs204247", bty="l", xaxt = "n", cex.axis=.8) 
axis(1, at=0:2, labels=c("A/A","A/G","G/G"), cex.axis=.8)
abline(coef=coef(summary(lm(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin["RANBP9", sampsWithImpGwasSnps]~numGenotypesMat_brca["rs204247", sampsWithImpGwasSnps])))[,1], col="#e31a1c")
dev.off()

#' For the eQTLs significant in the Bulk tissue, what do their p-values look like in the Cancer and Normal Components
print(sort(theData_bulk$pMatAll[1,sigBulkEqtls]))
print(sort(theData_Cancer$pMatAll[1,sigBulkEqtls]))
sort(theData_Normal$pMatAll[1,sigBulkEqtls])

#' What does the interaction term look like #### interaction p-value is the same for cancer and normal
print(sort(theData_Cancer$pMatAll[3,sigBulkEqtls]))

#' Create Figure 4a.
svg(theRootDir %&% "paper/figures/figure4/fig4a.svg", width=3, height=4)
plot(-log10(theData_bulk$pMatAll[1,sigBulkEqtls]), -log10(theData_Cancer$pMatAll[1,sigBulkEqtls]), pch=20, col="#00000099", las=1, bty="l", cex.axis=0.8, ylab=expression('-log'[10]*'('*italic("P")*')'*' Interaction model'), xlab=expression('-log'[10]*'('*italic("P")*')'*' Conventional model'))
abline(h=-log10(0.05), col="#4daf4a", lty=2)
abline(h=-log10(minSigPCancer), col="#984ea3", lty=2)
dev.off()


print(sessionInfo())



