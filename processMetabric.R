#' This script reads the METABRIC genotype matrix (output by Affymetrix Genotyping console) and expression data (obatained from EGA) and processes the data, preforming the required filtering, creating the required covaraites, purity estimates and saving out outputting these data.

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

# read the expression data and convert the illumina Ids to gene symbols.
# I did this locally as I couldnt' get the R package to wrok on bionimbus
library("illuminaHumanv4.db") # this needs to be installled. It gave an uninterpretable error on bionimbus
metabricExpressiondDisc <- read.delim(theRootDir %&% "Metabric/EGAD00010000210/_ega-box-04_discovery_ExpressionMatrix.txt")
ilToSym = toTable(illuminaHumanv4SYMBOL)
rownames(ilToSym) <- ilToSym[,1]
ilToSymNoDups <- ilToSym[!duplicated(ilToSym[,2]),]    # drop things that are duplicated
commonIds <- rownames(ilToSymNoDups)[rownames(ilToSymNoDups) %in% rownames(metabricExpressiondDisc)]
metabricExpressiondDiscNodups <- metabricExpressiondDisc[commonIds, ]
rownames(metabricExpressiondDiscNodups) <- ilToSymNoDups[commonIds, 2]
write.csv(metabricExpressiondDiscNodups, file="/home/pgeeleher/Downloads/metabricExpression/EGAD00010000210/metabricDiscovery_withGeneSyms.csv")

# Read the csv file (on bionimbus) that was create locally then transferred.
metabricExpressiondDiscNodups <- read.csv(theRootDir %&% "Metabric/metabricDiscovery_withGeneSyms.csv")
rownames(metabricExpressiondDiscNodups) <- metabricExpressiondDiscNodups[,1]
metabricExpressiondDiscNodups[, 1] <- NULL
metabricExpressiondDiscNodups <- data.matrix(metabricExpressiondDiscNodups)


## NOW I WOULD LIKELY NEED TO DO THE SAME THINGS TO THIS EXPRESSION MATRIX THAT I DID TO THE TCGA DATA. MAPPING TO NORMAL DISTRIBITON ETC.
# The normalization steps taken by GTEx were to map to a normal distribution and remove genes not expressed in a higher percentage of samples (>75%)
# Seeing as this is a comparative analysis, I should rmeove any genes that weren't also in TCGA, then map the remaining genes to a normal distribiton (as in the other analysis)

# Load the processed TCGA data
#' Load the filtered subsetted files needed for the GWAS analysis
load(file=theRootDir %&% "Results/rDatas/gwasFilteredInputData_WithQNexprAndCpe.RData") # tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin, geneLocation_filt_qnDat, theCpeProps_matched, covariatesVec_filt, expressionMatNorm_filt, expressionMat_filt, genotypesMat_filt, genotypePcs, snpsLocation_filt, geneLocation_filt, 
rm(list=c("expressionMatNorm_filt", "expressionMat_filt", "genotypePcs", "geneLocation_filt"))
genesInTcga <- rownames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin)
hist(metabricExpressiondDiscNodups[,1], breaks=500)
sum(genesInTcga %in% rownames(metabricExpressiondDiscNodups)) #[1] 13805 ; most of the 19201 metabric genes are shared with TCGA
metaTcgaCommonGenes <- genesInTcga[genesInTcga %in% rownames(metabricExpressiondDiscNodups)]

# These histograms are encouraging, i.e. it seems we lost many of the low expressed genes, but keep almost all the highly expressed genes.
# CAUTION: There are still quite a lot of "low" expressed genes remaining, which I guess would be expected, given these arrays are almost certainly not as sensitive as RNA-seq
hist(as.numeric(metabricExpressiondDiscNodups), breaks=100, col="#00ff0066")
hist(as.numeric(metabricExpressiondDiscNodups[metaTcgaCommonGenes,]), breaks=100, col="#0000ff66", add=T)


# Map these to a normal distribtion, as we did for TCGA (following the guidelines of GTEx), i.e. the distribution of every gene across samples should be normal.
expresssionMatMetaFin <- metabricExpressiondDiscNodups[metaTcgaCommonGenes,]
expressionMatNormMeta <- expresssionMatMetaFin
for(i in 1:nrow(expresssionMatMetaFin))
{
  expressionMatNormMeta[i,] <- qnorm((rank(expresssionMatMetaFin[i,]) / (length(expresssionMatMetaFin[i,])+1)))
  print(i)
}

# Get CPE estimates for metabric. Do this by lasso regression, test by a hold out set in TCGA too.
library(glmnet)
set.seed(12345)

# Test on a hold out set of 100 samples in TCGA
holdOutSamples <- sample(length(theCpeProps_matched), 100)
trainM = t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[, -holdOutSamples])
trainPtype = theCpeProps_matched[-holdOutSamples]
glmModel <- cv.glmnet(trainM, trainPtype) # lasso model
getLambdas <- glmModel$lambda.min # extract the beta values
colmn <-  which(glmModel$lambda == getLambdas)
a <- glmModel$glmnet.fit$beta[, colmn]
geneList <- a[a!=0]
predMat <- t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[, holdOutSamples]) # get the prediction
predictOut <- predict(glmModel, predMat, s = "lambda.min")
cor.test(predictOut, theCpeProps_matched[holdOutSamples]) # cor = 0.89, P < 2.2e-16

# Q: A double check: Are the CPE estimates correlated with the two cytolytic t-cell marker genes from Rooney et al? Yes, will this also be the case in METABRIC? Ans: yes it is.
cor(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin["GZMA", ], theCpeProps_matched) # [1] -0.7019058
cor(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin["PRF1", ], theCpeProps_matched) # [1] -0.6595709

# Now fit the model on all of TCGA and apply to METABRIC. 
# For learning, only keep the genes that are also METABRIC GENES.
set.seed(12345)
trainM = t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[metaTcgaCommonGenes, ])
trainPtype = theCpeProps_matched
glmModel <- cv.glmnet(trainM, trainPtype) # lasso model
getLambdas <- glmModel$lambda.min # extract the beta values
colmn <-  which(glmModel$lambda == getLambdas)
a <- glmModel$glmnet.fit$beta[, colmn]
geneList <- a[a!=0] # the genes used in the prediction...
predMat <- t(expressionMatNormMeta) # get the prediction
cpePredictionsMeta <- predict(glmModel, predMat, s = "lambda.min")[,1]
cor.test(cpePredictionsMeta, expressionMatNormMeta["GZMA", ]) # -0.7968994
cor.test(cpePredictionsMeta, expressionMatNormMeta["PRF1", ]) # -0.7828933

# For consistency, I would also like to map the predicted values to the same quantiles as the TCGA CPE data. (The lasso prediction has compressed the variance of the predicted distribution, and if the variance has been slightly underestimated this could slightly affect the downstream results)
quantilesVecTcgaSort <- sort(sample(theCpeProps_matched, length(cpePredictionsMeta), replace=T)) # a vector of suitable length, drawn from the TCGA cpe distribiton
names(quantilesVecTcgaSort) <- names(sort(cpePredictionsMeta))
cpePredictionsMetaFin <- quantilesVecTcgaSort[names(cpePredictionsMeta)]
cor.test(cpePredictionsMetaFin, expressionMatNormMeta["GZMA", ]) # -0.7956933, P = 7.047649e-219
cor.test(cpePredictionsMetaFin, expressionMatNormMeta["PRF1", ]) # -0.7831032, P = [1] 1.798996e-207
hist(cpePredictionsMeta, col="#ff000066", breaks=20)
hist(cpePredictionsMetaFin, add=T, col="#00ff0066", breaks=20)

# Load the genotypes etc., this will allow us to order the expression and proportions data and create a snpsToGenesList.
load(file=theRootDir %&% "data/genotypes/metaBric/mataBricGenotypesTxt/metabricGenotypesNumFin.RData") # genotypeMatNumDat, genotypePcs, 

# save the final expression and proportions data with same column names in same order as the genotype data
names(cpePredictionsMetaFin) <- gsub(".", "-", names(cpePredictionsMetaFin), fixed=T)
cpePredictionsMetaFinOrd <- cpePredictionsMetaFin[colnames(genotypeMatNumDat)]
colnames(expressionMatNormMeta) <- gsub(".", "-", colnames(expressionMatNormMeta), fixed=T)
expressionMatNormMetaOrd <- expressionMatNormMeta[, colnames(genotypeMatNumDat)]
save(cpePredictionsMetaFinOrd, expressionMatNormMetaOrd, file=theRootDir %&% "data/genotypes/metaBric/metabricExpressionAndProps.RData")

# Subset the metabric Snps matrix to be a subset of the TCGA snps, as we only care about the overlap for this analysis....
genotypeMatNumDat_inTcga <- genotypeMatNumDat[rownames(genotypeMatNumDat) %in% rownames(genotypesMat_filt), ]
save(genotypePcs, genotypeMatNumDat_inTcga, file=theRootDir %&% "data/genotypes/metaBric/metabricExpressionAndProps_inTcga.RData")

## Map the expression to genotype data, creating snpsToGenesListMeta (snpsToGenesList) for 
snpsLocation_filt_metaBric <- snpsLocation_filt[snpsLocation_filt[, 1] %in% rownames(genotypeMatNumDat_inTcga), ] # keep the snp location for the relevant metabric snps.
geneLocation_filt_metaBric <- geneLocation_filt_qnDat[geneLocation_filt_qnDat[,1] %in% rownames(expressionMatNormMetaOrd), ] 

snpsLocation_filt_metaBric[,4] <- snpsLocation_filt_metaBric[,3] - 500000 # lets do a megabase window around the gene.
snpsLocation_filt_metaBric[,4][snpsLocation_filt_metaBric[,4] < 0] <- 0
snpsLocation_filt_metaBric[,5] <- snpsLocation_filt_metaBric[,3] + 500000
colnames(snpsLocation_filt_metaBric)[c(4,5)] <- c("minusMB", "plusMB")
rownames(snpsLocation_filt_metaBric) <- snpsLocation_filt_metaBric[,1]

# We will match snps to genes using the genomic ranges package in R, which will allow us to do this very quickly.
library(GenomicRanges) # this library should be able to do this fairly efficiently I think.
geneRanges <- GRanges(seqnames=Rle(geneLocation_filt_metaBric[,2]), ranges=IRanges(geneLocation_filt_metaBric[,3], geneLocation_filt_metaBric[,4]))
names(geneRanges) <- geneLocation_filt_metaBric[,1]
snpLocations <- GRanges(seqnames=Rle(snpsLocation_filt_metaBric[,2]), ranges=IRanges(snpsLocation_filt_metaBric[,4], snpsLocation_filt_metaBric[,5]))
names(snpLocations) <- snpsLocation_filt_metaBric[,1]

countOut <- findOverlaps(geneRanges,snpLocations, type="any", ignore.strand=T)
countOut_mat <- as.matrix(countOut)
snpSplit <- split(countOut_mat[,1], countOut_mat[,2])

snpsToGenesMat <- countOut_mat
snpsToGenesMat[,1] <- geneLocation_filt_metaBric[,1][countOut_mat[,1]]
snpsToGenesMat[,2] <- snpsLocation_filt_metaBric[,1][countOut_mat[,2]]
snpsToGenesListMeta <- split(snpsToGenesMat[,1], snpsToGenesMat[,2])
length(snpsToGenesListMeta)

save(snpsToGenesListMeta, file=theRootDir %&% "data/genotypes/metaBric/snpsToGenesListMeta.RData")

# Now estimate PEER factors (need the genotype data for this....). SWAP OUT THE TCGA INFO, PUT IN THE METABRIC INFO.....
# First, do this for a model that includes the CPE tumor purity estimates (we will use these PEER factors with the interaction model).
library("peer")
peerExprMat <- t(expressionMatNormMetaOrd)
model = PEER()
PEER_setPhenoMean(model,as.matrix(peerExprMat))
mod <- model.matrix(t(expressionMatNormMetaOrd[snpsToGenesListMeta[[1]], ,drop=FALSE])~(genotypeMatNumDat_inTcga[names(snpsToGenesListMeta)[1], ]*cpePredictionsMetaFinOrd)+genotypePcs$x[,1:3]) # This line is only used to make it easier to include co-variates below
PEER_setCovariates(model, mod[,3:6]) # include the purity estimate and the PCs as co-variates.
PEER_setNk(model,35)
PEER_update(model)
factors_correctedAnalysis = PEER_getX(model)
theFactors_correctedAnaysis_metabric <- factors_correctedAnalysis[,5:39] #' drop the real covariates (Pcs, purity), which are included in this object that peer returns.

# Now create the PEER factors for the conventional model.
modelPeerNoInt = PEER()
PEER_setPhenoMean(modelPeerNoInt,as.matrix(peerExprMat))
PEER_setCovariates(modelPeerNoInt, mod[,4:6]) # use only the PCs from the genotype matrix as the co-variates.
PEER_setNk(modelPeerNoInt,35)
PEER_update(modelPeerNoInt)
factors_NoInteractionAnalysis = PEER_getX(modelPeerNoInt)
theFactors_NoInteractionAnalysis_metabric <- factors_NoInteractionAnalysis[, 4:38]

save(theFactors_correctedAnaysis_metabric, theFactors_NoInteractionAnalysis_metabric, file=theRootDir %&% "data/genotypes/metaBric/metabricPeer.RData")

## PEER won't converge, no documentation to state why this would be the case, nor any info online.... I will manually estimate hidden factors with PCA.
# Instead, I will calculate the hidden factors manually using PCA on an expression matrix where the variation from genetic PCs and tumor puirt has been regressed out.
# This should be be qualitatively almost identical to PEER.
pcsFullExprMat <- prcomp(t(expressionMatNormMetaOrd))
pcMatFullExprMat <- pcsFullExprMat$x

for(i in 1:10)
{
   print(cor(pcMatFullExprMat[,i], cpePredictionsMetaFinOrd))
#   print(cor.test(pcMatFullExprMat[,i], genotypePcs$x[,1])$p.value)
}

# Expression matrix with variation due to ansectry regressed out.
exprMatNoGeno <- expressionMatNormMetaOrd
for(i in 1:nrow(expressionMatNormMetaOrd))
{
  exprMatNoGeno[i,] <- residuals(lm(expressionMatNormMetaOrd[i,]~genotypePcs$x[,1:3]))
  print(i)
}

# Expression matix with genotype and cell type proportion regressed out
exprMatNoGenoNoCpe <- expressionMatNormMetaOrd
for(i in 1:nrow(expressionMatNormMetaOrd))
{
  exprMatNoGenoNoCpe[i,] <- residuals(lm(expressionMatNormMetaOrd[i,]~genotypePcs$x[,1:3]+cpePredictionsMetaFinOrd))
  print(i)
}

## Now calculate the first 35 PCs on these two matrices and use these as our hidden factors.
thePcFactors_NoInteractionAnalysis_metabric <- prcomp(t(exprMatNoGeno))$x[, 1:35]
thePcFactors_correctedAnaysis_metabric <- prcomp(t(exprMatNoGenoNoCpe))$x[, 1:35]
save(thePcFactors_NoInteractionAnalysis_metabric, thePcFactors_correctedAnaysis_metabric, file=theRootDir %&% "data/genotypes/metaBric/metabricPcFactors.RData")



