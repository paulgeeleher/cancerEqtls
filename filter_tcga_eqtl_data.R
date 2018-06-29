#' In the preliminary analysis I concluded that the following additional filtering steps were required, which largely follow the guildelines of GTEx:

#'  I should scale the expression data to a standard normal distribution, as per GTEx. This means a linear model can be used.
#'  I should remove genes not expressed in 75% of samples. 
#'  I should remove men (if there are any, there are).  
#'  I should remove SNPs where there is a low call rate (-1 in my data).
#'  I should remove SNPs on the Y Chromosome.
#'  I should correct for ancestry (GTEx used 3 PCs...)
#'  GTEx also used PEER on the expression data. 

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

#' Load the other relevant files that I have created. I.e. the SNP files, (created in previous script "tcga_expression_and_snp_data_matrixEqtl_format.R")
genotypesMat <- read.delim(theRootDir %&% "data/forMatrixEQTL/callsMatrix.txt", as.is=T)
rownames(genotypesMat) <- genotypesMat[,1]
genotypesMat[,1] <- NULL
genotypesMat <- data.matrix(genotypesMat)
expressionMat <- read.delim(theRootDir %&% "data/forMatrixEQTL/brcaExpressionMatrix.txt", as.is=T)
expressionMat <- expressionMat[-which(duplicated(expressionMat[,1])),]
rownames(expressionMat) <- expressionMat[,1]
expressionMat[,1] <- NULL
expressionMat <- data.matrix(expressionMat)
covariatesMat <- read.delim(theRootDir %&% "data/forMatrixEQTL/coVariates.txt", as.is=T)
rownames(covariatesMat) <- covariatesMat[,1]
covariatesMat[,1] <- NULL
covariatesMat <- data.matrix(covariatesMat)
covariatesVec <- as.numeric(covariatesMat)
names(covariatesVec) <- colnames(covariatesMat)

#' Load the locations of the snps and genes. I want to boot things that are on the Y chromsome.
geneLocation <- read.delim(theRootDir %&% "data/forMatrixEQTL/brcaGeneLocationInfo.txt", as.is=T)
geneLocation <- geneLocation[-which(duplicated(geneLocation[,1])),]
rownames(geneLocation) <- geneLocation[,1]

# identify (for removal) snps on Y or "XY"
snpsLocation <- read.delim(theRootDir %&% "data/forMatrixEQTL/snpLocMatrix.txt", as.is=T)
print(table(snpsLocation[,2])) # what the hell is the XY chromsome!?
snpsOnYChrom <- snpsLocation[,"snpid"][snpsLocation %in% c("Y", "XY")]

# identify (for removal) genes on Y
print(table(geneLocation[,2]))
genesOnYChrom <- geneLocation[,1][geneLocation[,2] %in% "Y"]


#' Transform the expression matrix to a matrix where the gene data is normally distributed.
expressionMatNorm <- expressionMat
for(i in 1:nrow(expressionMat))
{
  expressionMatNorm[i,] <- qnorm((rank(expressionMat[i,]) / (length(expressionMat[i,])+1)))
  print(i)
}

a <- unique(as.numeric(genotypesMat))
print(a) # [1]  2  1  0 -1


#' Load the original RAW breast cancer expression data and find if there are genes not expressed (i.e. don't have mapped reads) in many (>90%) samples. Remove these genes from further analysis.
brcaDataLoc <- paste(theRootDir %&% "dataIn/rnaSeq/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2015082100.0.0/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt", sep="")
tpmDatMat_bc <- read.delim(brcaDataLoc, as.is=T)
tpmDatMat_bc_tpm <- apply(tpmDatMat_bc[-1,which(tpmDatMat_bc[1,] == "scaled_estimate")], 2, as.numeric)
tpmDatMat_bc_tpm <- tpmDatMat_bc[-1,which(tpmDatMat_bc[1,] == "scaled_estimate")]
tpmDatMat_bc_tpm <- apply(tpmDatMat_bc_tpm, 2, as.numeric)
geneNames <- do.call(cbind, strsplit(tpmDatMat_bc[, "Hybridization.REF"], "|", fixed=TRUE))[1,][-1]
rownames(tpmDatMat_bc_tpm) <- geneNames
colnames(tpmDatMat_bc_tpm) <- substr(colnames(tpmDatMat_bc_tpm), 1, 28)
tpmDatMat_bc_tpm_logged <- log((tpmDatMat_bc_tpm*1000000)+1)
tpmDatMat_bc_tpm_logged_tumor <- tpmDatMat_bc_tpm_logged[, do.call(rbind, strsplit(colnames(tpmDatMat_bc_tpm_logged), ".", fixed=T))[,4] == "01A"]
numZeros <- apply(tpmDatMat_bc_tpm_logged_tumor, 1, function(vec)sum(vec == 0))
library(preprocessCore)
tpmDatMat_bc_tpm_logged_tumor_quantNorm <- normalize.quantiles(tpmDatMat_bc_tpm_logged_tumor)


#' This seems somewhat bimodal with peaks at 0% and 100%
hist(numZeros, breaks=100)
dev.off()

#' Identify (for removal) genes not expressed in > 90% of samples.
nonExpressedGenes <- rownames(tpmDatMat_bc_tpm_logged_tumor)[which(numZeros / ncol(tpmDatMat_bc_tpm_logged_tumor) > .9)]


#' Calculate the number of missing values for each genotype, keep track of Snps with an NA (i.e.  a -1 in the matrix)
numNas <- numeric(nrow(genotypesMat))
for(i in 1:nrow(genotypesMat))
{
  numNas[i] <- sum(genotypesMat[i,] == -1)
  # print(i)
}
snpsWithNas <- rownames(genotypesMat)[which(numNas != 0)]

#' Calculate the MAF for each genotype. This is the number of occuracnes of the less frequent allele divided by total number of alleles.
theMaf <- numeric(nrow(genotypesMat))
numSamps <- ncol(genotypesMat)
for(i in 1:nrow(genotypesMat))
{
  tab <- table(factor(genotypesMat[i,], levels=c(0,1,2))) # Count the number of 0, 1, 2 in each row of the genotype matrix.
  f1 <- ((tab["0"] * 2) + tab["1"]) / (sum(tab)*2) # count the total number of one of the alleles and divide by the total number alleles (total number of (0, 1, 2) x 2), to get the allele frequency for that allele
  f2 <- 1 - f1
  theMaf[i] <- min(c(f1, f2))
  print(i)
}

hist(theMaf, breaks=100)
dev.off()
print(sum(theMaf < 0.05)) # [1] 166150
hasLowMaf <- rownames(genotypesMat)[theMaf < 0.05]

#' Are there males in the dataset? There are.... I need to boot them too (this is breast cancer).
clinicalDataLocation <- paste(theRootDir %&% "dataIn/clinical/nationwidechildrens.org_clinical_patient_brca.txt", sep="")
clinDataBrca <- read.delim(clinicalDataLocation, as.is=T)[-(1:2),]
tcgaPatientIdsClin <- substring(clinDataBrca[, "bcr_patient_barcode"], 9, 12)
tcgaGenderClin <- clinDataBrca[, "gender"]
malePateints <- tcgaPatientIdsClin[tcgaGenderClin == "MALE"]
femalePateints <- tcgaPatientIdsClin[tcgaGenderClin == "FEMALE"]
theRace <- clinDataBrca[, "race"]
names(theRace) <- tcgaPatientIdsClin

print(sum(colnames(genotypesMat) %in% malePateints)) # [1] 11 Number with wrong genomtypes: [1] 12

#' I need to correct for Race, i.e. calculate PCs of the genotype matrix, plot these. 

#' filter the SNP data.
snpsToRemove <- unique(c(snpsWithNas, hasLowMaf, snpsOnYChrom))
individualsToRemove <- malePateints
genotypesMat_filt <- genotypesMat[which(!rownames(genotypesMat) %in% snpsToRemove), which(!colnames(genotypesMat) %in% individualsToRemove)]

#' Calculate the principal components of these SNP data. Plot it and color by race/ethnicity.
genotypePcs <- prcomp(t(genotypesMat_filt))

#' As expected, the samples heavily cluster by race.
theRace_ordFilt <- theRace[colnames(genotypesMat_filt)]
theRace_ordFilt[theRace_ordFilt == "ASIAN"] <- "#0000FF66" # blue
theRace_ordFilt[theRace_ordFilt == "AMERICAN INDIAN OR ALASKA NATIVE"] <- "#00FFFF66" # yellow
theRace_ordFilt[theRace_ordFilt == "BLACK OR AFRICAN AMERICAN"] <- "#00FF0066" # green
theRace_ordFilt[theRace_ordFilt == "[Not Evaluated]"] <- "#FFFF0066" # Purple
theRace_ordFilt[theRace_ordFilt == "[Not Available]"] <- "#FFFF0066" # Purple
theRace_ordFilt[theRace_ordFilt == "WHITE"] <- "#FF000066" # red
plot(genotypePcs$x[,1:2], col=theRace_ordFilt)


#' Filter the gene expression data. I need to remove these things from both the "normalized" and the "regular" gene exrpession data, and the male samples.
genesToRemove <- unique(c(genesOnYChrom, nonExpressedGenes)) # a list of genes to remove
expressionMatNorm_filt <- expressionMatNorm[which(!rownames(expressionMatNorm) %in% genesToRemove), which(!colnames(expressionMatNorm) %in% individualsToRemove)]
expressionMat_filt <- expressionMat[which(!rownames(expressionMat) %in% genesToRemove), which(!colnames(expressionMatNorm) %in% individualsToRemove)]


#' Calculate PCs of the gene expression matrix. Are these correlated with purity? 
expressionNormPcs <- prcomp(t(expressionMatNorm_filt))
expressionPcs <- prcomp(t(expressionMat_filt))

#' plot eigenvalues
plot((expressionNormPcs$sd^2) / sum(expressionNormPcs$sd^2))

#' Remove the males from the covariates matrix too.
covariatesVec_filt <- covariatesVec[which(!colnames(covariatesMat) %in% individualsToRemove)]

pcCorPs <- numeric()
pcCorRawExpr <- numeric()
for(i in 1:ncol(expressionNormPcs$x))
{
  pcCorPs[i] <- cor.test(expressionNormPcs$x[,i], covariatesVec_filt, method="pearson")$p.value
  pcCorRawExpr[i] <- cor.test(expressionPcs$x[,i], covariatesVec_filt, method="pearson")$p.value
  print(i)
}
print(pcCorPs[1:10]) # all of the first 10 PCs are correlated with the number of tumor cells.
plot(expressionNormPcs$x[,1], covariatesVec_filt)

#' Filter and order these files too, I will use them in downsteam analysis.
snpsLocation_filt <- snpsLocation[which(!rownames(genotypesMat) %in% snpsToRemove), ]
geneLocation_filt <- geneLocation[which(!rownames(expressionMat) %in% genesToRemove), ]
sum(snpsLocation_filt[,1] != rownames(genotypesMat_filt)) # [1] 0
sum(geneLocation_filt[,1] != rownames(expressionMat_filt)) # [1] 0


#' I need to save out these files and the eQTL analysis with the filtered, normalized, co-variates included data.
dir.create(theRootDir %&% "Results/rDatas/", recursive=TRUE, showWarnings = FALSE)
save(covariatesVec_filt, expressionMatNorm_filt, expressionMat_filt, genotypesMat_filt, genotypePcs, snpsLocation_filt, geneLocation_filt, file=theRootDir %&% "Results/rDatas/gwasFilteredInputData.RData")


#' 1. I also want to use the CPE estimates from the Nature Communications data as the primary method of estimating purity. Because its a consensus method. I will need to load these and impute the missing values using a principal component regression (there are only 7 missing values).
#' 2. For expression, I want to use the raw, quantiles normalized, filtered for non-expressed genes, mapped to normal distribution data.

#' Load the Nature Communications data and subset and order the same as the Cell Reports data.
nComsProps <- read.csv(theRootDir %&% "ncomms9971-s2.csv", as.is=T)
nComsProps_brca <- nComsProps[nComsProps[,2] == "BRCA", ]
nComsProps_brca_noCpe <- nComsProps_brca[!is.na(nComsProps_brca[,"CPE"]), ] # remove things where CPE wasn't calculated.
nComsProps_brca_noCpe_noDups <- nComsProps_brca_noCpe[(!duplicated(substring(nComsProps_brca_noCpe[,1], 9,12))), ]  # remove duplicates
nComsProps_brca_noCpe_noDups_onlyPrimary01 <- nComsProps_brca_noCpe_noDups[substring(nComsProps_brca_noCpe_noDups[,1], 14, 15) == "01", ]
theCpeProps <- nComsProps_brca_noCpe_noDups_onlyPrimary01[, "CPE"]
names(theCpeProps) <- substring(nComsProps_brca_noCpe_noDups_onlyPrimary01[, "Sample.ID"], 9, 12)

theCpeProps_matched <- theCpeProps[names(covariatesVec_filt)]
names(theCpeProps_matched) <- names(covariatesVec_filt)


#' Impute the missing values by leveraging the fact that PC1 of the gene expression matrix is highly correlated with purity estimates. There are only 7 missing values.
cor.test(expressionNormPcs$x[,1], theCpeProps_matched, method="pearson") # p e-118, cor -.63
missingValsCpe <- is.na(theCpeProps_matched)
theCpeProps_matched[missingValsCpe] <- mean(na.omit(theCpeProps_matched))
thepredpc <- expressionNormPcs$x[!missingValsCpe,1]
theMod <- lm(theCpeProps_matched[!missingValsCpe]~thepredpc)
plot(expressionNormPcs$x[!missingValsCpe,1], theCpeProps_matched[!missingValsCpe])
abline(coef=coef(lm(theCpeProps_matched[!missingValsCpe]~expressionNormPcs$x[!missingValsCpe,1])))
df <- data.frame(expressionNormPcs$x[missingValsCpe,1])
colnames(df) <- "thepredpc"
imputedCpes <- predict(theMod, df)
theCpeProps_matched[names(imputedCpes)] <- imputedCpes
cor.test(theCpeProps_matched, covariatesVec_filt) # 0.6186267 


#' Use the raw breast cancer expression data: 
#' 1. Quantile normalize this data. 
#' 2. Remove genes not expressed in 75% of samples. I want to be more strict on this than GTEx, because I want to avoid as much as I can introducing any sort of bias with this interaction term.
#' 3. Then map to normal distribution.

#' First, subset the raw data to keep only the patients for which we have genotype and purity estimates.
colnames(tpmDatMat_bc_tpm_logged_tumor) <- substring(colnames(tpmDatMat_bc_tpm_logged_tumor), 9, 12)
tpmDatMat_bc_tpm_logged_tumor_keepPatients <- tpmDatMat_bc_tpm_logged_tumor[, colnames(expressionMatNorm_filt)]

tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn <- normalize.quantiles(tpmDatMat_bc_tpm_logged_tumor_keepPatients)
rownames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn) <- rownames(tpmDatMat_bc_tpm_logged_tumor_keepPatients)
colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn) <- colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients)

#' Now normalize each gene to a normal distribtion, as was done by GTEx.
tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm <- tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn
for(i in 1:nrow(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn))
{
  tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm[i,] <- qnorm((rank(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn[i,]) / (length(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn[i,])+1)))
  print(i)
}

numNonExpressed <- apply(tpmDatMat_bc_tpm_logged_tumor_keepPatients, 1, function(vec)sum(vec == 0))
sum(numNonExpressed == 0) # [1] 12825; there are a lot of ubiquitously expressed genes in these samples....
sum(numNonExpressed == 1016) # [1] 298; a small number of genes never expressed.
genesNotExpressedEnough <- rownames(tpmDatMat_bc_tpm_logged_tumor)[which(numNonExpressed / ncol(tpmDatMat_bc_tpm_logged_tumor) > .25)] # genes should be expressed (i.e. have mapped reads) in 75% of samples. Removes 3654 genes.

genesToFilt <- unique(c(genesOnYChrom, genesNotExpressedEnough)) # a list of genes to remove
tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal <- tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm[which(!rownames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm) %in% genesToFilt), ]

#' test are these expression estimates still correlated with tumor purity? 
pcTest <- prcomp(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal))$x
cor.test(pcTest[,1], theCpeProps_matched)

#' It seems that the scaling differences between the samples (as measured by median expression) are correlated with tumor "purity"...!!!? 
sampMeds <- apply(tpmDatMat_bc_tpm_logged_tumor_keepPatients, 2, median)
cor.test(sampMeds, theCpeProps_matched)

#' Is purity correlated with the number of expressed genes!? It should be If there are MORE cell types.
numExpressedGenesPerSample <- apply(tpmDatMat_bc_tpm_logged_tumor_keepPatients, 2, function(vec)sum(vec != 0))
cor.test(numExpressedGenesPerSample, theCpeProps_matched)

plot(theCpeProps_matched, numExpressedGenesPerSample)

#' How many genes are correlated with tumor purity in the quantile normalized, versus non-quantile normalized data....
thePsNoQn <- numeric()
for(i in 1:nrow(tpmDatMat_bc_tpm_logged_tumor_keepPatients))
{
  thePsNoQn[i] <- cor.test(tpmDatMat_bc_tpm_logged_tumor_keepPatients[i,], theCpeProps_matched)$p.value
}
sum(na.omit(thePsNoQn) < 0.05)

thePsWithQn <- numeric()
for(i in 1:nrow(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn))
{
  thePsWithQn[i] <- cor.test(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn[i,], theCpeProps_matched)$p.value
}
sum(na.omit(thePsWithQn) < 0.05)


thePsWithQn_fin <- numeric()
theCorsWithQn_fin <- numeric()
for(i in 1:nrow(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal))
{
  thePsWithQn_fin[i] <- cor.test(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal[i,], theCpeProps_matched)$p.value
  theCorsWithQn_fin[i] <- cor.test(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal[i,], theCpeProps_matched)$estimate
}
sum(na.omit(thePsWithQn_fin) < 0.05)
hist(theCorsWithQn_fin, breaks=100)
dev.off()

thePsWithQn_orig <- numeric()
theCorsWithQn_orig <- numeric()
for(i in 1:nrow(expressionMatNorm_filt))
{
  thePsWithQn_orig[i] <- cor.test(expressionMatNorm_filt[i,], theCpeProps_matched)$p.value
  theCorsWithQn_orig[i] <- cor.test(expressionMatNorm_filt[i,], theCpeProps_matched)$estimate
}
sum(na.omit(thePsWithQn_orig) < 0.05)
mean(theCorsWithQn_orig) # [1] -0.08910712
hist(theCorsWithQn_orig, breaks=100) # distribtions of these correlation look very similar in quantiles and non-quantiles normalized data. I should stick to using the quantiles normalized data.
dev.off()

commonGenes <- rownames(geneLocation)[rownames(geneLocation) %in% rownames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal)]
tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin <- tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal[commonGenes, ]
geneLocation_filt_qnDat <- geneLocation[rownames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin), ]

#' Save all the same things as before, but inlcude the quantiles normalized expression data (with stronger filtering) and the CPE Nature Communications purity estimtes.
save(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin, geneLocation_filt_qnDat, theCpeProps_matched, covariatesVec_filt, expressionMatNorm_filt, expressionMat_filt, genotypesMat_filt, genotypePcs, snpsLocation_filt, geneLocation_filt, file=theRootDir %&% "Results/rDatas/gwasFilteredInputData_WithQNexprAndCpe.RData")




print(sessionInfo())