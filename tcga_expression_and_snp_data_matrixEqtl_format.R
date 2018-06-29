#' Load and subset and order the 3 datasets (expression, genotypes, co-variate). I need to match the patient IDs and order and subset the samples.

#' Set the root directory where the data will be stored. NB: this directory needs to be set / created based on your own system!! The file "theRootDir.R" is sourced by most of these scripts and should include the desired root directory used by your entire analysis.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

#' Load the genotype data, I want to keep the calls for the germline genotypes.
load(file=theRootDir %&% "data/genotypes/geneotypeMatrices/callsMatrix.RData") # callsMatrix

## This code would keep tumor, i.e. 01A samples. I do not want to do this.
# callsMatrix_tumor <- callsMatrix[, do.call(rbind, strsplit(colnames(callsMatrix), "-", fixed=T))[,4] == "01A"] # Keep only unique tumor samples.

callsMatrix_tumor <- callsMatrix[, do.call(rbind, strsplit(colnames(callsMatrix), "-", fixed=T))[,4] == "10A"] ## NB 10A are GERMLINE samples (from blood). 01A are Tumor!
callsMatrix_tumor_nodups <- callsMatrix_tumor[, !duplicated(substring(colnames(callsMatrix_tumor), 9, 12))]
colnames(callsMatrix_tumor_nodups) <- substring(colnames(callsMatrix_tumor_nodups), 9, 12)


#' Load the breast cancer gene expression data and put on patient IDs as column names
load(file=theRootDir %&% "dataIn/tenRuvNewStandardApproach.RData") # cancerTypesVec, tenRuvNewStandardApproach
brcaExpression <- tenRuvNewStandardApproach[, cancerTypesVec == "BRCA"]
print(table(do.call(rbind, strsplit(colnames(brcaExpression), ".", fixed=T))[,4]))
brcaExpression_tumor <- brcaExpression[, do.call(rbind, strsplit(colnames(brcaExpression), ".", fixed=T))[,4] == "01A"] # Keep only unique tumor samples.
colnames(brcaExpression_tumor) <- substring(colnames(brcaExpression_tumor), 9, 12) # make the column names be only the patient ID.


#' Load my "co-variates", i.e. the heterogeneity estimates for each sample. Lets use the Cell Reports estimates to begin with.
tcgaBreasDeconvData <- read.csv(theRootDir %&% "EstimatedCellTypeProportionsTCGA.csv", as.is=T)
cellRepsTumorProps <- apply(data.matrix(tcgaBreasDeconvData[, 2:6]), 1, sum)
tCellProps <- tcgaBreasDeconvData[, "Immune"]
stromalProps <- tcgaBreasDeconvData[, "Stromal"]
normalProps <- tcgaBreasDeconvData[, "Normal.epithelial"]
names(cellRepsTumorProps) <- tcgaBreasDeconvData[,1]
names(tCellProps) <- tcgaBreasDeconvData[,1]
names(stromalProps) <- tcgaBreasDeconvData[,1]
names(normalProps) <- tcgaBreasDeconvData[,1]
cellRepsTumorProps_tumor <- cellRepsTumorProps[substring(names(cellRepsTumorProps), 14, 15) == "01"]
tCellProps_tumor <- tCellProps[substring(names(tCellProps), 14, 15) == "01"]
stromalProps_tumor <- stromalProps[substring(names(stromalProps), 14, 15) == "01"]
normalProps_tumor <- normalProps[substring(names(stromalProps), 14, 15) == "01"]
cellRepsTumorProps_tumor_noDups <- cellRepsTumorProps_tumor[!duplicated(substring(names(cellRepsTumorProps_tumor), 9, 12))]
tCellProps_tumor_noDups <- tCellProps_tumor[!duplicated(substring(names(tCellProps_tumor), 9, 12))]
stromalProps_tumor_nodups <- stromalProps_tumor[!duplicated(substring(names(stromalProps_tumor), 9, 12))]
normalProps_tumor_nodups <- normalProps_tumor[!duplicated(substring(names(normalProps_tumor), 9, 12))]
names(cellRepsTumorProps_tumor_noDups) <- substring(names(cellRepsTumorProps_tumor_noDups), 9, 12)
names(tCellProps_tumor_noDups) <- substring(names(tCellProps_tumor_noDups), 9, 12)
names(stromalProps_tumor_nodups) <- substring(names(stromalProps_tumor_nodups), 9, 12)
names(normalProps_tumor_nodups) <- substring(names(normalProps_tumor_nodups), 9, 12)
# save(tCellProps_tumor_noDups, stromalProps_tumor_nodups, normalProps_tumor_nodups, file="/mnt/data_scratch/prediXcanProj/data/tCellProportions_cellReps.RData")

# stromal and immue signals seem reasonably independent.... (both also highly correlated with total number of cancer cells, cor = 0.65 and -.67)
cor.test(stromalProps_tumor_nodups, tCellProps_tumor_noDups)
# t = 2.6115, df = 1043, p-value = 0.009145
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.02005368 0.14055475
# sample estimates:
#        cor 
# 0.08059871




commonPatientsAll <- intersect(colnames(callsMatrix_tumor_nodups), intersect(names(cellRepsTumorProps_tumor_noDups), colnames(brcaExpression_tumor)))

#' Subset and order the three datasets.
callsMatrix_tumor_nodups_ord <- callsMatrix_tumor_nodups[-1, commonPatientsAll]
brcaExpression_tumor_ord <- brcaExpression_tumor[, commonPatientsAll]
cellRepsTumorProps_tumor_noDups_ord <- cellRepsTumorProps_tumor_noDups[commonPatientsAll]

#' Write out the subsetted and ordered data to tab-delimited text files for matrixEQTL.
#' Genotype calls
dir.create(theRootDir %&% "data/forMatrixEQTL")
callsMatrix_tumor_nodups_ord <- cbind(rownames(callsMatrix_tumor_nodups_ord), callsMatrix_tumor_nodups_ord)
callsMatrix_tumor_nodups_ord[1,] <- colnames(callsMatrix_tumor_nodups_ord)
callsMatrix_tumor_nodups_ord[1,1] <- "id"
write.table(callsMatrix_tumor_nodups_ord, file=theRootDir %&% "data/forMatrixEQTL/callsMatrix.txt", sep="\t", row.names=F, col.names=F, quote=F)
rsIds <- callsMatrix_tumor_nodups_ord[,1]
rsLocations <- read.table(theRootDir %&% "GenomeWideSNP_6.hg18.map", sep="\t", header=F, as.is=T) # this file is provided with Birdsuite and can be dowloaded from "https://www.broadinstitute.org/birdsuite/birdsuite-downloads", in the zip file "birdsuite_metadata_1.5.5.tgz". It contains mappings of the RSids to HG18 genonomic locations for the birdsuite calls. We will use this with matrixEQTL.
rownames(rsLocations) <- rsLocations[,2]
length(rsIds) # [1] 906599
sum(rsIds %in% rownames(rsLocations)) # [1] 906598: It looks like we have location info for all of the rsIds. So we can just order this and write it out.
rsLocationsOrd <- rsLocations[rsIds[2:length(rsIds)], c(2, 1, 4)] # Order the rows and columns as necessary (to match the rsId matrix and the correct format)
rsLocationsOrd_fin <- rbind(c("snpid", "chr", "pos"), rsLocationsOrd)
write.table(rsLocationsOrd_fin, file=theRootDir %&% "data/forMatrixEQTL/snpLocMatrix.txt", sep="\t", row.names=F, col.names=F, quote=F)

#' Expression Data
brcaExpression_tumor_ord <- cbind(rownames(brcaExpression_tumor_ord), brcaExpression_tumor_ord)
geneNamesAll <- brcaExpression_tumor_ord[,1]
write.table(geneNamesAll[-1], file=theRootDir %&% "data/forMatrixEQTL/genesIds.txt", sep="\t", row.names=F, col.names=F, quote=F) # Get the start and end location of each gene from biomaRt hg18 (http://may2009.archive.ensembl.org/biomart/martview/). (this is ncbi36 / hg18)
geneLocs <- read.table(file=theRootDir %&% "data/forMatrixEQTL/biomart_geneLocs_hg18.txt", sep="\t", header=T, as.is=T)
geneLocsNoDups <- geneLocs[!duplicated(geneLocs[,3]), ]
geneLocsNoDups_theChrs <- geneLocsNoDups[geneLocsNoDups[, 2] %in% c(as.character(1:22), "X", "Y"), ] # Only keep genes on the main chromosomes (ignore the various contigs)
rownames(geneLocsNoDups_theChrs) <- geneLocsNoDups_theChrs[,3]
genesWithLocs <- rownames(brcaExpression_tumor_ord)[rownames(brcaExpression_tumor_ord) %in% rownames(geneLocsNoDups_theChrs)] # the genes for which we have location info. The gene expression data and the gene locations files may have to be identical and ordered the same, I don't know, but I might as well make them this way to be safe, cause that's how they are in the matrixEqtl documentation.
geneLocsNoDups_formatted_ord <- geneLocsNoDups_theChrs[genesWithLocs, c(3,2,1,4)] # data is now organised to be in same order as gene expression data.
geneLocsNoDups_formatted_fin <- rbind(c("geneid", "chr", "left", "right"), geneLocsNoDups_formatted_ord)
brcaExpression_tumor_ord_ord <- brcaExpression_tumor_ord[genesWithLocs, ] # this data is now ordered and subsetted to exactly match the snp data (for samples) and the gene location data (for gene ids)
brcaExpression_tumor_ord_fin <- rbind(colnames(brcaExpression_tumor_ord_ord), brcaExpression_tumor_ord_ord)
brcaExpression_tumor_ord_fin[1,1] <- "id"
write.table(brcaExpression_tumor_ord_fin, file=theRootDir %&% "data/forMatrixEQTL/brcaExpressionMatrix.txt", sep="\t", row.names=F, col.names=F, quote=F)
write.table(geneLocsNoDups_formatted_fin, file=theRootDir %&% "data/forMatrixEQTL/brcaGeneLocationInfo.txt", sep="\t", row.names=F, col.names=F, quote=F)

#' Proportions
cellRepsTumorProps_tumor_noDups_mat <- rbind(names(cellRepsTumorProps_tumor_noDups_ord), cellRepsTumorProps_tumor_noDups_ord)
cellRepsTumorProps_tumor_noDups_mat <- cbind(c("id", "cancerProp"), cellRepsTumorProps_tumor_noDups_mat)
write.table(cellRepsTumorProps_tumor_noDups_mat, file=theRootDir %&% "data/forMatrixEQTL/coVariates.txt", sep="\t", row.names=F, col.names=F, quote=F)



print(sessionInfo())