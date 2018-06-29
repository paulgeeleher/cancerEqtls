#' I need to get the GTEx data in a format that I can use. I.e. the same format as the other TCGA results.

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

#' Read the annotations files and create a map of variant ID -> RS ID and ensemble ID to gene name.
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

#' I need to get a list of snps and genes to keep from my TCGA data.
load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheDatas.RData") # theData_conv_Peer, theData_int_selectPeers, theData_conv_noPeer, theData_int_noPeer, theData_int_Peer, theData_intInv_noPeer, theData_intInv_Peer, theData_intInv_selectPeers
# snpsAndGenesAllTCGA <- colnames(theData_conv_Peer[[1]])
rm(list=c("theData_int_selectPeers", "theData_conv_noPeer", "theData_int_noPeer", "theData_intInv_noPeer", "theData_intInv_selectPeers"))
load(file=theRootDir %&% "data/forMatrixEQTL/output/snpsAndGenesAllTCGA.RData") # snpsAndGenesAllTCGA; I saved this here.
snpsAndGenesMatTcga <- do.call(rbind, strsplit(snpsAndGenesAllTCGA, "~"))
allTcgaSnps <- unique(snpsAndGenesMatTcga[,2])
allTcgaGenes <- unique(snpsAndGenesMatTcga[,1])
length(allTcgaSnps) # [1] 615568
length(allTcgaGenes) # [1] 14995

#' This will be a lot faster if I filter to only the things I need before mapping the gene / variant IDs:
tcgaGenesEnsembl <- names(geneSyms[geneSyms %in% allTcgaGenes])
tcgaSnpsGtex <- names(geneEns[geneEns %in% allTcgaSnps])


#' Read the data for blood, breast, lcl and fibroblast. (lets skip whole blood for now)
wholeBloodResultsGtex <- read.delim(theRootDir %&% "data/Whole_Blood_Analysis.v6p.all_snpgene_pairs.txt", as.is=T)

# Keep the rows with genes and variants that are also in GTEx, and convert the gene and snps ids to symbols and rs ids.
gtexBrcaRowsInTcga <- intersect(which(wholeBloodResultsGtex$gene_id %in% tcgaGenesEnsembl), which(wholeBloodResultsGtex$variant_id %in% tcgaSnpsGtex))
wholeBloodResultsGtexInTcga <- wholeBloodResultsGtex[gtexBrcaRowsInTcga, ]
wholeBloodResultsGtexInTcga$geneIdSym <- geneSyms[wholeBloodResultsGtexInTcga$gene_id]
wholeBloodResultsGtexInTcga$rsId <- geneEns[wholeBloodResultsGtexInTcga$variant_id]
wholeBloodResultsGtexInTcga$geneToSnps <- paste(wholeBloodResultsGtexInTcga$geneIdSym, "~", wholeBloodResultsGtexInTcga$rsId, sep="")
wholeBloodResultsGtexInTcga_noDupEqtls <- wholeBloodResultsGtexInTcga[!duplicated(wholeBloodResultsGtexInTcga$geneToSnps), ]
rownames(wholeBloodResultsGtexInTcga_noDupEqtls) <- wholeBloodResultsGtexInTcga_noDupEqtls$geneToSnps

snpGenePairsInBothDatasets <- rownames(wholeBloodResultsGtexInTcga_noDupEqtls)[rownames(wholeBloodResultsGtexInTcga_noDupEqtls) %in% colnames(theData_int_Peer[[1]])]
ncol(theData_int_Peer[[1]]) # [1] 3602220
length(snpGenePairsInBothDatasets) # [1] 2992823
tcgaBetasNormal <- theData_int_Peer$betaMatAll[1, snpGenePairsInBothDatasets]
tcgaBetasCancer <- theData_intInv_Peer$betaMatAll[1, snpGenePairsInBothDatasets]
gtexBetasBlood <- wholeBloodResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasets, "slope"]
gtexPsBlood <- wholeBloodResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasets, "pval_nominal"]
names(gtexPsBlood) <- rownames(wholeBloodResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasets,])
names(gtexBetasBlood) <- rownames(wholeBloodResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasets,])

#' There is no evidence of similarity of the beta values (systematically). But the signal here could easily be obscured by noise. (NB these results are BS, alleles are flipped between studies, thus betas not comparable)
cor.test(tcgaBetasNormal, gtexBetasBlood, method="spearman") 
cor.test(tcgaBetasCancer, gtexBetasBlood, method="spearman") # 0.0003468816
cor.test(tcgaBetasNormal, tcgaBetasCancer, method="spearman") # -0.70565

#' The reference allele is flipped randomly for these...
cor.test(abs(tcgaBetasNormal), abs(gtexBetasBlood), method="spearman") # 0.08242306; p-value < 2.2e-16
cor.test(abs(tcgaBetasCancer), abs(gtexBetasBlood), method="spearman") # 0.09422852; p-value < 2.2e-16

#' Are the top hits in tcga "cancer" or "normal" enriched in the GTEx data. 
topHitsCancerTcga <- names(sort(theData_intInv_Peer$pMatAll[1, snpGenePairsInBothDatasets])[1:1000])
topHitsNormalTcga <- names(sort(theData_int_Peer$pMatAll[1, snpGenePairsInBothDatasets])[1:1000])

median(gtexPsBlood[topHitsNormalTcga]) # [1] 9.66238e-08
median(gtexPsBlood[topHitsCancerTcga]) # [1] 7.56049e-11
median(rank(gtexPsBlood)[topHitsNormalTcga]) # [1] 19944.5
median(rank(gtexPsBlood)[topHitsCancerTcga]) # [1] 10883.5


#' Load the GTEx LCL and fibroblast data. Get the P-values for all of these.
lclGtex <- read.delim(theRootDir %&% "data/Cells_EBV-transformed_lymphocytes_Analysis.v6p.all_snpgene_pairs.txt", as.is=T)

# Keep the rows with genes and variants that are also in TCGA, and convert the gene and snps ids to symbols and rs ids.
gtexLclRowsInTcga <- intersect(which(lclGtex$gene_id %in% tcgaGenesEnsembl), which(lclGtex$variant_id %in% tcgaSnpsGtex))
lclResultsGtexInTcga <- lclGtex[gtexLclRowsInTcga, ]
lclResultsGtexInTcga$geneIdSym <- geneSyms[lclResultsGtexInTcga$gene_id]
lclResultsGtexInTcga$rsId <- geneEns[lclResultsGtexInTcga$variant_id]
lclResultsGtexInTcga$geneToSnps <- paste(lclResultsGtexInTcga$geneIdSym, "~", lclResultsGtexInTcga$rsId, sep="")
lclResultsGtexInTcga_noDupEqtls <- lclResultsGtexInTcga[!duplicated(lclResultsGtexInTcga$geneToSnps), ]
rownames(lclResultsGtexInTcga_noDupEqtls) <- lclResultsGtexInTcga_noDupEqtls$geneToSnps

snpGenePairsInBothDatasetsLcl <- rownames(lclResultsGtexInTcga_noDupEqtls)[rownames(lclResultsGtexInTcga_noDupEqtls) %in% colnames(theData_int_Peer[[1]])]
ncol(theData_int_Peer[[1]]) # [1] 3602220
length(snpGenePairsInBothDatasetsLcl) # [1] 2992823
tcgaBetasNormalLcl <- theData_int_Peer$betaMatAll[1, snpGenePairsInBothDatasetsLcl]
tcgaBetasCancerLcl <- theData_intInv_Peer$betaMatAll[1, snpGenePairsInBothDatasetsLcl]
gtexBetasLcl <- lclResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsLcl, "slope"]
gtexPsLcl <- lclResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsLcl, "pval_nominal"]
names(gtexPsLcl) <- rownames(lclResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsLcl,])
names(gtexBetasLcl) <- rownames(lclResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsLcl,])

# There is no evidence of similarity of the beta values (systematically). But the signal here could easily be obscured by noise. (actually the allleles are flipped.)
cor.test(tcgaBetasNormalLcl, gtexBetasLcl, method="spearman") 
cor.test(tcgaBetasCancerLcl, gtexBetasLcl, method="spearman") # 0.0003468816
cor.test(tcgaBetasNormalLcl, tcgaBetasCancerLcl, method="spearman") # -0.70565


#' Load the fibroblast GTEx data.
fibroblastGtex <- read.delim(theRootDir %&% "data/Cells_Transformed_fibroblasts_Analysis.v6p.all_snpgene_pairs.txt", as.is=T)

gtexFibroblastRowsInTcga <- intersect(which(fibroblastGtex$gene_id %in% tcgaGenesEnsembl), which(fibroblastGtex$variant_id %in% tcgaSnpsGtex))
fibroblastResultsGtexInTcga <- fibroblastGtex[gtexFibroblastRowsInTcga, ]
fibroblastResultsGtexInTcga$geneIdSym <- geneSyms[fibroblastResultsGtexInTcga$gene_id]
fibroblastResultsGtexInTcga$rsId <- geneEns[fibroblastResultsGtexInTcga$variant_id]
fibroblastResultsGtexInTcga$geneToSnps <- paste(fibroblastResultsGtexInTcga$geneIdSym, "~", fibroblastResultsGtexInTcga$rsId, sep="")
fibroblastResultsGtexInTcga_noDupEqtls <- fibroblastResultsGtexInTcga[!duplicated(fibroblastResultsGtexInTcga$geneToSnps), ]
rownames(fibroblastResultsGtexInTcga_noDupEqtls) <- fibroblastResultsGtexInTcga_noDupEqtls$geneToSnps

snpGenePairsInBothDatasetsFibroblast <- rownames(fibroblastResultsGtexInTcga_noDupEqtls)[rownames(fibroblastResultsGtexInTcga_noDupEqtls) %in% colnames(theData_int_Peer[[1]])]
ncol(theData_int_Peer[[1]]) # [1] 3602220
length(snpGenePairsInBothDatasetsFibroblast) # [1] 2992823
tcgaBetasNormalFibroblast <- theData_int_Peer$betaMatAll[1, snpGenePairsInBothDatasetsFibroblast]
tcgaBetasCancerFibroblast <- theData_intInv_Peer$betaMatAll[1, snpGenePairsInBothDatasetsFibroblast]
gtexBetasFibroblast <- fibroblastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsFibroblast, "slope"]
gtexPsFibroblast <- fibroblastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsFibroblast, "pval_nominal"]
names(gtexPsFibroblast) <- rownames(fibroblastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsFibroblast,])
names(gtexBetasFibroblast) <- rownames(fibroblastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsFibroblast,])


#' Compare to breast data betas and p-values for TCGA and GTEx breast data.
breastResultsGtex <- read.delim(theRootDir %&% "data/Breast_Mammary_Tissue_Analysis.v6p.all_snpgene_pairs.txt", as.is=T)

gtexBrcaRowsInTcga <- intersect(which(breastResultsGtex$gene_id %in% tcgaGenesEnsembl), which(breastResultsGtex$variant_id %in% tcgaSnpsGtex))
breastResultsGtexInTcga <- breastResultsGtex[gtexBrcaRowsInTcga, ]
breastResultsGtexInTcga$geneIdSym <- geneSyms[breastResultsGtexInTcga$gene_id]
breastResultsGtexInTcga$rsId <- geneEns[breastResultsGtexInTcga$variant_id]
breastResultsGtexInTcga$geneToSnps <- paste(breastResultsGtexInTcga$geneIdSym, "~", breastResultsGtexInTcga$rsId, sep="")
breastResultsGtexInTcga_noDupEqtls <- breastResultsGtexInTcga[!duplicated(breastResultsGtexInTcga$geneToSnps), ]
rownames(breastResultsGtexInTcga_noDupEqtls) <- breastResultsGtexInTcga_noDupEqtls$geneToSnps

snpGenePairsInBothDatasetsBreast <- rownames(breastResultsGtexInTcga_noDupEqtls)[rownames(breastResultsGtexInTcga_noDupEqtls) %in% colnames(theData_int_Peer[[1]])]
ncol(theData_int_Peer[[1]]) # [1] 3602220
length(snpGenePairsInBothDatasetsBreast) #

#' get the "normal" tcga summary stats
tcgaBetasNormal <- theData_int_Peer$betaMatAll[1, snpGenePairsInBothDatasetsBreast]
tcgaPsNormal <- theData_int_Peer$pMatAll[1, snpGenePairsInBothDatasetsBreast]
tcgaSeNormal <- theData_int_Peer$errMatAll[1, snpGenePairsInBothDatasetsBreast]

#' get the "bulk" tcga summary stats
tcgaBetasBulk <- theData_conv_Peer$betaMatAll[1, snpGenePairsInBothDatasetsBreast]
tcgaSeBulk <- theData_conv_Peer$errMatAll[1, snpGenePairsInBothDatasetsBreast]
tcgaPsNormal <- theData_conv_Peer$pMatAll[1, snpGenePairsInBothDatasetsBreast]
tcgaSeNormal <- theData_conv_Peer$errMatAll[1, snpGenePairsInBothDatasetsBreast]

#' get the "cancer" tcga summary stats
tcgaBetasCancer <- theData_intInv_Peer$betaMatAll[1, snpGenePairsInBothDatasetsBreast]
tcgaSeCancer <- theData_intInv_Peer$errMatAll[1, snpGenePairsInBothDatasetsBreast]
tcgaPsCancer <- theData_intInv_Peer$pMatAll[1, snpGenePairsInBothDatasetsBreast]
gtexBetasBreast <- breastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsBreast, "slope"]
gtexPsBreast <- breastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsBreast, "pval_nominal"]
gtexSeBreast <- breastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsBreast, "slope_se"]
names(gtexPsBreast) <- rownames(breastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsBreast,])
names(gtexBetasBreast) <- rownames(breastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsBreast,])
names(gtexSeBreast) <- rownames(breastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsBreast,])
topPvalsGTExBreast <- sort(gtexPsBreast)[1:1000]
topHitsGTExBreast <- names(topPvalsGTExBreast)

#' There is no evidence of similarity of the beta values (systematically). But the signal here could easily be obscured by noise. (NB the beta values are flipped.)
cor.test(tcgaBetasNormal, gtexBetasBreast, method="spearman") 
cor.test(tcgaBetasCancer, gtexBetasBreast, method="spearman") # 
cor.test(tcgaBetasNormal, tcgaBetasCancer, method="spearman") # 

#' Actually the reference alleles have been flipped for a random bunch of these samples.Which ones are flipped? Flip them back (below). 
cor.test(abs(tcgaBetasCancer), abs(gtexBetasBreast), method="spearman") # 0.09315999 ; p-value < 2.2e-16
cor.test(abs(tcgaBetasNormal), abs(gtexBetasBreast), method="spearman") # 0.07742939 ; p-value < 2.2e-16


#' Are the top hits in tcga "cancer" or "normal" enriched in the GTEx data.
topHitsCancerTcga <- names(sort(theData_intInv_Peer$pMatAll[1, snpGenePairsInBothDatasetsBreast])[1:1000])
topHitsNormalTcga <- names(sort(theData_int_Peer$pMatAll[1, snpGenePairsInBothDatasetsBreast])[1:1000])
topHitsBulkTcga <- names(sort(theData_conv_Peer$pMatAll[1, snpGenePairsInBothDatasetsBreast])[1:1000])

middleHitsBulkTcga <- names(sort(theData_conv_Peer$pMatAll[1, snpGenePairsInBothDatasetsBreast])[40000:41000])

#' Rank-wise, this seems to be an improvement over the blood samples, despite smaller sample size for the breast dataset.
median(gtexPsBreast[topHitsNormalTcga]) # [1] [1] 5.940345e-06
median(gtexPsBreast[topHitsCancerTcga]) # [1] 1.311815e-11
median(rank(gtexPsBreast)[topHitsNormalTcga]) # [1] 13430
median(rank(gtexPsBreast)[topHitsCancerTcga]) # [1] 3131.5

cor.test(gtexBetasBreast[topHitsCancerTcga], tcgaBetasCancer[topHitsCancerTcga])

plot(gtexBetasBreast[topHitsCancerTcga], tcgaBetasCancer[topHitsCancerTcga])
plot(abs(gtexBetasBreast[topHitsCancerTcga]), abs(tcgaBetasCancer[topHitsCancerTcga]))
dev.off()

#' The directionalities have been flipped for many of these eQTLs. Reference Allele has changed!
cor.test(gtexBetasBreast[topHitsCancerTcga], tcgaBetasCancer[topHitsCancerTcga])
cor.test(abs(gtexBetasBreast[topHitsCancerTcga]), abs(tcgaBetasCancer[topHitsCancerTcga]))



#' NB: I want to flip for all the LCL, fibro, whole, breast. First, pull out the reference allele and snps names for all of these datasets, we will then combine these.

#' Alleles in the GTEx whole blood data. I guess the first allele is the reference allele. 
findAllels <- do.call(rbind, strsplit(wholeBloodResultsGtexInTcga_noDupEqtls[, "variant_id"], "_"))
refAllelGtex <- findAllels[,3]
names(refAllelGtex) <- wholeBloodResultsGtexInTcga_noDupEqtls[, "rsId"]
refAllelGtex_noDups <- refAllelGtex[!duplicated(names(refAllelGtex))]

#' Breast
findAllelsBreast <- do.call(rbind, strsplit(breastResultsGtexInTcga_noDupEqtls[, "variant_id"], "_"))
refAllelGtexBreast <- findAllelsBreast[,3]
names(refAllelGtexBreast) <- breastResultsGtexInTcga_noDupEqtls[, "rsId"]
refAllelGtexBreast_noDups <- refAllelGtexBreast[!duplicated(names(refAllelGtexBreast))]

#' LCL
findAllelsLcl <- do.call(rbind, strsplit(lclResultsGtexInTcga_noDupEqtls[, "variant_id"], "_"))
refAllelGtexLcl <- findAllelsLcl[,3]
names(refAllelGtexLcl) <- lclResultsGtexInTcga_noDupEqtls[, "rsId"]
refAllelGtexLcl_noDups <- refAllelGtexLcl[!duplicated(names(refAllelGtexLcl))]

#' Fibroblast 
findAllelsFibroblast <- do.call(rbind, strsplit(fibroblastResultsGtexInTcga_noDupEqtls[, "variant_id"], "_"))
refAllelGtexFibroblast <- findAllelsFibroblast[,3]
names(refAllelGtexFibroblast) <- fibroblastResultsGtexInTcga_noDupEqtls[, "rsId"]
refAllelGtexFibroblast_noDups <- refAllelGtexFibroblast[!duplicated(names(refAllelGtexFibroblast))]

#' Reference allele is the same for all common snps across datasets in GTEx, this is good.
commonSnpsAll <- Reduce(intersect, list(names(refAllelGtex_noDups), names(refAllelGtexBreast_noDups), names(refAllelGtexLcl_noDups), names(refAllelGtexFibroblast_noDups)))
length(commonSnpsAll) # [1] 569089
sum(refAllelGtex_noDups[commonSnpsAll] == refAllelGtexBreast_noDups[commonSnpsAll]) # [1] 569089
sum(refAllelGtex_noDups[commonSnpsAll] == refAllelGtexLcl_noDups[commonSnpsAll]) # [1] 569089
sum(refAllelGtex_noDups[commonSnpsAll] == refAllelGtexFibroblast_noDups[commonSnpsAll]) # [1] 569089


#' Combine blood, breast, lcl and fibroblast snps with reference allele. We will use TCGA below and figure out what needs to be flipped.
gtexSnpsWithRefAllele_all <- c(refAllelGtex_noDups, refAllelGtexBreast_noDups, refAllelGtexLcl_noDups, refAllelGtexFibroblast_noDups)
gtexSnpsWithRefAllele <- gtexSnpsWithRefAllele_all[!duplicated(names(gtexSnpsWithRefAllele_all))]

#' What has TCGA used as the reference allele? Load the data and compare to GTEx reference alleles above.
theMap <- read.delim(theRootDir %&% "GenomeWideSNP_6.rs_snp_map", as.is=T)
rsIds <- theMap[,1]
names(rsIds) <- theMap[,2]

#' This file contains the actual genotypes
theAlleles <- read.delim(theRootDir %&% "GenomeWideSNP_6_alleles.csv", as.is=T)

#' Write out a dataframe that has affy snp ids
alleleData <- cbind(theAlleles, rsIds[theAlleles[,1]])
alleleData_fin <- alleleData[-which(is.na(alleleData[,4])), c(4,2,3)]
rownames(alleleData_fin) <- alleleData_fin[,1]

#' Get the reference and alternate alleles for TCGA.
refAlleleTcga <- alleleData_fin[, "Allele.A"]
names(refAlleleTcga) <- alleleData_fin[, 1]
altAlleleTcga <- alleleData_fin[, "Allele.B"]
names(altAlleleTcga) <- alleleData_fin[, 1]

rsIdsAllelesInBoth <- names(gtexSnpsWithRefAllele)[names(gtexSnpsWithRefAllele) %in% names(refAlleleTcga)]
sum(gtexSnpsWithRefAllele[rsIdsAllelesInBoth] == refAlleleTcga[rsIdsAllelesInBoth]) # [1] 169541
sum(gtexSnpsWithRefAllele[rsIdsAllelesInBoth] != refAlleleTcga[rsIdsAllelesInBoth]) # [1] 417896
sum(gtexSnpsWithRefAllele[rsIdsAllelesInBoth] == altAlleleTcga[rsIdsAllelesInBoth]) # [1] 168975
sum(gtexSnpsWithRefAllele[rsIdsAllelesInBoth] != altAlleleTcga[rsIdsAllelesInBoth]) # [1] 418462
length(rsIdsAllelesInBoth) # [1] 587437 (all GTEx snps are in TCGA too.)

snpsWhereRefAlleleDoesntMatch <- rsIdsAllelesInBoth[gtexSnpsWithRefAllele[rsIdsAllelesInBoth] != refAlleleTcga[rsIdsAllelesInBoth]]
snpsWhereGtexRefMatchesTcgaAlt <- rsIdsAllelesInBoth[gtexSnpsWithRefAllele[rsIdsAllelesInBoth] == altAlleleTcga[rsIdsAllelesInBoth]]

#' Find alleles where the reference allele in GTEx matches neither the reference or alternate allele in TCGA.
gtexAlleleMathesNeither <- which((gtexSnpsWithRefAllele[rsIdsAllelesInBoth] != altAlleleTcga[rsIdsAllelesInBoth]) & (gtexSnpsWithRefAllele[rsIdsAllelesInBoth] != refAlleleTcga[rsIdsAllelesInBoth]))


#' Lets flip the beta value in GTEx if it has used a different reference allele to TCGA.
snpWithStrand <- read.csv(theRootDir %&% "affy6arrayStrandInfo.csv") # This file was part of the birdsuite metadata file (https://www.broadinstitute.org/ftp/pub/mpg/birdsuite/birdsuite_metadata_1.5.5.tgz)
strandVec <- as.character(snpWithStrand[, "strand"])
names(strandVec) <- snpWithStrand[, "dbsnp_rs_id"]
snpsOnReverseStrand <- snpWithStrand[, "dbsnp_rs_id"][snpWithStrand[, "strand"] == "-"]
snpsOnForwardStrand <- snpWithStrand[, "dbsnp_rs_id"][snpWithStrand[, "strand"] == "+"]

sum(snpsWhereRefAlleleDoesntMatch %in% snpsOnReverseStrand) # [1] 270737 : Presumably SNPs that don't match because of standedness? And possibly also reversal of reference allele? (need to flip the beta value for the ones that have had it flipped).
sum(snpsWhereRefAlleleDoesntMatch %in% snpsOnForwardStrand) # [1] 146474 : Presumably SNPs that don't match because of reversal of reference allele?
nonMatchingSnpsOnForwardStrand <- snpsWhereRefAlleleDoesntMatch[snpsWhereRefAlleleDoesntMatch %in% snpsOnForwardStrand]
sum(nonMatchingSnpsOnForwardStrand %in% snpsWhereGtexRefMatchesTcgaAlt) # [1] 146474 : Yep, in all cases these have had the reference allele flipped. Need to flip the beta value for ALL of these.
sum(!nonMatchingSnpsOnForwardStrand %in% snpsWhereGtexRefMatchesTcgaAlt) # [1] 0

#' Flip the SNPs on Allele A, where the snps was annotated on the reverse strand in TCGA
alleleA_forward <- as.character(alleleData_fin[, "Allele.A"])
snpsOnReverseStrand_index <- rownames(alleleData_fin) %in% snpsOnReverseStrand # index of snps on the reverse strand (that need to be flipped)
table(alleleA_forward)
alleleA_forward[alleleA_forward == "A" & snpsOnReverseStrand_index] <- "T" # flip all As to Ts (for snps on the reverse strand)
table(alleleA_forward)
alleleA_forward[alleleA_forward == "C" & snpsOnReverseStrand_index] <- "X" # flip all Cs to Xs (for snps on the reverse strand)
table(alleleA_forward)
alleleA_forward[alleleA_forward == "G" & snpsOnReverseStrand_index] <- "C" # flip all Gs to Cs (for snps on the reverse strand)
table(alleleA_forward)
alleleA_forward[alleleA_forward == "X" & snpsOnReverseStrand_index] <- "G" # flip all Xs to Gs (for snps on the reverse strand)
table(alleleA_forward)

#' Flip Allele B, where the SNP was annotated on the reverse strand in TCGA.
alleleB_forward <- as.character(alleleData_fin[, "Allele.B"])
table(alleleB_forward)
alleleB_forward[alleleB_forward == "T" & snpsOnReverseStrand_index] <- "A" # # flip all Ts to As (for snps on the reverse strand)
table(alleleB_forward)
alleleB_forward[alleleB_forward == "C" & snpsOnReverseStrand_index] <- "X" # Flip all Cs to Xs (for snps on the reverse strand)
table(alleleB_forward)
alleleB_forward[alleleB_forward == "G" & snpsOnReverseStrand_index] <- "C" # flip all Gs to Cs (for snps on the reverse strand)
table(alleleB_forward)
alleleB_forward[alleleB_forward == "X" & snpsOnReverseStrand_index] <- "G" # all flip Xs to Gs (for snps on the reverse strand)
table(alleleB_forward)

#' TCGA allele data WITH strand info and strand flipped snps so everything is forward strand.
strandInfoOrd <- strandVec[rownames(alleleData_fin)]
alleleData_fin_withFORWARD <- cbind(alleleData_fin, alleleA_forward, alleleB_forward, strandInfoOrd)

refAlleleTcga_forward <- alleleData_fin_withFORWARD[, "alleleA_forward"]
names(refAlleleTcga_forward) <- alleleData_fin_withFORWARD[, 1]
altAlleleTcga_forward <- alleleData_fin_withFORWARD[, "alleleB_forward"]
names(altAlleleTcga_forward) <- alleleData_fin_withFORWARD[, 1]

rsIdsAllelesInBoth <- names(gtexSnpsWithRefAllele)[names(gtexSnpsWithRefAllele) %in% names(refAlleleTcga_forward)]
sum(gtexSnpsWithRefAllele[rsIdsAllelesInBoth] == refAlleleTcga_forward[rsIdsAllelesInBoth]) # [1] 293250 
sum(gtexSnpsWithRefAllele[rsIdsAllelesInBoth] != refAlleleTcga_forward[rsIdsAllelesInBoth]) # [1] 294187
sum(gtexSnpsWithRefAllele[rsIdsAllelesInBoth] == altAlleleTcga_forward[rsIdsAllelesInBoth]) # [1] 293783
sum(gtexSnpsWithRefAllele[rsIdsAllelesInBoth] != altAlleleTcga_forward[rsIdsAllelesInBoth]) # [1] 293654

#' Now for TCGA Snps, where the reference allele still doesn't match GTEx, flip the AlleleA and AlleleB.
refAllelesNeedSwapping <- rsIdsAllelesInBoth[gtexSnpsWithRefAllele[rsIdsAllelesInBoth] != refAlleleTcga_forward[rsIdsAllelesInBoth]] ## NB: These should be the SNPs where the beta value needs to be flipped!!


#' Lets flip these and see if they now match, they should almost all match with the correct flipping.
alleleData_fin_withFORWARD_andSwaps <- cbind(alleleData_fin_withFORWARD, alleleA_forwardSwapped=alleleA_forward, alleleB_forwardSwapped=alleleB_forward)
for(i in 1:length(refAllelesNeedSwapping))
{
  newVals <- alleleData_fin_withFORWARD_andSwaps[refAllelesNeedSwapping[i], c("alleleB_forwardSwapped", "alleleA_forwardSwapped")] # b, then A
  alleleData_fin_withFORWARD_andSwaps[refAllelesNeedSwapping[i], c("alleleA_forwardSwapped", "alleleB_forwardSwapped")] <- newVals
  print(i)
}


#' Flip the GTEx BREAST betas, where needed. 
reOrder_breastStuff <- breastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsBreast, "rsId"]
flipIndex <- which(reOrder_breastStuff %in% refAllelesNeedSwapping)
gtexBetasBreast_flipped <- gtexBetasBreast
gtexBetasBreast_flipped[flipIndex] <- (gtexBetasBreast[flipIndex] * -1)

#' Flip LCLs.
reOrder_lclStuff <- lclResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsLcl, "rsId"]
flipIndexLcl <- which(reOrder_lclStuff %in% refAllelesNeedSwapping)
gtexBetasLcl_flipped <- gtexBetasLcl
gtexBetasLcl_flipped[flipIndexLcl] <- (gtexBetasLcl[flipIndexLcl] * -1)

#' Flip Fibroblasts
reOrder_fibroblastStuff <- fibroblastResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasetsFibroblast, "rsId"]
flipIndexFibroblast <- which(reOrder_fibroblastStuff %in% refAllelesNeedSwapping)
gtexBetasFibroblast_flipped <- gtexBetasFibroblast
gtexBetasFibroblast_flipped[flipIndexFibroblast] <- (gtexBetasFibroblast[flipIndexFibroblast] * -1)

#' Flip whole blood.
reOrder_wholeBloodStuff <- wholeBloodResultsGtexInTcga_noDupEqtls[snpGenePairsInBothDatasets, "rsId"]
flipIndex <- which(reOrder_wholeBloodStuff %in% refAllelesNeedSwapping)
gtexBetasBlood_flipped <- gtexBetasBlood
gtexBetasBlood_flipped[flipIndex] <- (gtexBetasBlood[flipIndex] * -1)


#' Now that the alleles are flipped, these are very correlated.
cor.test(gtexBetasBreast_flipped[topHitsCancerTcga], tcgaBetasCancer[topHitsCancerTcga]) #
# data:  gtexBetasBreast_flipped[topHitsCancerTcga] and tcgaBetasCancer[topHitsCancerTcga]
# t = 57.991, df = 998, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8631471 0.8916057
# sample estimates:
#      cor 
# 0.878151 

plot(gtexBetasBreast_flipped[topHitsCancerTcga], tcgaBetasCancer[topHitsCancerTcga], col="#00000055")
cor.test(gtexBetasBreast_flipped[topHitsCancerTcga], tcgaBetasBulk[topHitsCancerTcga])
cor.test(gtexBetasBreast_flipped[topHitsBulkTcga], tcgaBetasBulk[topHitsBulkTcga])
cor.test(gtexBetasBreast_flipped[middleHitsBulkTcga], tcgaBetasBulk[middleHitsBulkTcga])
cor.test(gtexBetasBreast_flipped[middleHitsBulkTcga], tcgaBetasNormal[middleHitsBulkTcga])
cor.test(gtexBetasBreast_flipped[middleHitsBulkTcga], tcgaBetasCancer[middleHitsBulkTcga], col="#00000055")
cor.test(gtexBetasBreast_flipped[topHitsCancerTcga], tcgaBetasNormal[topHitsCancerTcga])
cor.test(gtexBetasLcl_flipped[topHitsCancerTcga], tcgaBetasCancer[topHitsCancerTcga])
cor.test(gtexBetasLcl_flipped[middleHitsBulkTcga], tcgaBetasBulk[middleHitsBulkTcga])
cor.test(gtexBetasFibroblast_flipped[topHitsCancerTcga], tcgaBetasCancer[topHitsCancerTcga], method="spearman")
cor.test(gtexBetasFibroblast_flipped[topHitsCancerTcga], tcgaBetasBulk[topHitsCancerTcga], method="spearman")
cor.test(gtexBetasFibroblast_flipped[middleHitsBulkTcga], tcgaBetasBulk[middleHitsBulkTcga])


#' Also load the top hits for random CPE estimates
load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheDatas_theData_intInv_Peer_RandomCpe.RData") # theData_intInv_Peer_RandomCpe
topHitsRandom <- names(sort(theData_intInv_Peer_RandomCpe$pMatAll[1, ])[1:8833])


#' TCGA sig hits.
sigHitsCancer <- names(which(p.adjust(theData_intInv_Peer$pMatAll[1, ], method="BH") < 0.05))
length(sigHitsCancer) # [1] 8833
topHitsCancer <- names(sort(theData_intInv_Peer$pMatAll[1, ])[1:57189])
sigHitsBulk <- names(which(p.adjust(theData_conv_Peer$pMatAll[1, ], method="BH") < 0.05))
length(sigHitsBulk) # [1] 57189
topHitsBulk <- sort(theData_conv_Peer$pMatAll[1, ])[1:8833]
sigHitsBulkBetas <- theData_intInv_Peer$betaMatAll[1, sigHitsBulk]
topHitsBulk <- sort(p.adjust(theData_conv_Peer$pMatAll[1, ], method="BH"))[1:8833]
sum(sigHitsCancer %in% names(topHitsBulk)) # [1] 5279
lostEqtls <- sigHitsBulk[!sigHitsBulk %in% sigHitsCancer]
length(lostEqtls) # [1] 49647

gtexSigEqtlsBreast <- rownames(breastResultsGtexInTcga_noDupEqtls)[p.adjust(breastResultsGtexInTcga_noDupEqtls[, "pval_nominal"], method="BH") < 0.05]
length(gtexSigEqtlsBreast) # [1] 27683                                                                                                                                                                                                      
nominallySigBreast <- rownames(breastResultsGtexInTcga_noDupEqtls)[breastResultsGtexInTcga_noDupEqtls[, "pval_nominal"] < 0.05]
length(nominallySigBreast)
gtexSigEqtlsLcl <- rownames(lclResultsGtexInTcga_noDupEqtls)[p.adjust(lclResultsGtexInTcga_noDupEqtls[, "pval_nominal"], method="BH") < 0.05]
length(gtexSigEqtlsLcl) # [1] 16932 
gtexSigEqtlsFibroblast <- rownames(fibroblastResultsGtexInTcga_noDupEqtls)[p.adjust(fibroblastResultsGtexInTcga_noDupEqtls[, "pval_nominal"], method="BH") < 0.05]
length(gtexSigEqtlsFibroblast) # [1] 86070
gtexSigEqtlsWholeBlood <- rownames(wholeBloodResultsGtexInTcga_noDupEqtls)[p.adjust(wholeBloodResultsGtexInTcga_noDupEqtls[, "pval_nominal"], method="BH") < 0.05]
length(gtexSigEqtlsWholeBlood)

#' How many bulk eQTLs are called in either tissue type, versus how many could have been called...
sum(sigHitsBulk %in% gtexSigEqtlsBreast) # [1] 12153, breast eQTLs sig in TCGA bulk.
sum(sigHitsBulk %in% rownames(breastResultsGtexInTcga_noDupEqtls)) # [1] 51137, breast eQTLs could have been called in TCGA bulk.
sum(sigHitsBulk %in% gtexSigEqtlsLcl) # [1] 6847
sum(sigHitsBulk %in% rownames(lclResultsGtexInTcga_noDupEqtls)) # [1] 44094
sum(sigHitsBulk %in% gtexSigEqtlsFibroblast) # [1] 20491
sum(sigHitsBulk %in% rownames(fibroblastResultsGtexInTcga_noDupEqtls)) # [1] 45680
sum(sigHitsBulk %in% gtexSigEqtlsWholeBlood) # [1] 14862
sum(sigHitsBulk %in% rownames(wholeBloodResultsGtexInTcga_noDupEqtls)) # [1] 46532

#' How many lost eQTLs are 
sum(lostEqtls %in% gtexSigEqtlsBreast) # [1] 8769
sum(lostEqtls %in% gtexSigEqtlsLcl) # [1] 4719
sum(lostEqtls %in% gtexSigEqtlsFibroblast) # [1] 16484

#' For how many of these is directionality the same as for the conventional model?
lostAndSigBreast <- lostEqtls[lostEqtls %in% gtexSigEqtlsBreast]
lostAndSigLcl <- lostEqtls[lostEqtls %in% gtexSigEqtlsLcl]
lostAndSigFibro <- lostEqtls[lostEqtls %in% gtexSigEqtlsFibroblast]
sameBreast <- gtexBetasBreast_flipped[lostAndSigBreast] * sigHitsBulkBetas[lostAndSigBreast]
sum(sameBreast > 0) # [1] 8536
sum(sameBreast < 0) # [1] 233
sameBreastConsistent <- sameBreast[sameBreast > 0]
sameBreastNotConsistent <- sameBreast[sameBreast < 0]
sameLcl <- gtexBetasLcl_flipped[lostAndSigLcl] * sigHitsBulkBetas[lostAndSigLcl]
sum(sameLcl > 0) # [1] 4531
sum(sameLcl < 0) # [1] 188
sameLclConsistent <- sameLcl[sameLcl > 0]
sameLclNotConsistent <- sameLcl[sameLcl < 0]
sameFibro <- gtexBetasFibroblast_flipped[lostAndSigFibro] * sigHitsBulkBetas[lostAndSigFibro]
sum(sameFibro > 0) # [1] 15810
sum(sameFibro < 0) # [1] 674
sameFibroConsistent <- sameFibro[sameFibro > 0]
sameFibroNotConsistent <- sameFibro[sameFibro < 0]


#' How many of these "lost" eQTLs are significant with concordant directionality in at least one of breast, lcl or fibroblast.
allConsistenteQTLs <- union(union(names(sameBreastConsistent), names(sameLclConsistent)), names(sameFibroConsistent))
length(allConsistenteQTLs) # [1] 18595

allNotConsistenteQTLs <- union(union(names(sameBreastNotConsistent), names(sameLclNotConsistent)), names(sameFibroNotConsistent))
length(allNotConsistenteQTLs) # [1] 897

#' How many of these significant bulk eQTLs are significant with concordant directionality in breast, lcl or fibroblast!?
sum(sigHitsCancer %in% gtexSigEqtlsBreast) # [1] 3397
sum(sigHitsCancer %in% rownames(breastResultsGtexInTcga_noDupEqtls)) # [1] 7795
sum(sigHitsCancer %in% gtexSigEqtlsLcl) # [1] 2158
sum(sigHitsCancer %in% rownames(lclResultsGtexInTcga_noDupEqtls)) # [1] 6858
sum(sigHitsCancer %in% gtexSigEqtlsFibroblast) # [1] 4064
sum(sigHitsCancer %in% rownames(fibroblastResultsGtexInTcga_noDupEqtls)) # [1] 7085
sum(sigHitsCancer %in% gtexSigEqtlsWholeBlood) # [1] [1] 3267
sum(sigHitsCancer %in% rownames(wholeBloodResultsGtexInTcga_noDupEqtls)) # [1] 7108

sigInLclsAndFibroblasts <- union(gtexSigEqtlsLcl, gtexSigEqtlsFibroblast)
sigInLclFibro_butNotBreast <- sigInLclsAndFibroblasts[!sigInLclsAndFibroblasts %in% nominallySigBreast]
sigInBreast_butNotLclAndFibro <- nominallySigBreast[!nominallySigBreast %in% sigInLclsAndFibroblasts]
sigInBreastLclAndFibro <- Reduce(intersect, list(gtexSigEqtlsLcl, gtexSigEqtlsFibroblast, gtexSigEqtlsBreast))
sum(sigHitsBulk %in% sigInLclsAndFibroblasts)

length(sigInLclFibro_butNotBreast) # [1] 51641
sum(colnames(theData_intInv_Peer$pMatAll) %in% sigInLclFibro_butNotBreast) # [1] 47196
a <- colnames(theData_intInv_Peer$pMatAll)[colnames(theData_intInv_Peer$pMatAll) %in% sigInLclFibro_butNotBreast]

sum(sigHitsCancer %in% sigInLclFibro_butNotBreast) # [1] 572
sum(sigHitsBulk %in% sigInLclFibro_butNotBreast) # [1] 5440

sum(sigHitsCancer %in% sigInBreastLclAndFibro) # [1] 572
sum(sigHitsBulk %in% sigInBreastLclAndFibro) # [1] 5440

fisher.test(matrix(c(5440, 51749, 572, 8261), nrow=2))
fisher.test(matrix(c(572, 8261, 5440, 51749), nrow=2))$p.value # [1] 8.065267e-22

#' Make Figure 3(b)
barplotData <- c(572/8261, 5440/51749)
names(barplotData) <- c("Cancer", "Bulk Tumor")
dir.create(theRootDir %&% "paper/figures/figure3/", recursive=T, showWarnings=F)
svg(theRootDir %&% "paper/figures/figure3/figBarplot.svg", width=3, height=3)
barplot(barplotData, las=1, cex.axis=.8, col=c("#1f78b4", "#1a9641"), ylab="Proportion of likely false positives", bty="l")
dev.off()


##' BELOW HERE IS SOME CODE TO COMPARE THE EQtl profiles in GTEx and TCGA, i.e. what is the same, what is different, genome wide....
# This Z score test is proposed in the following stackexchange threads  and several papers:
# https://stats.stackexchange.com/questions/93540/testing-equality-of-coefficients-from-two-different-regressions
# https://stats.stackexchange.com/questions/55501/test-a-significant-difference-between-two-slope-values
# Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859–866.
# Clogg, C. C., Petkova, E., & Haritou, A. (1995). Statistical methods for comparing regression coefficients between models. American Journal of Sociology, 100(5), 1261-1293.
zScore <- numeric(length(tcgaBetasCancer))
for(i in 1:length(tcgaBetasCancer)) # length: [1] 3270829
{
  zScore[i] <- (tcgaBetasCancer[i] - gtexBetasBreast_flipped[i]) / sqrt(tcgaSeCancer[i]^2 + gtexSeBreast[i]^2)
  print(i)
}
names(zScore) <- names(tcgaBetasCancer)


#' The vector of p-values for differences between GTEx breast and TCGA breast *cancer* cells.
pValues = (2*pnorm(-abs(zScore)))
names(pValues) <- names(tcgaBetasCancer)
hist(pValues, breaks=100)
dev.off()
pSort <- sort(pValues)
pAdj <- p.adjust(pSort, method="BH")
pSort[1:10]
sum(pAdj < 0.05) # [1] 3885
changingEqtls <- names(which(pAdj < 0.05))

#' This is a Supplementary Figure
svg(theRootDir %&% "paper/figures/histForDifferences_tcgaCancer_and_gtexBreast.svg", width=4, height=4)
hist(pValues, main="", cex.axis=.8, col="#8dd3c7", las=1, xlab="P-Value")
dev.off()


#' A supplementary Table, with the eQTLs that are different between GTEx and TCGA. 
# Format: eQTL, ZScore, Pvalue difference, FDR difference, PValue Cancer, Beta Cancer, SE Cancer, Pvalue GTEx, Beta GTEx, SE Gtex.
suppTabsDiffs <- data.frame(changingEqtls, pValues[changingEqtls], pAdj[changingEqtls], tcgaPsCancer[changingEqtls], tcgaBetasCancer[changingEqtls], tcgaSeCancer[changingEqtls], gtexPsBreast[changingEqtls], gtexBetasBreast_flipped[changingEqtls], gtexSeBreast[changingEqtls])
colnames(suppTabsDiffs) <- c("eQTL", "P-value (change)", "FDR (change)", "P-value (cancer)", "Effect size (cancer)", "Standard Error (cancer)", "P-value (GTEx breast)", "Effect Size (GTEx breast)", "Standard Error (GTEx breast)")
write.csv(suppTabsDiffs, file=theRootDir %&% "/paper/suppTab_gtex_and_cancer_diffs.csv", row.names=F)

#' How many of these "changing eQTLs" are significant in TCGA and how many in GTEx? Do we tend to lose / gain eQTLs in TCGA?
tcgaFdrsCancer <- p.adjust(tcgaPsCancer, method="BH")
gtexFdrsBreast <- p.adjust(gtexPsBreast, method="BH")
sum(tcgaFdrsCancer < 0.05) # [1] 7710
sum(gtexFdrsBreast < 0.05) # [1] 31547

sigEqtlsTcgaCancer <- names(which(tcgaFdrsCancer < 0.05))
sigEqtlsGtexBreast <- names(which(gtexFdrsBreast < 0.05))

sum(gtexFdrsBreast[changingEqtls] < 0.05) # [1] 3040
sum(tcgaFdrsCancer[changingEqtls] < 0.05) # [1] 628

changingEqtls_sigInTcgaCancer <- changingEqtls[changingEqtls %in% sigEqtlsTcgaCancer]
changingEqtls_sigInGtexBreast <- changingEqtls[changingEqtls %in% sigEqtlsGtexBreast]
changingEqtls_sigInBoth <- changingEqtls_sigInTcgaCancer[changingEqtls_sigInTcgaCancer %in% changingEqtls_sigInGtexBreast]
length(changingEqtls_sigInBoth) # [1] 205

eQtlsLostInCancer <- changingEqtls_sigInGtexBreast[!changingEqtls_sigInGtexBreast %in% changingEqtls_sigInBoth] # eQTLs significant in GTEx, but not in TCGA cancer, and with statistical evidecen of a difference.
eQtlsGainedInCancer <- changingEqtls_sigInTcgaCancer[!changingEqtls_sigInTcgaCancer %in% changingEqtls_sigInBoth] # eQTLs significant in Cancer, but not in GTEx, and with statistical evidence of a difference.


#' What are these genes and their characteristics, figure out wiht GOseq
pAdjSortAll <- p.adjust(pSort, method="BH")
sum(pAdjSortAll < 0.05) # [1] 3885
sigDiffEqtls <- names(pAdjSortAll[pAdjSortAll < 0.05])
sigDiffEqtls_greaterCancer <- sigDiffEqtls[abs(tcgaBetasCancer[sigDiffEqtls]) > abs(gtexBetasBreast_flipped[sigDiffEqtls])]
length(sigDiffEqtls_greaterCancer) # [1] 797
sigDiffEqtls_greaterNormal <- sigDiffEqtls[abs(tcgaBetasCancer[sigDiffEqtls]) < abs(gtexBetasBreast_flipped[sigDiffEqtls])]
length(sigDiffEqtls_greaterNormal) # [1] 3068
genesMoreEQtlCancer <- unique(do.call(cbind, strsplit(sigDiffEqtls_greaterCancer, "~"))[1,]) # Genes with greater eQTL in cancer.
genesMoreEQtlNormal <- unique(do.call(cbind, strsplit(sigDiffEqtls_greaterNormal, "~"))[1,]) # Genes with greater eQTL in normal.
sum(genesMoreEQtlCancer %in% genesMoreEQtlNormal) # [1] 45; There's some overlap here....

write.table(genesMoreEQtlNormal, quote=F, row.names=F)
write.table(genesMoreEQtlCancer, quote=F, row.names=F)

#' Sig genes for CANCER profile
sigGenes <- unique(do.call(cbind, strsplit(sigDiffEqtls, "~"))[1,])

#' Sig genes for BULK profile.
pAdjSortAll_bulk <- p.adjust(pSort_bulk, method="BH")
sum(pAdjSortAll_bulk < 0.05) # [1] 3885
sigDiffEqtls_bulk <- names(pAdjSortAll_bulk[pAdjSortAll_bulk < 0.05])
sigGenes_bulk <- unique(do.call(cbind, strsplit(sigDiffEqtls_bulk, "~"))[1,])


#' Use GOseq to adjust for the number of SNPs targetting each gene.
genesToSnpsList <- do.call(rbind, strsplit(names(tcgaBetasCancer), "~"))
genesToSnpsList <- split(genesToSnpsList[,2], genesToSnpsList[,1])
genesToSnpsNum <- sapply(genesToSnpsList, length)
hist(genesToSnpsNum)



#' The code to run GOseq will look something like this (where "deVec" is a vector of 1s and 0s indicating which genes are DM and "numProbesVec" is the number of probes associated with each gene's promoter region):
library("goseq")

#' Stronger eQTL in cancer (Supp fig)
numProbesVec <- genesToSnpsNum
deVec <- rep(0, length(genesToSnpsNum))
names(deVec) <- names(genesToSnpsNum)
deVec[genesMoreEQtlCancer] <- 1
wf = nullp(deVec, "hg18", "geneSymbol", bias.data=numProbesVec)
GOout = goseq(wf, "hg18", "geneSymbol", test.cats=c("GO:BP"))
library("GO.db")
goterms <- Term(GOTERM)
GOout_new <- cbind(GOout, goterms[GOout[,1]])
GOout_new_filt <- GOout_new[GOout_new[, "numInCat"] > 15, ]
GOout_new_filt[1:10, ]

library(Hmisc)
dim(GOout_new_filt) # [1] 3679    6
barplotData <- rev(-log10(GOout_new_filt[1:8, "over_represented_pvalue"]))
names(barplotData) <- rev(capitalize(as.character(GOout_new_filt[1:8, 6])))
svg(theRootDir %&% "paper/figures/figGoTerms_moreEqtlInCancer.svg", width=6, height=3.5)
par(mar=c(5.1, 17, 4.1, 2.1))
barplot(barplotData, las=1, cex.axis=.8, col=c("#1f78b4"), xlab="-log10(P-Value)", bty="l", horiz=T)
dev.off()

#' Normal (Supp Fig.)
deVecNormal <- rep(0, length(genesToSnpsNum))
names(deVecNormal) <- names(genesToSnpsNum)
deVecNormal[genesMoreEQtlNormal] <- 1
wf = nullp(deVecNormal, "hg18", "geneSymbol", bias.data=numProbesVec)
GOout = goseq(wf, "hg18", "geneSymbol", test.cats=c("GO:BP"))
goterms <- Term(GOTERM)
GOout_new <- cbind(GOout, goterms[GOout[,1]])
GOout_new_filt <- GOout_new[GOout_new[, "numInCat"] > 15, ]
GOout_new_filt[1:10, ]

barplotData <- rev(-log10(GOout_new_filt[1:8, "over_represented_pvalue"]))
names(barplotData) <- rev(capitalize(as.character(GOout_new_filt[1:8, 6])))
svg(theRootDir %&% "paper/figures/figGoTerms_moreEqtlInNormal.svg", width=7, height=3.5)
par(mar=c(5.1, 23, 4.1, 1.1))
barplot(barplotData, las=1, cex.axis=.8, col=c("#1f78b4"), xlab="-log10(P-Value)", bty="l", horiz=T)
dev.off()


#' Everything combined (Fig. 3(c))
deVecAll <- rep(0, length(genesToSnpsNum))
names(deVecAll) <- names(genesToSnpsNum)
deVecAll[sigGenes] <- 1
wf = nullp(deVecAll, "hg18", "geneSymbol", bias.data=numProbesVec)
GOout = goseq(wf, "hg18", "geneSymbol", test.cats=c("GO:BP"))
library("GO.db")
goterms <- Term(GOTERM)
GOout_new <- cbind(GOout, goterms[GOout[,1]])
GOout_new_filt <- GOout_new[GOout_new[, "numInCat"] > 15, ]
GOout_new_filt[1:10,]


#' Export this to a supplementary Table:
colnames(GOout_new_filt)[6] <- "GO term name"
write.csv(GOout_new_filt, file=theRootDir %&% "paper/GO_results_all.csv", row.names=F)

#' Make a plot outta the top 8.
library(Hmisc)
dim(GOout_new_filt) # [1] 3679    6
barplotData <- rev(-log10(GOout_new_filt[1:8, "over_represented_pvalue"]))
names(barplotData) <- rev(capitalize(as.character(GOout_new_filt[1:8, 6])))
dir.create()
svg(theRootDir %&% "paper/figures/figure3/figGoTerms.svg", width=3, height=3.5)
barplot(barplotData, las=1, cex.axis=.8, col=c("#1f78b4"), xlab="-log10(P-Value)", bty="l", horiz=T)
dev.off()


print(sessionInfo())







