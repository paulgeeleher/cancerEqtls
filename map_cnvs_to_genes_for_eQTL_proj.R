#' An R script to map CNV regions to genes."

#' Intersect the gene with the overlapping CNV level. If the gene isn't fully captured by the CNV, an NA will be assigned.

#' Load libraries
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

#' Set the root directory. 
# theRootDir <- "/home/pgeeleher/postdoc_stuff/drug_repurposing/data/tcga_cnv_subtracted/" # the location to where the raw data was downloaded by the "download_tcga_data.R" script.

#' Set root directiory for this project.
source("theRootDir.R")

dir.create(paste(theRootDir, "/dataIn/tcga_cnv_subtracted/cnvsMappedToGenes/", sep=""), showWarnings = FALSE)


#' Load the gene ranges for HG19 using
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(GenomicFeatures)
geneRanges <- genes(txdb)
library(org.Hs.eg.db)
e2s = toTable(org.Hs.egSYMBOL)
syms <- e2s[, "symbol"]
names(syms) <- e2s[, "gene_id"]
theGeneSymsOrd <- syms[as.character(geneRanges$gene_id)]
# load(file="/home/pgeeleher/postdoc_stuff/drug_repurposing/data/geneRangesHg19.RData") # geneRanges

#' For each cancer type, load the CNV data and map to genes.
diseaseAbbrvs <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "STES", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
for(j in 1:length(diseaseAbbrvs))
{
  fname <- paste(theRootDir, "dataIn/tcga_cnv_subtracted/gdac.broadinstitute.org_", diseaseAbbrvs[j],".Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2015082100.0.0/", diseaseAbbrvs[j],".snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", sep="")
  Cnvs <- read.table(fname, header=TRUE)
  CnvsList <- split(Cnvs, Cnvs[, "Sample"])

  cnvEdGenesList <- list()
  ampGenesList <- list()
  delGenesList <- list()
  numGenesQuantifid <- numeric()
  theCnvQuantVecList <- list()
  for(i in 1:length(CnvsList))
  {
    chrs <- paste("chr", CnvsList[[i]]$Chromosome, sep="")
    starts <- CnvsList[[i]]$Start
    ends <- CnvsList[[i]]$End
    grCnvs <- GRanges(seqnames=Rle(chrs),ranges=IRanges(starts, ends), segMeans=CnvsList[[i]]$Segment_Mean)
    
    # Amp or del
    grCnvs_ampDel <- grCnvs[grCnvs$segMeans > 1 | grCnvs$segMeans < -1]
    cnvedGenes <- subsetByOverlaps(geneRanges, grCnvs_ampDel, type="within")
    cnvEdGenesList[[i]] <- cnvedGenes$gene_sym
    
    # Amps
    grCnvs_amp <- grCnvs[grCnvs$segMeans > 1]
    ampedGenes <- subsetByOverlaps(geneRanges, grCnvs_amp, type="within")
    ampGenesList[[i]] <- ampedGenes$gene_sym
    
    # Dels
    grCnvs_Del <- grCnvs[grCnvs$segMeans < -1]
    deledGenes <- subsetByOverlaps(geneRanges, grCnvs_Del, type="within")
    delGenesList[[i]] <- deledGenes$gene_sym
    
    # Continuous gene level 
    # Use count overlaps to find genes that unambiguously overlap a single peak. Give it an NA it it doesn't overlap a single peak. Assign it the value of the peak if it unambiguously overlaps a peak. PC.
    numOverlaps <- countOverlaps(geneRanges, grCnvs)
    numGenesQuantifid[i] <- sum(numOverlaps == 1)
    inCnv <- which(numOverlaps == 1) # take only gene unambiguously overlaping a peak, this is usually most genes.
    
    theCnvQuantVec <- rep(NA, length(geneRanges))
    olaps <- findOverlaps(geneRanges, grCnvs, type="within")
    theCnvQuantVec[olaps@queryHits] <- grCnvs$segMeans[olaps@subjectHits]
    theCnvQuantVecList[[i]] <- theCnvQuantVec
    # names(theCnvQuantVecList[[i]]) <- geneRanges$gene_sym
    names(theCnvQuantVecList[[i]]) <- theGeneSymsOrd
    print(i)
  }
  
  names(theCnvQuantVecList) <- names(CnvsList)
  theCnvQuantVecList_mat <- do.call(cbind, theCnvQuantVecList)
  siteVec <- sapply(strsplit(names(CnvsList), "-"), function(l)return(l[4]))
  tumorSamps <- which(siteVec == "01A")
  save(theCnvQuantVecList_mat, tumorSamps, file=paste(theRootDir, "dataIn/tcga_cnv_subtracted/cnvsMappedToGenes/", diseaseAbbrvs[j], ".RData", sep="")) # Save these RData files for use by other scripts.
  print(j)
}

print(sessionInfo())