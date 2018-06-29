#' This script process the VCF files that were returned by the Michigan imputation server and provides the imputed genotypes for the Breast Cancer GWAS significant hits. We are not permitted to re-post these germline genotype data as they are protected and authorization must be obtained from TCGA. 

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

#' It is assumed that the VCF files are here.
setwd(theRootDir %&% "vcfFiles/")

theFiles <- dir()
theSnpNames <- unlist(strsplit(theFiles[grep(".tfam", theFiles, fixed=T)], ".tfam"))
theFamFiles <- theFiles[grep(".tfam", theFiles, fixed=T)]
thePedFiles <- theFiles[grep(".tped", theFiles, fixed=T)]


#' Are the sample names in the fam files? Yes these are all identical and in identical order. Makes code below easier.
famFileSampNamesList <- read.delim(paste(theRootDir, "vcfFiles/", theSnpNames[i], ".tfam", sep=""), as.is=T, header=F, stringsAsFactors=FALSE, colClasses = c("character"))[,1]
numequal <- numeric()
for(i in 2:length(theSnpNames))
{
  test <- read.delim(paste(theRootDir, "vcfFiles/", theSnpNames[i], ".tfam", sep=""), as.is=T, header=F, stringsAsFactors=FALSE, colClasses = c("character"))[,1]
  numequal[i] <- sum(famFileSampNamesList == test)
}
unique(numequal)
# [1]   NA 8950

 
#' For each file, get the data in matrix format.
genotypesList <- list()
numGenotypesList <- list()
pedRowIndex <- rep(1, length(theSnpNames)) # For 3 snps, Rs1432679, rs2046210, rs7904519, their snp info is on the 2nd row of the ped file. This is because vcftools also gave back snps one base pair up/down stream of the snp of interest. I need to handle this here and make sure for these snps it pulls out the 2nd row
pedRowIndex[theSnpNames %in% c("rs1432679", "rs2046210", "rs7904519")] <- 2
hg19BrcaGwasHits <- read.delim(theRootDir %&% "brcaGwasSigSnps_hg19Locs_fromBiomart.txt", as.is=T)[-92, ] # the resulting tsv file (from biomart) is saved here, load it. The last row contains duplicated snp on a contig: drop it

for(i in 1:length(theSnpNames))
{
  thePedFile <- read.delim(paste(theRootDir, "vcfFiles/", theSnpNames[i], ".tped", sep=""), as.is=T, header=F, stringsAsFactors=FALSE, colClasses = c("character"))
  theFamFile <- read.delim(paste(theRootDir, "vcfFiles/", theSnpNames[i], ".tfam", sep=""), as.is=T, header=F, stringsAsFactors=FALSE, colClasses = c("character"))
  
  theCharData <- as.character(thePedFile[pedRowIndex[i],])
  allele1 <- theCharData[seq(5, length(theCharData)-1, 2)]
  allele2 <- theCharData[seq(6, length(theCharData), 2)]
  genotype <- paste(allele1, "/", allele2, sep="")
  names(genotype) <- theFamFile[,1]

  # find which are homozygotes, heterozygotes
  uniqueGenos <- unique(genotype)
  uniqueGenosMat <- do.call(rbind, strsplit(uniqueGenos, "/"))
  theHomos <- uniqueGenos[uniqueGenosMat[,1] == uniqueGenosMat[,2]]
  theHeteros <- uniqueGenos[uniqueGenosMat[,1] != uniqueGenosMat[,2]]
  genotype[genotype == theHeteros[2]] <- theHeteros[1] # make sure heterozygoes are only represented by 1 unique entry.
  genotypesList[[i]] <- genotype
  
  # and encode genotypes as 0,1, 2
  genotypes012 <- genotype
  genotypes012[genotype %in% theHomos[1]] <- 0
  genotypes012[genotype %in% theHomos[2]] <- 2
  genotypes012[genotype %in% theHeteros] <- 1
  genotypes012 <- as.numeric(genotypes012)

  numGenotypesList[[i]] <- genotypes012
  print(i)
}
names(numGenotypesList) <- theSnpNames
names(genotypesList) <- theSnpNames
numGenotypesMat <- do.call(rbind, numGenotypesList)
charGenotypesMat <- do.call(rbind, genotypesList)
rownames(charGenotypesMat) <- theSnpNames
colnames(charGenotypesMat) <- names(genotypesList[[1]])
rownames(numGenotypesMat) <- theSnpNames
colnames(numGenotypesMat) <- names(genotypesList[[1]])

#' Save the output, used by imputed_snps_with_gtex_analysis.R. Note that we are not permitted to re-post these data.
save(charGenotypesMat, numGenotypesMat, file="brcaGwasSigSnpsImputed.RData")














