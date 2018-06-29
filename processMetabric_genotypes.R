# Load the metabric genotypes files and convert to a matrix of 0,1,2, with the heterozygotes encoded as 1, so I can analyse these data with a linear model. This file will also remove SNPs with low call rate and/or low minor allele frequency.


#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

#' The genotypes are here:
genoFiles <- dir(theRootDir %&% "data/genotypes/metaBric/mataBricGenotypesTxt/")

a <- read.delim(theRootDir %&% "data/genotypes/metaBric/mataBricGenotypesTxt/test.txt", comment.char="#")

genotypeList <- list()
for(i in 1:length(genoFiles)) # Read the genotype calls (the confident calls, not the forced calls, I should get rid of SNPs with a "low" (< 99%?) call rate?)
{
  thisFile <- read.delim(theRootDir %&% "data/genotypes/metaBric/mataBricGenotypesTxt/" %&% genoFiles[i], comment.char="#", as.is=T)
  genotypeList[[i]] <- thisFile[, 2]
  names(genotypeList[[i]]) <- thisFile[, 5]
  print(i)
}
names(genotypeList) <- genoFiles

genotypeMat <- do.call(cbind, genotypeList)
colnames(genotypeMat) <- substring(genoFiles, 13, 19)
# save(genotypeMat, file=theRootDir %&% "data/genotypes/metaBric/mataBricGenotypesTxt/metabricGenotypesChar.RData")
load(file=theRootDir %&% "data/genotypes/metaBric/mataBricGenotypesTxt/metabricGenotypesChar.RData")
genotypeMatNum <- genotypeMat

# First unphase the genotype matrix, i.e. make CA and AC both CA etc.
genotypeMatNum[genotypeMatNum == "CA"] <- "AC"
genotypeMatNum[genotypeMatNum == "GA"] <- "AG"
genotypeMatNum[genotypeMatNum == "TA"] <- "AT"
genotypeMatNum[genotypeMatNum == "GC"] <- "CG"
genotypeMatNum[genotypeMatNum == "TC"] <- "CT"
genotypeMatNum[genotypeMatNum == "TG"] <- "GT"
genotypeMatNum[genotypeMatNum == "---"] <- NA # use NA for the missing values.
reductedHets <- c("AC", "AG", "AT", "CG", "CT", "GT")
homozygotes <- c("AA", "CC", "TT", "GG")

# Make all Hets "1"
genotypeMatNum[genotypeMatNum %in% reductedHets] <- 1


# Find the call rate for each SNP
numNas <- apply(genotypeMatNum, 1, function(row){sum(is.na(row))}) # how many missing values are there per SNP 
callRate <- 1 - (numNas / ncol(genotypeMatNum))
hist(callRate, breaks=100)
sum(callRate > .95) / length(callRate) # [1] 0.8337265

# Find the call rate for each sample
numNasSamples <- apply(genotypeMatNum, 2, function(row){sum(is.na(row))}) # how many missing values are there per SNP 
callRateSamples <- 1 - (numNasSamples / nrow(genotypeMatNum))
hist(callRateSamples, breaks=100) # based on this plot I don't think there's a justification in removing samples because of call rate, there aren't any obvious outliers in the distribution.

# Keep SNPs with a > .95 call rate
genotypeMatNum <- genotypeMatNum[callRate > .95, ]

# Keep SNPs that have at least one heterozygote.
tabList <- apply(genotypeMatNum, 1, table)# Get a table of genotypes for each snp 
numGtypes <- sapply(tabList, length)
table(numGtypes)
hasHets <- sapply(tabList, function(item){return(names(item)[1] == "1")})

# Remove rows with no heterozygotes
genotypeMatNum <- genotypeMatNum[hasHets, ]
tabListNames <- lapply(tabList[hasHets], names)

# Convert the homozygotes to 1 or 2 (depending which appears first alphabetically).
for(i in 1:nrow(genotypeMatNum))
{
  genotypeMatNum[i, genotypeMatNum[i,] == tabListNames[[i]][2]] <- 0
  genotypeMatNum[i, genotypeMatNum[i,] == tabListNames[[i]][3]] <- 2
  print(i)
}
table(c(genotypeMatNum)) ## this hasn't worked fully, NB: there are still HETs that haven't been replaced, presumably there are > 2 genotypes for some of these, DELETE THESE ROWS.
#         0         1         2 
# 319983032 180206808 223866932

genotypeMatNumDat <- genotypeMatNum
class(genotypeMatNumDat) <- "numeric"

# Now remove anythign with a MAF < 0.05...
theMaf <- numeric(nrow(genotypeMatNumDat))
numSamps <- ncol(genotypeMatNumDat)
for(i in 1:nrow(genotypeMatNumDat))
{
  tab <- table(factor(genotypeMatNumDat[i,], levels=c(0,1,2))) # Count the number of 0, 1, 2 in each row of the genotype matrix.
  f1 <- ((tab["0"] * 2) + tab["1"]) / (sum(tab)*2) # count the total number of one of the alleles and divide by the total number alleles (total number of (0, 1, 2) x 2), to get the allele frequency for that allele
  f2 <- 1 - f1
  theMaf[i] <- min(c(f1, f2))
  print(i)
}

hist(theMaf, breaks=100)
print(sum(theMaf < 0.05)) # [1] 181368
hasAcceptableMaf <- rownames(genotypeMatNumDat)[theMaf > 0.05]
genotypeMatNumDat <- genotypeMatNumDat[hasAcceptableMaf, ]


# Calculate PCs of the Genotype matrix, to estimate anscestry (requires a complete matrix, so we'll use only SNPs that were called in every sample.)
numRowNas <- apply(genotypeMatNumDat, 1, function(x){sum(is.na(x))})
noNaSnps <- which(numRowNas == 0)
table(numRowNas)
genotypePcs <- prcomp(t(genotypeMatNumDat[noNaSnps , ]))

# save this finalized genotypeMatNumDat for METABRIC.
save(genotypeMatNumDat, genotypePcs, file=theRootDir %&% "data/genotypes/metaBric/mataBricGenotypesTxt/metabricGenotypesNumFin.RData")



