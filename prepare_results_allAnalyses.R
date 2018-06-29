##' This script will take the results of the eQTL analysis, which are arranged in coefficient matrices,
##' and will arrange the p-values, beta values, standard errors, etc in easy to use matrices.

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

dir.create(theRootDir %&% "data/forMatrixEQTL/output/", recursive=TRUE, showWarnings = FALSE)

#' Function to convert the reuslts of the linear model fits to matrices of p-values and effect sizes, to make downstream analysis easier.
returnMatices <- function(theCoefList, rowsToGet, snpsToGenesList) # rowsToGet is the rows of the coef matrix that we want, e.g. the genotype, interaction term etc.
{
  #' I'd like a matrix of snp x gene with p-values, estimtes, which I will have to extract from "coefsList".
  numTargetGenes <- sapply(snpsToGenesList, length)
  withMultipleTargets <- which(numTargetGenes != 1) # I will have to run this separately for snps with 1 or more than 1 target gene, because the function has returned objects of different classes ("listof", or "matrix" if only one response variable)
  singleTarget <- which(numTargetGenes == 1)

  #' First run for snps with many possible gene targets
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
    
    #' get the std Errors in a matrix
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
    print(i)
  }
  colnames(pMat_sing) <- colnamesVec_sing
  colnames(betaMat_sing) <- colnamesVec_sing
  colnames(errMat_sing) <- colnamesVec_sing
  pMatAll <- cbind(pMat, pMat_sing)
  betaMatAll <- cbind(betaMat, betaMat_sing)
  errMatAll <- cbind(errMat, errMat_sing)
  
  return(list(pMatAll=pMatAll, betaMatAll=betaMatAll, errMatAll=errMatAll))
}



#' Function to convert the reuslts of the linear model fits to matrices of p-values and effect sizes, to make downstream analysis easier.
#' Note I'm using this function to pull out the results of the CNV analysis, because the data was returned in a slightly different format (for the cases where a snp was mapped toa  single gene.)
returnMatices2 <- function(theCoefList, rowsToGet, snpsToGenesList) # rowsToGet is the rows of the coef matrix that we want, e.g. the genotype, interaction term etc.
{
  #' I'd like a matrix of snp x gene with p-values, estimtes, which I will have to extract from "coefsList".
  numTargetGenes <- sapply(snpsToGenesList, length)
  withMultipleTargets <- which(numTargetGenes != 1) # I will have to run this separately for snps with 1 or more than 1 target gene, because the function has returned objects of different classes ("listof", or "matrix" if only one response variable)
  singleTarget <- which(numTargetGenes == 1)

  #' First run for snps with many possible gene targets
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
    
    #' get the std Errors in a matrix
    errs <- array(unlist(theCoefList[[i]]), dim=c(theMatSize,4,length(theCoefList[[i]])))[rowsToGet, 2, ]
    errMat[1:length(rowsToGet), (cumCols+1):(cumColsNew)] <- errs
    
    # colnamesVec[(cumCols+1):(cumColsNew)] <- paste(substring(names(theCoefList[[i]]), 10), "~", snpNamesOrd[i], sep="")
    colnamesVec[(cumCols+1):(cumColsNew)] <- paste(names(theCoefList[[i]]), "~", snpNamesOrd[i], sep="") # somehow this line needs to be changed....
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
    pMat_sing[,colNum] <- theCoefList[[i]][[1]][rowsToGet, 4]
    betaMat_sing[,colNum] <- theCoefList[[i]][[1]][rowsToGet, 1]
    errMat_sing[,colNum] <- theCoefList[[i]][[1]][rowsToGet, 2]
    colNum <- colNum + 1
    print(i)
  }
  colnames(pMat_sing) <- colnamesVec_sing
  colnames(betaMat_sing) <- colnamesVec_sing
  colnames(errMat_sing) <- colnamesVec_sing
  pMatAll <- cbind(pMat, pMat_sing)
  betaMatAll <- cbind(betaMat, betaMat_sing)
  errMatAll <- cbind(errMat, errMat_sing)
  
  return(list(pMatAll=pMatAll, betaMatAll=betaMatAll, errMatAll=errMatAll))
}




load(file=theRootDir %&% "Results/rDatas/snpsToGenesList.RData") # snpsToGenesList
load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheCoefLists.RData") # coefsList_conv_noPeer, coefsList_conv_Peer, coefsList_int_noPeer, coefsList_int_Peer, coefsList_int_selectPeers

rowsToGet <- 2
theData_conv_Peer <- returnMatices(coefsList_conv_Peer, rowsToGet)
rm(coefsList_conv_Peer)

rowsToGet <- c(2,3, 34)
theData_int_selectPeers <- returnMatices(coefsList_int_selectPeers, rowsToGet)
rm(coefsList_int_selectPeers)

rowsToGet <- 2
theData_conv_noPeer <- returnMatices(coefsList_conv_noPeer, rowsToGet)
rm(coefsList_conv_noPeer)

rowsToGet <- c(2, 3, 7)
theData_int_noPeer <- returnMatices(coefsList_int_noPeer, rowsToGet)
rm(coefsList_int_noPeer)

rowsToGet <- c(2, 3, 42)
theData_int_Peer <- returnMatices(coefsList_int_Peer, rowsToGet)
rm(coefsList_int_Peer)

load(file=theRootDir %&% "data/forMatrixEQTL/output/allInvCoefLists.RData") # coefsList_intInv_noPeer, coefsList_intInv_Peer, coefsList_intInv_selectPeers, 

rowsToGet <- c(2, 3, 7)
theData_intInv_noPeer <- returnMatices(coefsList_intInv_noPeer, rowsToGet)
rm(coefsList_intInv_noPeer)

rowsToGet <- c(2, 3, 42)
theData_intInv_Peer <- returnMatices(coefsList_intInv_Peer, rowsToGet)
rm(coefsList_intInv_Peer)

rowsToGet <- c(2, 3, 34)
theData_intInv_selectPeers <- returnMatices(coefsList_intInv_selectPeers, rowsToGet)
rm(coefsList_intInv_selectPeers)

save(theData_conv_Peer, theData_int_selectPeers, theData_conv_noPeer, theData_int_noPeer, theData_int_Peer, theData_intInv_noPeer, theData_intInv_Peer, theData_intInv_selectPeers, file=theRootDir %&% "data/forMatrixEQTL/output/allTheDatas.RData")

rowsToGet <- c(2, 3, 42)
theData_intInv_Peer_RandomCpe <- returnMatices(coefList_intInv_Peer_RandomCpe, rowsToGet)
save(theData_intInv_Peer_RandomCpe, file=theRootDir %&% "data/forMatrixEQTL/output/allTheDatas_theData_intInv_Peer_RandomCpe.RData")


#' The models that were fit on the 10% of the data with the most cancer cell content.....:
load(file=theRootDir %&% "data/forMatrixEQTL/output/coefList10pct.RData") # coefsList_conv_Peer_10pct
rowsToGet <- 2
theData_conv_10pct <- returnMatices(coefsList_conv_Peer_10pct, rowsToGet) # length: 615568
save(theData_conv_10pct, file=theRootDir %&% "data/forMatrixEQTL/output/theData_conv_10pct.RData")

#' The models (conventional and interaction) that control for CNVs
load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheCoefLists_cnvs.RData") # coefsList_conv_cnvs, coefsList_int_cnvs.
load(file=theRootDir %&% "data/forMatrixEQTL/output/coefsList_intInv_cnvs.RData") # coefsList_intInv_cnvs
load(file=theRootDir %&% "data/forMatrixEQTL/output/snpsToGenesList_cnvData.RData") # snpsToGenesList_cnvData ## I need the version of this for the CNV data.

rowsToGet <- c(2, 41)
theData_conv_cnvs <- returnMatices2(coefsList_conv_cnvs, rowsToGet, snpsToGenesList_cnvData)
sum(p.adjust(theData_conv_cnvs$pMatAll[1,], method="BH") < 0.05) # [1] 61198

#' The results for the "cancer specific" eQTLs, also increases very slightly.
rowsToGet <- c(2, 42, 43)
theData_intInv_cnvs <- returnMatices2(coefsList_intInv_cnvs, rowsToGet, snpsToGenesList_cnvData)
sum(p.adjust(theData_intInv_cnvs$pMatAll[1,], method="BH") < 0.05) # [1] 9716

#' Normal specific eqtls with CNV.
rowsToGet <- c(2, 42, 43)
theData_int_cnvs <- returnMatices2(coefsList_int_cnvs, rowsToGet, snpsToGenesList_cnvData)
sum(p.adjust(theData_int_cnvs$pMatAll[1,], method="BH") < 0.05) # [1] 468

save(theData_conv_cnvs, theData_intInv_cnvs, theData_int_cnvs, file=theRootDir %&% "data/forMatrixEQTL/output/theData_cnvs.RData")


# Breast cancer subtype specific results, from the PAM50 signatures.
# load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheCoefLists.RData") # coefsList_conv_LuminalA, coefsList_conv_LuminalB, coefsList_conv_Basal, coefsList_conv_Her2, coefsList_int_LuminalA, coefsList_int_LuminalB, coefsList_int_Basal, coefsList_int_Her2
load(file=theRootDir %&% "data/forMatrixEQTL/output/subtypeIntInvCoefLists.RData") # coefsList_intInv_LuminalA, coefsList_intInv_LuminalB, coefsList_intInv_Basal, coefsList_intInv_Her2, 

# rowsToGet <- c(2, 3, 27)
# theData_intInv_Her2 <- returnMatices(coefsList_intInv_Her2, rowsToGet) # for some reason this bombed out on i = 2086 cause the interaction never got calculated for that combo, WHY!!!? No variability in genotype causes this to drop out.

# New values from PAM50 from that paper
#  Basal   Her2   LumA   LumB Normal 
#    183     78    534    203     37

# HER2
rowsToGet <- c(2, 3)
theData_intInv_Her2 <- returnMatices(coefsList_intInv_Her2, rowsToGet)
sum(p.adjust(theData_intInv_Her2$pMatAll[1, ], method="BH") < 0.05) # # new sig results = 154 (this seems more reasonable!). Was: [1] 1155, this is probably just because of outliers (not enough samples)......!!!!?????

# Luminal A
rowsToGet <- c(2, 3)
theData_intInv_LuminalA <- returnMatices(coefsList_intInv_LuminalA, rowsToGet)
sum(p.adjust(theData_intInv_LuminalA$pMatAll[1, ], method="BH") < 0.05) # new sig results = 3041. Was: [1] 352

# Luminal B
rowsToGet <- c(2, 3)
theData_intInv_LuminalB <- returnMatices(coefsList_intInv_LuminalB, rowsToGet)
sum(p.adjust(theData_intInv_LuminalB$pMatAll[1, ], method="BH") < 0.05) # new num sig results = 219; Was: [1] 25

# Basal
rowsToGet <- c(2, 3)
theData_intInv_Basal <- returnMatices(coefsList_intInv_Basal, rowsToGet)
sum(p.adjust(theData_intInv_Basal$pMatAll[1, ], method="BH") < 0.05) # new num sig results = [1] 14

# Save the subtype specific data, for the intInv data (i.e. the cancer/interaction model)
save(theData_intInv_Her2, theData_intInv_LuminalA, theData_intInv_LuminalB, theData_intInv_Basal, file=theRootDir %&% "data/forMatrixEQTL/output/theData_subtypeIntInvs.RData")



#' Do the same as above for the conventional subtype analysis.
load(file=theRootDir %&% "data/forMatrixEQTL/output/subtype_conv_CoefLists.RData") # coefsList_conv_LuminalA, coefsList_conv_LuminalB, coefsList_conv_Basal, coefsList_conv_Her2
rowsToGet <- c(2)
theData_conv_Her2 <- returnMatices(coefsList_conv_Her2, rowsToGet)
sum(p.adjust(theData_conv_Her2$pMatAll[1, ], method="BH") < 0.05) # [1] 273

# Luminal A conventional
rowsToGet <- c(2)
theData_conv_LuminalA <- returnMatices(coefsList_conv_LuminalA, rowsToGet)
sum(p.adjust(theData_conv_LuminalA$pMatAll[1, ], method="BH") < 0.05) # [1] 33883

# Luminal B conventional
rowsToGet <- c(2)
theData_conv_LuminalB <- returnMatices(coefsList_conv_LuminalB, rowsToGet)
sum(p.adjust(theData_conv_LuminalB$pMatAll[1, ], method="BH") < 0.05) # [1] 3941

# Basal conventional
rowsToGet <- c(2)
theData_conv_Basal <- returnMatices(coefsList_conv_Basal, rowsToGet)
sum(p.adjust(theData_conv_Basal$pMatAll[1, ], method="BH") < 0.05) # [1] 2119

# Save the subtype specific data, for the conv data (i.e. the cancer/interaction model)
save(theData_conv_Her2, theData_conv_LuminalA, theData_conv_LuminalB, theData_conv_Basal, file=theRootDir %&% "data/forMatrixEQTL/output/theData_subtypeConv.RData") 



## Finally, they methylation analysis. Unsurprisingly, this causes a loss of power.
rowsToGet <- c(2)
theData_intInv_methyl <- returnMatices2(coefsList_intInv_methys, rowsToGet, snpsToGenesList_methyl)
sum(p.adjust(theData_intInv_methyl$pMatAll[1,], method="BH") < 0.05) # 

#' cancer specific" eQTLs
rowsToGet <- c(2)
theData_conv_methyl <- returnMatices2(coefsList_conv_methys, rowsToGet, snpsToGenesList_methyl)
sum(p.adjust(theData_conv_methyl$pMatAll[1,], method="BH") < 0.05) #

save(theData_intInv_methyl, theData_conv_methyl, file=theRootDir %&% "data/forMatrixEQTL/output/theData_methyl.RData")


## Even more finally, the METABRIC data....
load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheCoefLists_metabric.RData") # coefsList_conv_Peer_metabric, coefsList_intInv_Peer_metabric, snpsToGenesListMeta
rowsToGet <- c(2)
theData_intInv_metabric <- returnMatices(coefsList_intInv_Peer_metabric, rowsToGet, snpsToGenesListMeta)
sum(p.adjust(theData_intInv_metabric$pMatAll[1,], method="BH") < 0.05) # [1] 9235

#' The results for the "cancer specific" eQTLs, also increases very slightly.
rowsToGet <- c(2)
theData_conv_metabric <- returnMatices(coefsList_conv_Peer_metabric, rowsToGet, snpsToGenesListMeta)
sum(p.adjust(theData_conv_metabric$pMatAll[1,], method="BH") < 0.05) # [1] 47354

save(theData_intInv_metabric, theData_conv_metabric, file=theRootDir %&% "data/forMatrixEQTL/output/theData_metabric.RData")


print(sessionInfo())




