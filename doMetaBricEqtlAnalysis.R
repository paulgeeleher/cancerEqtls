# This script runs the eQTL analysis for METABRIC, and outputs coefList data, which is processed  y "prepare_results_allAnalyses.R"

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")


# Load the required metabric data...
load(file=theRootDir %&% "data/genotypes/metaBric/metabricExpressionAndProps.RData") # cpePredictionsMetaFinOrd, expressionMatNormMetaOrd, 
load(file=theRootDir %&% "data/genotypes/metaBric/metabricExpressionAndProps_inTcga.RData") # genotypePcs, genotypeMatNumDat_inTcga
load(file=theRootDir %&% "data/genotypes/metaBric/metabricPcFactors.RData") # thePcFactors_NoInteractionAnalysis_metabric, thePcFactors_correctedAnaysis_metabric
load(file=theRootDir %&% "data/genotypes/metaBric/snpsToGenesListMeta.RData") # snpsToGenesListMeta
cpePredictionsMetaFinOrdInv <- 1 - cpePredictionsMetaFinOrd # I need to use the inverse of these proportions to recover the "cancer" component

coefsList_conv_Peer_metabric <- list(length(snpsToGenesListMeta))
coefsList_intInv_Peer_metabric <- list(length(snpsToGenesListMeta))
for(i in 1:length(snpsToGenesListMeta)) # [1] 615568
{
  coefsList_conv_Peer_metabric[[i]] <- coef(summary(lm(t(expressionMatNormMetaOrd[snpsToGenesListMeta[[i]], ,drop=FALSE])~genotypeMatNumDat_inTcga[names(snpsToGenesListMeta)[i], ]+genotypePcs$x[,1:3]+thePcFactors_NoInteractionAnalysis_metabric)))
  
  coefsList_intInv_Peer_metabric[[i]] <- coef(summary(lm(t(expressionMatNormMetaOrd[snpsToGenesListMeta[[i]], ,drop=FALSE])~(genotypeMatNumDat_inTcga[names(snpsToGenesListMeta)[i], ]*cpePredictionsMetaFinOrdInv)+genotypePcs$x[,1:3]+thePcFactors_correctedAnaysis_metabric)))
  print(i)
}
names(coefsList_conv_Peer_metabric) <- names(snpsToGenesListMeta)
names(coefsList_intInv_Peer_metabric) <- names(snpsToGenesListMeta)

save(coefsList_conv_Peer_metabric, coefsList_intInv_Peer_metabric, snpsToGenesListMeta, file=theRootDir %&% "data/forMatrixEQTL/output/allTheCoefLists_metabric.RData")





