
#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

# Load the METABRIC results
load(file=theRootDir %&% "data/forMatrixEQTL/output/theData_metabric.RData") # theData_intInv_metabric, theData_conv_metabric, 
psMetabricConv <- theData_conv_metabric$pMatAll[1,]
fdrSigMetabricConv <- names(which(p.adjust(psMetabricConv, method="BH") < 0.05))
nomSigMetabricConv <- names(which(psMetabricConv < 0.05))
psMetabricIntInv <- theData_intInv_metabric$pMatAll[1,]
fdrSigMetabricIntInv <- names(which(p.adjust(psMetabricIntInv, method="BH") < 0.05))
nomSigMetabricIntInv <- names(which(psMetabricIntInv < 0.05))
length(fdrSigMetabricConv) # [1] [1] 47354
length(fdrSigMetabricIntInv) # [1] [1] 9235
sum(fdrSigMetabricIntInv %in% fdrSigMetabricConv) # [1] 8142



# Load the results for PEER, conventional approach
load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheDatas.RData") # theData_conv_Peer, theData_int_selectPeers, theData_conv_noPeer, theData_int_noPeer, theData_int_Peer, theData_intInv_noPeer, theData_intInv_Peer, theData_intInv_selectPeers, 
rm(list=(c("theData_int_selectPeers", "theData_conv_noPeer", "theData_int_noPeer", "theData_int_Peer", "theData_intInv_noPeer", "theData_intInv_selectPeers")))
psBulk <- theData_conv_Peer$pMatAll[1,]
names(psBulk) <- colnames(theData_conv_Peer$pMatAll)
betasBulk <- theData_conv_Peer$betaMatAll[1,]
names(betasBulk) <- colnames(theData_conv_Peer$betaMatAll)
fdrSigAllBrcaConv <- names(which(p.adjust(psBulk , method="BH") < 0.05))
nomSigAllBrcaConv <- names(which(psBulk < 0.05))

# Load the results for the cancer specific model.
psBulkIntInv <- theData_intInv_Peer$pMatAll[1,]
names(psBulkIntInv) <- colnames(theData_intInv_Peer$pMatAll)
betasBulkIntInv <- theData_intInv_Peer$betaMatAll[1,]
names(betasBulkIntInv) <- colnames(theData_intInv_Peer$betaMatAll)
fdrSigAllBrcaintInv <- names(which(p.adjust(psBulkIntInv, method="BH") < 0.05))
nomSigAllBrcaintInv <- names(which(psBulkIntInv < 0.05))

# Find the eQTLs that were assessed by BOTH TCGA and METABRIC.... (metabric was already subsetted to things in TCGA before the analysis)
commonEqtls <- names(psBulkIntInv)[names(psBulkIntInv) %in% colnames(theData_conv_metabric$pMatAll)]
fdrSigAllBrcaConv_com <- fdrSigAllBrcaConv[fdrSigAllBrcaConv %in% commonEqtls]
nomSigAllBrcaConv_com <- nomSigAllBrcaConv[nomSigAllBrcaConv %in% commonEqtls]
fdrSigAllBrcaintInv_com <- fdrSigAllBrcaintInv[fdrSigAllBrcaintInv %in% commonEqtls]
nomSigAllBrcaintInv_com <- nomSigAllBrcaintInv[nomSigAllBrcaintInv %in% commonEqtls]

# Create a vector indicating whether there is consistent direction of effect between TCGA and METABRIC eQTLs
# Note: there isn't any point in doing this because directionality wasn't presereved in the eQTL analysis....
betaVecConvMetabric <- theData_conv_metabric$betaMatAll[1,commonEqtls]
betaVecIntInvMetabric <- theData_intInv_metabric$betaMatAll[1,commonEqtls]
betaVecConvTcga <- theData_conv_Peer$betaMatAll[1,commonEqtls]
betaVecIntInvTcga <- theData_intInv_Peer$betaMatAll[1,commonEqtls]
directionVecConv <- betaVecConvMetabric * betaVecConvTcga # this will be positive if direction is the same.
directionVecIntInv <- betaVecIntInvMetabric * betaVecIntInvTcga

# How many of the TCGA results are also significant in METABRIC
sum(fdrSigAllBrcaintInv_com %in% nomSigMetabricIntInv) / length(fdrSigAllBrcaintInv_com) # [1] 0.5240725
sum(fdrSigAllBrcaConv_com %in% nomSigMetabricConv) / length(fdrSigAllBrcaConv_com) # [1] 0.5736259
sum(fdrSigAllBrcaintInv_com %in% fdrSigMetabricIntInv) / length(fdrSigAllBrcaintInv_com) # [1] 0.3147541
sum(fdrSigAllBrcaConv_com %in% fdrSigMetabricConv) / length(fdrSigAllBrcaConv_com) # [1] 0.3936987







