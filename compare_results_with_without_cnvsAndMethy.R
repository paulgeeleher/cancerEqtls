
#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

# Load the results controlled for CNV
load(file=theRootDir %&% "data/forMatrixEQTL/output/theData_cnvs.RData") #theData_conv_cnvs, theData_intInv_cnvs, theData_int_cnvs, 
psCnvConv <- theData_conv_cnvs$pMatAll[1,]
fdrSigCnvConv <- names(which(p.adjust(psCnvConv, method="BH") < 0.05))
psCnvIntInv <- theData_intInv_cnvs$pMatAll[1,]
fdrSigCnvIntInv <- names(which(p.adjust(psCnvIntInv, method="BH") < 0.05))
length(fdrSigCnvConv) # [1] 61198
length(fdrSigCnvIntInv) # [1] 9716

# Load the results controlled for methylation
load(file=theRootDir %&% "data/forMatrixEQTL/output/theData_methyl.RData") # theData_intInv_methyl, theData_conv_methyl, 
psMethylConv <- theData_conv_methyl$pMatAll[1,]
fdrSigMethylConv <- names(which(p.adjust(psMethylConv, method="BH") < 0.05))
psMethylIntInv <- theData_intInv_methyl$pMatAll[1,]
fdrSigMethylIntInv <- names(which(p.adjust(psMethylIntInv, method="BH") < 0.05))
length(fdrSigMethylConv) # [1] 
length(fdrSigMethylIntInv) # [1] 


# Load the results for PEER, conventional approach
load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheDatas.RData") # theData_conv_Peer, theData_int_selectPeers, theData_conv_noPeer, theData_int_noPeer, theData_int_Peer, theData_intInv_noPeer, theData_intInv_Peer, theData_intInv_selectPeers, 
rm(list=(c("theData_int_selectPeers", "theData_conv_noPeer", "theData_int_noPeer", "theData_int_Peer", "theData_intInv_noPeer", "theData_intInv_selectPeers")))
psBulk <- theData_conv_Peer$pMatAll[1,]
names(psBulk) <- colnames(theData_conv_Peer$pMatAll)
betasBulk <- theData_conv_Peer$betaMatAll[1,]
names(betasBulk) <- colnames(theData_conv_Peer$betaMatAll)
rm(theData_conv_Peer)
fdrSigAllBrcaConv <- names(which(p.adjust(psBulk , method="BH") < 0.05))
nomSigAllBrcaConv <- names(which(psBulk < 0.05))

# Load the results for the cancer specific model.
psBulkIntInv <- theData_intInv_Peer$pMatAll[1,]
names(psBulkIntInv) <- colnames(theData_intInv_Peer$pMatAll)
betasBulkIntInv <- theData_intInv_Peer$betaMatAll[1,]
names(betasBulkIntInv) <- colnames(theData_intInv_Peer$betaMatAll)
rm(theData_intInv_Peer)
fdrSigAllBrcaintInv <- names(which(p.adjust(psBulkIntInv, method="BH") < 0.05))
nomSigAllBrcaintInv <- names(which(psBulkIntInv < 0.05))

# PEER vs CNV
# There is strong overlap here, i.e. 45263 of a max of 57189 (80%) eQTL association are identical. 61198 associations were identified using the CNV approach
sum(fdrSigAllBrcaConv %in% fdrSigCnvConv) # [1] 45263
sum(nomSigAllBrcaConv %in% fdrSigCnvConv) # [1] 59882

# Also strong overlap here, 6461 of a max of 8833 remain significant following CNV correction
sum(fdrSigAllBrcaintInv %in% fdrSigCnvIntInv) # [1] 6461
sum(nomSigAllBrcaintInv %in% fdrSigCnvIntInv) # [1] 9575

# PEER vs Methylation, very strong overlap
sum(fdrSigAllBrcaConv %in% fdrSigMethylConv) # [1] 30552
sum(nomSigAllBrcaConv %in% fdrSigMethylConv) # [1] 33762

# Also strong overlap here, 6461 of a max of 8833 remain significant following CNV correction
sum(fdrSigAllBrcaintInv %in% fdrSigMethylIntInv) # 3864
sum(nomSigAllBrcaintInv %in% fdrSigMethylIntInv) # 4716

## VENN DIAGRAMs for comparison of results, include methylation also 
# Get a list of eQTLs that were included in all 3 analyses
commoneQTLsAll <- intersect(intersect(colnames(theData_conv_methyl$pMatAll), names(psBulk)), colnames(theData_conv_cnvs$pMatAll))

require(VennDiagram)
# conventional model (all 3): 
venn.diagram(list("PEER" = fdrSigAllBrcaConv[fdrSigAllBrcaConv %in% commoneQTLsAll], "CNVs" = fdrSigCnvConv[fdrSigCnvConv %in% commoneQTLsAll], "Methylation"=fdrSigMethylConv[fdrSigMethylConv %in% commoneQTLsAll]),fill = c("red", "green", "blue"), alpha = c(0.5, 0.5, .5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3, filename = "SUPPVENN_conv.tiff")

# conventional model, PEER only vs PEER + CNV
venn.diagram(list("PEER" = fdrSigAllBrcaConv[fdrSigAllBrcaConv %in% commoneQTLsAll], "CNVs" = fdrSigCnvConv[fdrSigCnvConv %in% commoneQTLsAll]),fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3, filename = "SUPPVENN_conv_peerVscnv.tiff")

# conventional model, PEER only vs PEER + methylation
venn.diagram(list("PEER" = fdrSigAllBrcaConv[fdrSigAllBrcaConv %in% commoneQTLsAll], "Methylation" = fdrSigMethylConv[fdrSigMethylConv %in% commoneQTLsAll]),fill = c("red", "blue"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3, filename = "SUPPVENN_conv_peerVsMethyl.tiff")

### interaction model (all 3):
venn.diagram(list("PEER" = fdrSigAllBrcaintInv[fdrSigAllBrcaintInv %in% commoneQTLsAll], "CNVs" = fdrSigCnvIntInv[fdrSigCnvIntInv %in% commoneQTLsAll], "Methylation"=fdrSigMethylIntInv[fdrSigMethylIntInv %in% commoneQTLsAll]),fill = c("red", "green", "blue"), alpha = c(0.5, 0.5, .5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3, filename = "SUPPVENN_intInv.tiff")

# interaction model (Peer vs cnv)
venn.diagram(list("PEER" = fdrSigAllBrcaintInv[fdrSigAllBrcaintInv %in% commoneQTLsAll], "CNVs" = fdrSigCnvIntInv[fdrSigCnvIntInv %in% commoneQTLsAll]),fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3, filename = "SUPPVENN_intInv_peerVsCnv.tiff")

# interaction model (Peer vs methyl)
venn.diagram(list("PEER" = fdrSigAllBrcaintInv[fdrSigAllBrcaintInv %in% commoneQTLsAll], "Methylation"=fdrSigMethylIntInv[fdrSigMethylIntInv %in% commoneQTLsAll]),fill = c("red", "blue"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3, filename = "SUPPVENN_intInv_peerVsMethylayion.tiff")
 








