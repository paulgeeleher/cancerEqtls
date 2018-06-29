#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")


load(file=theRootDir %&% "data/forMatrixEQTL/output/theData_subtypeIntInvs.RData") # theData_intInv_Her2, theData_intInv_LuminalA, theData_intInv_LuminalB, theData_intInv_Basal

# I could plot a table with for each subtype, the proportions of FDR < 0.05 at least nominally significant in the other subtypes!? 
# And/or the entire data?
fdrSigHer2 <- names(which(p.adjust(theData_intInv_Her2$pMatAll[1, ], method="BH") < 0.05))
nomSigHer2 <- names(which(theData_intInv_Her2$pMatAll[1, ] < 0.05))
fdrSigLumA <- names(which(p.adjust(theData_intInv_LuminalA$pMatAll[1, ], method="BH") < 0.05))
nomSigLumA <- names(which(theData_intInv_LuminalA$pMatAll[1, ] < 0.05))
fdrSigLumB <- names(which(p.adjust(theData_intInv_LuminalB$pMatAll[1, ], method="BH") < 0.05))
nomSigLumB <- names(which(theData_intInv_LuminalB$pMatAll[1, ] < 0.05))
fdrSigBasal <- names(which(p.adjust(theData_intInv_Basal$pMatAll[1, ], method="BH") < 0.05))
nomSigBasal <- names(which(theData_intInv_Basal$pMatAll[1, ] < 0.05))

length(fdrSigHer2) # [1] 154
length(nomSigHer2) # [1] 190628
length(fdrSigLumA) # [1] 3041
length(nomSigLumA) # [1] 247684
length(fdrSigLumB) # [1] 219
length(nomSigLumB) # [1] 214525
length(fdrSigBasal) # [1] 44
length(nomSigBasal) # [1] 215565

# comparisons for her2
# her2 (numSig =  vs other
sum(fdrSigHer2 %in% nomSigLumA) # [1] 5
sum(fdrSigHer2 %in% nomSigLumB) # [1] 7
sum(fdrSigHer2 %in% nomSigBasal) # [1] 11

# Lum A (numsig = ) vs other
sum(fdrSigLumA %in% nomSigHer2) # [1] 474
sum(fdrSigLumA %in% nomSigLumB) # [1] 1394
sum(fdrSigLumA %in% nomSigBasal) # [1] 806
 
# Lum B (numsig = ) vs other
sum(fdrSigLumB %in% nomSigHer2) # [1] 46
sum(fdrSigLumB %in% nomSigLumA) # [1] 144
sum(fdrSigLumB %in% nomSigBasal) # [1] 107
 
# Basal (numsig = ) vs other
sum(fdrSigBasal %in% nomSigHer2) # [1] 26
sum(fdrSigBasal %in% nomSigLumA) # [1] 32
sum(fdrSigBasal %in% nomSigLumB) # [1] 31

# Numbers of significant associations in each subtype
length(fdrSigHer2) # [1] 154
length(nomSigHer2) # [1] 190628
length(fdrSigLumA) # [1] 3041
length(nomSigLumA) # [1] 247684
length(fdrSigLumB) # [1] 219
length(nomSigLumB) # [1] 214525
length(fdrSigBasal) # [1] 44
length(nomSigBasal) # [1] 215565
 



# The HER2 results are still having problems becuase of the MAF at the small sample size... (the outlier points will cause false positives because the assumptions of our statistical test are being broken)
# Plot the MAF of the "significant" hits against the MAF of everything.
# I will need to recalculate the MAFs in each subtype.....!!!!???
# This is going to be a little memory intensive.....
# Load the genotype data
load(file=theRootDir %&% "Results/rDatas/gwasFilteredInputData_WithQNexprAndCpe.RData") # tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin, geneLocation_filt_qnDat, theCpeProps_matched, covariatesVec_filt, expressionMatNorm_filt, expressionMat_filt, genotypesMat_filt, genotypePcs, snpsLocation_filt, geneLocation_filt, 
rm(list=c("tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin", "geneLocation_filt_qnDat", "theCpeProps_matched", "covariatesVec_filt", "expressionMatNorm_filt", "expressionMat_filt", "genotypePcs", "snpsLocation_filt", "geneLocation_filt"))

# Load PAM 50 subtypes and get a list of samples
pam50SubtypesFull <- read.csv(theRootDir %&% "tcgaPam50_fromNetanelyEtAl.csv", as.is=T, header=F, skip=1)
pam50SubtypesFull01 <- pam50SubtypesFull[substring(pam50SubtypesFull[,1], 14, 15) == "01", ]
table(pam50SubtypesFull01[, 7])
pam50SubtypesFullNorm <- pam50SubtypesFull01[, 7]
names(pam50SubtypesFullNorm) <- substring(pam50SubtypesFull01[,1], 1, 12)

# get them for different genotpes
aSamples <- intersect(substring(names(pam50SubtypesFullNorm)[which(pam50SubtypesFullNorm == "LumA")], 9, 12), colnames(genotypesMat_filt)) # samples for which there's both expression and PAM50 annotations
bSamples <- intersect(substring(names(pam50SubtypesFullNorm)[which(pam50SubtypesFullNorm == "LumB")], 9, 12), colnames(genotypesMat_filt))
basalSamples <- intersect(substring(names(pam50SubtypesFullNorm)[which(pam50SubtypesFullNorm == "Basal")], 9, 12), colnames(genotypesMat_filt))
her2Samples <- intersect(substring(names(pam50SubtypesFullNorm)[which(pam50SubtypesFullNorm == "Her2")], 9, 12), colnames(genotypesMat_filt))


getMafs <- function(mat) {
  theMaf <- numeric(nrow(mat))
  numSamps <- ncol(mat)
  for(i in 1:nrow(mat))
  {
    tab <- table(factor(mat[i,], levels=c(0,1,2))) # Count the number of 0, 1, 2 in each row of the genotype matrix.
    f1 <- ((tab["0"] * 2) + tab["1"]) / (sum(tab)*2) # count the total number of one of the alleles and divide by the total number alleles (total number of (0, 1, 2) x 2), to get the allele frequency for that allele
    f2 <- 1 - f1
    theMaf[i] <- min(c(f1, f2))
    print(i)
  }
  return(theMaf)
}

mafsHer2 <- getMafs(genotypesMat_filt[, her2Samples])
mafsLumA <- getMafs(genotypesMat_filt[, aSamples])
mafsLumB <- getMafs(genotypesMat_filt[, bSamples])
mafsBasal <- getMafs(genotypesMat_filt[, basalSamples])

# In truth the problem here is when you've an outlier more so than the MAF, constrain to at least 5% MAF the minor allele!!!!????' 
names(mafsHer2)  <- rownames(genotypesMat_filt)
names(mafsLumA) <- rownames(genotypesMat_filt)
names(mafsLumB) <- rownames(genotypesMat_filt)
names(mafsBasal) <- rownames(genotypesMat_filt)
fdrSigGenotypesHer2 <- unique(do.call(cbind, strsplit(fdrSigHer2, "~"))[2,])
fdrSigGenotypesLumA <- unique(do.call(cbind, strsplit(fdrSigLumA, "~"))[2,])
fdrSigGenotypesLumB <- unique(do.call(cbind, strsplit(fdrSigLumB, "~"))[2,])
fdrSigGenotypesBasal <- unique(do.call(cbind, strsplit(fdrSigBasal, "~"))[2,])
par(mfrow=c(4,2))
hist(mafsHer2, breaks=100)
hist(mafsHer2[fdrSigGenotypesHer2], breaks=100) # This is a big problem here, there's no way the things in this HER2 spike are valid
hist(mafsLumA, breaks=100)
hist(mafsLumA[fdrSigGenotypesLumA], breaks=100) # This is a big problem here, there's no way the things in this HER2 spike are valid
hist(mafsLumB, breaks=100)
hist(mafsLumB[fdrSigGenotypesLumB], breaks=100) # This is a big problem here, there's no way the things in this HER2 spike are valid
hist(mafsBasal, breaks=100)
hist(mafsBasal[fdrSigGenotypesBasal], breaks=100) # This is a big problem here, there's no way the things in this HER2 spike are valid

 
# We need to filter the SNPs that have a MAF of < 5% in the HER2 postive group (i.e. the smallest samples size), as these will
# cause us to break the assumptions of our linear model and thus produce a huge number of false postives at low MAF.
mafOk <- names(mafsHer2)[mafsHer2 > 0.05]
genotypeMapVec <- do.call(cbind, strsplit(colnames(theData_intInv_Her2$pMatAll), "~"))[2,] # This will allow me to index the snps that have a suitable MAF

# Calculate the updated FDRs, with the eQTLs very low MAF in the HER2 (smallest) group removed.
fdrSigHer2 <- names(which(p.adjust(theData_intInv_Her2$pMatAll[1, genotypeMapVec %in% mafOk], method="BH") < 0.05))
nomSigHer2 <- names(which(theData_intInv_Her2$pMatAll[1, ] < 0.05))
fdrSigLumA <- names(which(p.adjust(theData_intInv_LuminalA$pMatAll[1, genotypeMapVec %in% mafOk], method="BH") < 0.05))
nomSigLumA <- names(which(theData_intInv_LuminalA$pMatAll[1, ] < 0.05))
fdrSigLumB <- names(which(p.adjust(theData_intInv_LuminalB$pMatAll[1, genotypeMapVec %in% mafOk], method="BH") < 0.05))
nomSigLumB <- names(which(theData_intInv_LuminalB$pMatAll[1, ] < 0.05))
fdrSigBasal <- names(which(p.adjust(theData_intInv_Basal$pMatAll[1, genotypeMapVec %in% mafOk], method="BH") < 0.05))
nomSigBasal <- names(which(theData_intInv_Basal$pMatAll[1, ] < 0.05))

length(fdrSigHer2) # [1] 0
length(fdrSigLumA) # [1] 3038
length(fdrSigLumB) # [1] 211
length(fdrSigBasal) # [1] 43

# compare the results between subtypes when we remove very low MAF from HER2+
# Lum A (nsig = 3038)
sum(fdrSigLumA %in% nomSigHer2)  # [1] 476
sum(fdrSigLumA %in% nomSigLumB)  # [1] 1389
sum(fdrSigLumA %in% nomSigBasal)  # [1] 798

# Lum B (nsig = 211) vs other
sum(fdrSigLumB %in% nomSigHer2)  # [1] 45
sum(fdrSigLumB %in% nomSigLumA) # [1] 142
sum(fdrSigLumB %in% nomSigBasal) # [1] 104
 
# Basal (n = 43) vs other
sum(fdrSigBasal %in% nomSigHer2)  # [1] 26
sum(fdrSigBasal %in% nomSigLumA)  # [1] 31
sum(fdrSigBasal %in% nomSigLumB)  # [1] 30

# FDR sig in both subtypes
sum(fdrSigLumA %in% fdrSigHer2)  # [1] 0
sum(fdrSigLumA %in% fdrSigLumB)  # [1] 116
sum(fdrSigLumA %in% fdrSigBasal)  # [1] 24
 
# Lum B (nsig = 211) vs other
sum(fdrSigLumB %in% fdrSigHer2)  # [1] 0
sum(fdrSigLumB %in% fdrSigLumA) # [1] 116
sum(fdrSigLumB %in% fdrSigBasal) # [1] 8
  
# Basal (n = 43) vs other
sum(fdrSigBasal %in% fdrSigHer2)  # [1] 0
sum(fdrSigBasal %in% fdrSigLumA)  # [1] 24
sum(fdrSigBasal %in% fdrSigLumB)  # [1] 8




# Each versus All breast cancer pooled.
# I'm having some memory issues here, I'll have to get rid of sutff after its loaded.
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


# Are the eQTLs identified in the various subtypes also identified in the bulk data? Ans: Yes, almost all are.
# FDR sig in subtype AND FDR sig in ALL BRCA
sum(fdrSigHer2 %in% fdrSigAllBrcaintInv) # [1] 0
sum(fdrSigLumA %in% fdrSigAllBrcaintInv) # [1] 2372
sum(fdrSigLumB %in% fdrSigAllBrcaintInv) # [1] 142
sum(fdrSigBasal %in% fdrSigAllBrcaintInv) # [1] 30

# FDR sig in subtype and Nominally significant in All BRCA
sum(fdrSigHer2 %in% nomSigAllBrcaintInv) # [1] 0
sum(fdrSigLumA %in% nomSigAllBrcaintInv) # [1] 2911
sum(fdrSigLumB %in% nomSigAllBrcaintInv) # [1] 188
sum(fdrSigBasal %in% nomSigAllBrcaintInv) # [1] 41


# Output these results to a table for format: 
# Num Samples (n), Num Sig. (FDR < 0.05), num FDR Sig in Her2, num FDR Sig in LumA, num FDR Sig Lum B, num FDR Sig in Basal, num Nom Sig Her2, num Nom Sig LumA, num Nom sig LumB, num nom sig Basal.
ns <- table(pam50SubtypesFull01[, 7])

# YOU ARE HERE, FILL IN THE TABLE BELOW AND ADD TO THE MANUSCRIPT, CONSIDER A VENN DIAGRAM.

Her2_InteractionModel <- list("Num Samples"=ns["Her2"], "Num FDR < 0.05"=length(fdrSigHer2), "Also FDR < 0.05 HER2"= "NA", "Also FDR < 0.05 LumA"=0 , "Also FDR < 0.05 LumB"=0 , "Also FDR < 0.05 Basal"=0, "Also FDR < 0.05 BRCA"=sum(fdrSigHer2 %in% fdrSigAllBrcaintInv), "Also P < 0.05 HER2"= "NA", "Also P < 0.05 LumA"=0 , "Also P < 0.05 LumB"=0 , "Also P < 0.05 Basal"=0, "Also P < 0.05 BRCA"=sum(fdrSigHer2 %in% nomSigAllBrcaintInv))
LuminalA_InteractionModel <- list("Num Samples"=ns["LumA"], "Num FDR < 0.05"=length(fdrSigLumA), "Also FDR < 0.05 HER2"=0, "Also FDR < 0.05 LumA"= "NA", "Also FDR < 0.05 LumB"=sum(fdrSigLumA %in% fdrSigLumB), "Also FDR < 0.05 Basal"=sum(fdrSigLumA %in% fdrSigBasal), "Also FDR < 0.05 BRCA"=sum(fdrSigLumA %in% fdrSigAllBrcaintInv), "Also P < 0.05 HER2"=sum(fdrSigLumA %in% nomSigHer2), "Also P < 0.05 LumA"="NA" , "Also P < 0.05 LumB"=sum(fdrSigLumA %in% nomSigLumB), "Also P < 0.05 Basal"=sum(fdrSigLumA %in% nomSigBasal), "Also P < 0.05 BRCA"=sum(fdrSigLumA %in% nomSigAllBrcaintInv))
LuminalB_InteractionModel <- list("Num Samples"=ns["LumB"], "Num FDR < 0.05"=length(fdrSigLumB), "Also FDR < 0.05 HER2"=0, "Also FDR < 0.05 LumA"=sum(fdrSigLumB %in% fdrSigLumA), "Also FDR < 0.05 LumB"="NA", "Also FDR < 0.05 Basal"=sum(fdrSigLumB %in% fdrSigBasal), "Also FDR < 0.05 BRCA"=sum(fdrSigLumB %in% fdrSigAllBrcaintInv), "Also P < 0.05 HER2"=sum(fdrSigLumB %in% nomSigHer2), "Also P < 0.05 LumA"=sum(fdrSigLumB %in% nomSigLumA), "Also P < 0.05 LumB"="NA" , "Also P < 0.05 Basal"=sum(fdrSigLumB %in% nomSigBasal), "Also P < 0.05 BRCA"=sum(fdrSigLumB %in% nomSigAllBrcaintInv))
Basal_InteractionModel <- list("Num Samples"=ns["Basal"], "Num FDR < 0.05"=length(fdrSigBasal), "Also FDR < 0.05 HER2"=0, "Also FDR < 0.05 LumA"=sum(fdrSigBasal %in% fdrSigLumA), "Also FDR < 0.05 LumB"=sum(fdrSigBasal %in% fdrSigLumB), "Also FDR < 0.05 Basal"="NA", "Also FDR < 0.05 BRCA"=sum(fdrSigBasal %in% fdrSigAllBrcaintInv) , "Also P < 0.05 HER2"=sum(fdrSigBasal %in% nomSigHer2), "Also P < 0.05 LumA"=sum(fdrSigBasal %in% nomSigLumA), "Also P < 0.05 LumB"=sum(fdrSigBasal %in% nomSigLumB), "Also P < 0.05 Basal"="NA", "Also P < 0.05 BRCA"=sum(fdrSigBasal %in% nomSigAllBrcaintInv))

subtypeDfOut <- rbind(Her2_InteractionModel, LuminalA_InteractionModel, LuminalB_InteractionModel, Basal_InteractionModel)
write.csv(subtypeDfOut, theRootDir %&% "paper/suppTabSubtypes.csv")


################### DO THE SAME ANALYSIS AS ABOVE, BUT FOR A CONVETIONAL MODEL.
load(file=theRootDir %&% "data/forMatrixEQTL/output/theData_subtypeConv.RData") # theData_conv_Her2, theData_conv_LuminalA, theData_conv_LuminalB, theData_conv_Basal
fdrSigHer2_conv <- names(which(p.adjust(theData_conv_Her2$pMatAll[1, genotypeMapVec %in% mafOk], method="BH") < 0.05))
nomSigHer2_conv <- names(which(theData_conv_Her2$pMatAll[1, ] < 0.05))
fdrSigLumA_conv <- names(which(p.adjust(theData_conv_LuminalA$pMatAll[1, genotypeMapVec %in% mafOk], method="BH") < 0.05))
nomSigLumA_conv <- names(which(theData_conv_LuminalA$pMatAll[1, ] < 0.05))
fdrSigLumB_conv <- names(which(p.adjust(theData_conv_LuminalB$pMatAll[1, genotypeMapVec %in% mafOk], method="BH") < 0.05))
nomSigLumB_conv <- names(which(theData_conv_LuminalB$pMatAll[1, ] < 0.05))
fdrSigBasal_conv <- names(which(p.adjust(theData_conv_Basal$pMatAll[1, genotypeMapVec %in% mafOk], method="BH") < 0.05))
nomSigBasal_conv <- names(which(theData_conv_Basal$pMatAll[1, ] < 0.05))

length(fdrSigHer2_conv) # [1] 258
length(fdrSigLumA_conv) # [1] 33579
length(fdrSigLumB_conv) # [1] 3934
length(fdrSigBasal_conv) # [1] 2102


# compare the results between subtypes when we remove very low MAF from HER2+

# COMPARE THOSE THAT ARE NOMINALLY SIGNIFICANT p < 0.05
# HER2 (nsig = 258)
sum(fdrSigHer2_conv %in% nomSigLumA_conv)  # [1] 242
sum(fdrSigHer2_conv %in% nomSigLumB_conv)  # [1] 241
sum(fdrSigHer2_conv %in% nomSigBasal_conv)  # [1] 232
sum(fdrSigHer2_conv %in% nomSigAllBrcaConv) # 247

# Lum A (nsig = 33579)
sum(fdrSigLumA_conv %in% nomSigHer2_conv)  # [1] 6576
sum(fdrSigLumA_conv %in% nomSigLumB_conv)  # [1] 15979
sum(fdrSigLumA_conv %in% nomSigBasal_conv)  # [1] 11151
sum(fdrSigLumA_conv %in% nomSigAllBrcaConv) # [1] 32317

# Lum B (nsig = 3934) vs other
sum(fdrSigLumB_conv %in% nomSigHer2_conv)  # [1] 2005
sum(fdrSigLumB_conv %in% nomSigLumA_conv) # [1] 3588
sum(fdrSigLumB_conv %in% nomSigBasal_conv) #  [1] 2758
sum(fdrSigLumB_conv %in% nomSigAllBrcaConv) # [1] 3783

# Basal (n = 2102) vs other
sum(fdrSigBasal_conv %in% nomSigHer2_conv)  # [1] 1398
sum(fdrSigBasal_conv %in% nomSigLumA_conv)  # [1] 1902
sum(fdrSigBasal_conv %in% nomSigLumB_conv)  # [1] 1799
sum(fdrSigBasal_conv %in% nomSigAllBrcaConv) # [1] 2033


# COMPARE THOSE THAT ARE FDR SIGNIFINCANT
# HER2
sum(fdrSigHer2_conv %in% fdrSigLumA_conv)  # [1] 236
sum(fdrSigHer2_conv %in% fdrSigLumB_conv)  # [1] 146
sum(fdrSigHer2_conv %in% fdrSigBasal_conv)  # [1] 210
sum(fdrSigHer2_conv %in% fdrSigAllBrcaConv) # [1] 242

# Lum A (nsig = 33579)
sum(fdrSigLumA_conv %in% fdrSigHer2_conv)  # [1] 236
sum(fdrSigLumA_conv %in% fdrSigLumB_conv)  # [1] 3329
sum(fdrSigLumA_conv %in% fdrSigBasal_conv)  #  [1] 1799
sum(fdrSigLumA_conv %in% fdrSigAllBrcaConv) # [1] 27495

# Lum B (nsig = 3934) vs other
sum(fdrSigLumB_conv %in% fdrSigHer2_conv)  # [1] 146
sum(fdrSigLumB_conv %in% fdrSigLumA_conv) # [1] 3329
sum(fdrSigLumB_conv %in% nomSigBasal_conv) # [1] 2758
sum(fdrSigLumB_conv %in% fdrSigAllBrcaConv)# [1] 3616

# Basal (n = 2102) vs other
sum(fdrSigBasal_conv %in% fdrSigHer2_conv)  # [1] 210
sum(fdrSigBasal_conv %in% fdrSigLumA_conv)  # [1] 1799
sum(fdrSigBasal_conv %in% fdrSigLumB_conv)  # [1] 1166
sum(fdrSigBasal_conv %in% fdrSigAllBrcaConv) # [1] 1962

Her2_ConventionalModel <- list("Num Samples"=ns["Her2"], "Num FDR < 0.05"=length(fdrSigHer2_conv), "Also FDR < 0.05 HER2"= "NA", "Also FDR < 0.05 LumA"=sum(fdrSigHer2_conv %in% fdrSigLumA_conv) , "Also FDR < 0.05 LumB"=sum(fdrSigHer2_conv %in% fdrSigLumB_conv) , "Also FDR < 0.05 Basal"=sum(fdrSigHer2_conv %in% fdrSigBasal_conv), "Also FDR < 0.05 BRCA"=sum(fdrSigHer2_conv %in% fdrSigAllBrcaConv), "Also P < 0.05 HER2"= "NA", "Also P < 0.05 LumA"=sum(fdrSigHer2_conv %in% nomSigLumA_conv) , "Also P < 0.05 LumB"=sum(fdrSigHer2_conv %in% nomSigLumB_conv) , "Also P < 0.05 Basal"=sum(fdrSigHer2_conv %in% nomSigBasal_conv), "Also P < 0.05 BRCA"=sum(fdrSigHer2_conv %in% nomSigAllBrcaintInv))
LuminalA_ConventionalModel <- list("Num Samples"=ns["LumA"], "Num FDR < 0.05"=length(fdrSigLumA_conv), "Also FDR < 0.05 HER2"=sum(fdrSigLumA_conv %in% fdrSigHer2_conv), "Also FDR < 0.05 LumA"= "NA", "Also FDR < 0.05 LumB"=sum(fdrSigLumA_conv %in% fdrSigLumB_conv), "Also FDR < 0.05 Basal"=sum(fdrSigLumA_conv %in% fdrSigBasal_conv), "Also FDR < 0.05 BRCA"=sum(fdrSigLumA_conv %in% fdrSigAllBrcaConv), "Also P < 0.05 HER2"=sum(fdrSigLumA_conv %in% nomSigHer2_conv), "Also P < 0.05 LumA"="NA" , "Also P < 0.05 LumB"=sum(fdrSigLumA_conv %in% nomSigLumB_conv), "Also P < 0.05 Basal"=sum(fdrSigLumA_conv %in% nomSigBasal_conv), "Also P < 0.05 BRCA"=sum(fdrSigLumA_conv %in% nomSigAllBrcaConv))
LuminalB_ConventionalModel <- list("Num Samples"=ns["LumB"], "Num FDR < 0.05"=length(fdrSigLumB_conv), "Also FDR < 0.05 HER2"=sum(fdrSigLumB_conv %in% fdrSigHer2_conv), "Also FDR < 0.05 LumA"=sum(fdrSigLumB_conv %in% fdrSigLumA_conv), "Also FDR < 0.05 LumB"="NA", "Also FDR < 0.05 Basal"=sum(fdrSigLumB_conv %in% fdrSigBasal_conv), "Also FDR < 0.05 BRCA"=sum(fdrSigLumB_conv %in% fdrSigAllBrcaConv), "Also P < 0.05 HER2"=sum(fdrSigLumB_conv %in% nomSigHer2_conv), "Also P < 0.05 LumA"=sum(fdrSigLumB_conv %in% nomSigLumA_conv), "Also P < 0.05 LumB"="NA" , "Also P < 0.05 Basal"=sum(fdrSigLumB_conv %in% nomSigBasal_conv), "Also P < 0.05 BRCA"=sum(fdrSigLumB_conv %in% nomSigAllBrcaConv))
Basal_ConventionalModel <- list("Num Samples"=ns["Basal"], "Num FDR < 0.05"=length(fdrSigBasal_conv), "Also FDR < 0.05 HER2"=sum(fdrSigBasal_conv %in% fdrSigHer2_conv), "Also FDR < 0.05 LumA"=sum(fdrSigBasal_conv %in% fdrSigLumA_conv), "Also FDR < 0.05 LumB"=sum(fdrSigBasal_conv %in% fdrSigLumB_conv), "Also FDR < 0.05 Basal"="NA", "Also FDR < 0.05 BRCA"=sum(fdrSigBasal_conv %in% fdrSigAllBrcaConv) , "Also P < 0.05 HER2"=sum(fdrSigBasal_conv %in% nomSigHer2_conv), "Also P < 0.05 LumA"=sum(fdrSigBasal_conv %in% nomSigLumA_conv), "Also P < 0.05 LumB"=sum(fdrSigBasal_conv %in% nomSigLumB_conv), "Also P < 0.05 Basal"="NA", "Also P < 0.05 BRCA"=sum(fdrSigBasal_conv %in% nomSigAllBrcaConv))

subtypeDfOut_conv <- rbind(Her2_ConventionalModel, LuminalA_ConventionalModel, LuminalB_ConventionalModel, Basal_ConventionalModel)
write.csv(subtypeDfOut, theRootDir %&% "paper/suppTabSubtypes_conv.csv")

write.csv(subtypeDfOut_conv, "/home/ubuntu/Dropbox/predixcanProj/paper/SupplementaryTables/suppTabSubtypes_conv.csv")
write.csv(subtypeDfOut, "/home/ubuntu/Dropbox/predixcanProj/paper/SupplementaryTables/suppTabSubtypes.csv") 











