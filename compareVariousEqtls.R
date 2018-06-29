##' In this script I will attempt to compare the results that I got between the different methods amd present some results. 

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

dir.create(theRootDir %&% "paper/figures/figure2", recursive=TRUE, showWarnings = FALSE)

#' Load the results of the various eQTL analysis. (these we3re created in: prepare_results_allAnalyses.R)
load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheDatas.RData") # theData_conv_Peer, theData_int_selectPeers, theData_conv_noPeer, theData_int_noPeer, theData_int_Peer, theData_intInv_noPeer, theData_intInv_Peer, theData_intInv_selectPeers
load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheDatas_theData_intInv_Peer_RandomCpe.RData") # theData_intInv_Peer_RandomCpe


#' The results of a conventional eQTL analysis with PEER factors
adjs_conv_Peer <- p.adjust(theData_conv_Peer[[1]][1,], method="BH")
sum(adjs_conv_Peer < 0.05) 
sig_conv_peer <- names(adjs_conv_Peer)[adjs_conv_Peer < 0.05]

# Make a QQ-plot of these P-values.
png(paste(theRootDir, "paper/figures/figure2/qqPlotConventional.png", sep=""), height=900, width=900, res=300)
plot(-log10(1:length(theData_conv_Peer[[1]][1,]) / length(theData_conv_Peer[[1]][1,])), sort(-log10(theData_conv_Peer[[1]][1,]), decreasing=T), xlab="-Log10 P-Value (uniform distribution)", ylab="-Log10 P-Value (measured)", cex.axis=0.8, las=1, col="#00000066", pch=20, bty="l")
abline(0,1, col="red")
dev.off()

#' The results of a conventional eQTL analysis without PEER factors
adjs_conv_noPeer <- p.adjust(theData_conv_noPeer[[1]][1,], method="BH")
sum(adjs_conv_noPeer < 0.05) 
sig_conv_noPeer <- names(adjs_conv_noPeer)[adjs_conv_noPeer < 0.05]

#' The results of an interaction model eQTL analysis without PEER factors. This model is set up so the the main effect recovers an eQTL in "normal", i.e. "normal" is encoded as 0.
adjs_int_noPeer <- p.adjust(theData_int_noPeer[[1]][1,], method="BH")
sum(adjs_int_noPeer < 0.05) 
sig_int_noPeer <- names(adjs_int_noPeer)[adjs_int_noPeer < 0.05]

#' The results of an interaction model eQTL analysis with selected PEER factors. This model is set up so the the main effect recovers an eQTL in "normal", i.e. "normal" is encoded as 0.
adjs_int_selectPeers <- p.adjust(theData_int_selectPeers[[1]][1,], method="BH")
sum(adjs_int_selectPeers < 0.05) 
sig_int_selectPeers <- names(adjs_int_selectPeers)[adjs_int_selectPeers < 0.05]

#' The results of an interaction model eQTL analysis with all PEER factors. This model is set up so the the main effect recovers an eQTL in "normal", i.e. "normal" is encoded as 0.
adjs_int_Peer <- p.adjust(theData_int_Peer[[1]][1,], method="BH")
sum(adjs_int_Peer < 0.05) 
sig_int_Peer <- names(adjs_int_Peer)[adjs_int_Peer < 0.05]

#' The results of an interaction model eQTL analysis with all PEER factors. This model is set up so the the main effect recovers an eQTL in "cancer", i.e. "cancer" is encoded as 0.
adjs_intInv_noPeer <- p.adjust(theData_intInv_noPeer[[1]][1,], method="BH")
sum(adjs_intInv_noPeer < 0.05) 
sig_intInv_noPeer <- names(adjs_intInv_noPeer)[adjs_intInv_noPeer < 0.05]

#' The results of an interaction model eQTL analysis with all PEER factors. This model is set up so the the main effect recovers an eQTL in "cancer", i.e. "cancer" is encoded as 0.
adjs_intInv_Peer <- p.adjust(theData_intInv_Peer[[1]][1,], method="BH")
sum(adjs_intInv_Peer < 0.05) 
sig_intInv_Peer <- names(adjs_intInv_Peer)[adjs_intInv_Peer < 0.05]

#' The results of the model with random tumor purity estimates
adjs_intInvRandom_Peer <- p.adjust(theData_intInv_Peer_RandomCpe[[1]][1,], method="BH")
sum(adjs_intInvRandom_Peer < 0.05) 
sigRandomStuff <- names(adjs_intInvRandom_Peer)[adjs_intInvRandom_Peer < 0.05]



# Make a QQ-plot of these P-values.
png(paste(theRootDir, "paper/figures/figure2/qqPlotInteraction.png", sep=""), height=900, width=900, res=300)
plot(-log10(1:length(theData_intInv_Peer[[1]][1,]) / length(theData_intInv_Peer[[1]][1,])), sort(-log10(theData_intInv_Peer[[1]][1,]), decreasing=T), xlab="-Log10 P-Value (uniform distribution)", ylab="-Log10 P-Value (measured)", cex.axis=0.8, las=1, col="#00000066", pch=20, bty="l")
abline(0,1, col="red")
dev.off()

#' Make a QQ-plot that overlays the Bulk and Interaction model results.
convSort <- sort(-log10(theData_conv_Peer[[1]][1,]), decreasing=T)
intSort <- sort(-log10(theData_intInv_Peer[[1]][1,]), decreasing=T)
png(paste(theRootDir, "paper/figures/figure2/qqPlotBoth.png", sep=""), height=900, width=1500, res=300)
plot(-log10(1:length(theData_conv_Peer[[1]][1,]) / length(theData_conv_Peer[[1]][1,])), convSort, xlab="-Log10 P-Value (uniform distribution)", ylab="-Log10 P-Value (measured)", cex.axis=0.8, las=1, col="#1a964166", pch=20, bty="l")
points(-log10(1:length(theData_intInv_Peer[[1]][1,]) / length(theData_intInv_Peer[[1]][1,])), intSort, pch=20, col="#1f78b466")
abline(0,1, col="red")
dev.off()


#' Make a QQ-plot that overlays the Bulk and Interaction model results. And the reuslts wiht random CPE estimates....
convSort <- sort(-log10(theData_conv_Peer[[1]][1,]), decreasing=T)
intSort <- sort(-log10(theData_intInv_Peer[[1]][1,]), decreasing=T)
intSortRand <- sort(-log10(theData_intInv_Peer_RandomCpe[[1]][1,]), decreasing=T)
png(paste(theRootDir, "paper/figures/figure2/qqPlot_three.png", sep=""), height=900, width=1500, res=300)
plot(-log10(1:length(theData_conv_Peer[[1]][1,]) / length(theData_conv_Peer[[1]][1,])), convSort, xlab="-Log10 P-Value (uniform distribution)", ylab="-Log10 P-Value (measured)", cex.axis=0.8, las=1, col="#1a964166", pch=20, bty="l")
points(-log10(1:length(theData_intInv_Peer[[1]][1,]) / length(theData_intInv_Peer[[1]][1,])), intSort, pch=20, col="#1f78b466")
points(-log10(1:length(theData_intInv_Peer_RandomCpe[[1]][1,]) / length(theData_intInv_Peer_RandomCpe[[1]][1,])), intSortRand, pch=20, col="#ff7f0066")
abline(0,1, col="red")
dev.off()


svg(paste(theRootDir, "paper/figures/figure2/points.svg", sep=""), width=3, height=3)
plot(c(1,2), c(2,2), col=c("#1f78b4", "#1a9641"), pch=20)
dev.off()

#' The results of an interaction model eQTL analysis with select PEER factors. This model is set up so the the main effect recovers an eQTL in "cancer", i.e. "cancer" is encoded as 0.
adjs_intInv_selectPeers <- p.adjust(theData_intInv_selectPeers[[1]][1,], method="BH")
sum(adjs_intInv_selectPeers < 0.05) 
sig_intInv_selectPeers <- names(adjs_intInv_selectPeers)[adjs_intInv_selectPeers < 0.05]

#' Most of the things identified specific to cancer or normal cells are also identified in the bulk data.
sum(sig_intInv_Peer %in% sig_conv_peer) 
sum(sig_int_Peer %in% sig_conv_peer)

# What do the interaction term P-values look like for the bulk eQTLs? About 10% have an interaction P-value of < 0.05...
interactionsPsForSigBulk <- theData_intInv_Peer[[1]][3,sig_conv_peer]
sum(interactionsPsForSigBulk > 0.05) 
sum(interactionsPsForSigBulk > 0.05)
hist(interactionsPsForSigBulk, breaks=100)
sum(p.adjust(interactionsPsForSigBulk, method="BH") < 0.05) # [1] 61



# What do the purity P-values look like for the bulk eQTLs?
purityPsForSigBulk <- theData_intInv_Peer[[1]][2,sig_conv_peer]
sum(purityPsForSigBulk < 0.05) # [1] 6053
sum(purityPsForSigBulk > 0.05) # [1] 52639
hist(purityPsForSigBulk, breaks=100)

# How many of the eQTLs that are significant in the bulk analysis remain significant if we only consider these eQTLs in the "Cancer" analysis. Is Directionality the same (answer: it always is for the things that remain significant)?
cancerPsForSigBulk <- theData_intInv_Peer[[1]][1,sig_conv_peer]
sum(p.adjust(cancerPsForSigBulk, method="BH") < 0.05) # [1] 29330
stillSigInteractionModel <- which(p.adjust(cancerPsForSigBulk, method="BH") < 0.05)
effectProduct <- theData_intInv_Peer[["betaMatAll"]][1,sig_conv_peer] * theData_conv_Peer[["betaMatAll"]][1,sig_conv_peer] 
effectNotConcordant <- which(effectProduct < 0)
sum(effectNotConcordant %in% stillSigInteractionModel) # [1] 0

# Overall, only slightly over 5% have an interaction P-value of < 0.05.
sum(theData_int_Peer[[1]][3,] < 0.05) / length(theData_int_Peer[[1]][3,] < 0.05) # [1] 0.05374158

# 61 of these interaction p-values are significant following correction for multiple testing.
sum(p.adjust(interactionsPsForSigBulk, method="BH") < 0.05) # [1] 61
eqtlsWithSigInt <- names(which(p.adjust(interactionsPsForSigBulk, method="BH") < 0.05))

# What's the tendency for directionality in these interactions? Loss or Gain of effect size?
effectSizesForIntEqtls <- theData_intInv_Peer[[2]][3,eqtlsWithSigInt]
sum(effectSizesForIntEqtls < 0) 
sum(effectSizesForIntEqtls > 0)

#' In general, what are these genes? Are they simply genes that gain/lose expression, thus gain / lose eQTL? Or maybe they are regulated by a TF switched on in cancer etc.
genesWithCancerSpecificEqtl <- unique(do.call(cbind, strsplit(names(effectSizesForIntEqtls[effectSizesForIntEqtls < 0]), "~"))[1,])
genesWithNormalSpecificEqtl <- unique(do.call(cbind, strsplit(names(effectSizesForIntEqtls[effectSizesForIntEqtls > 0]), "~"))[1,])
write.table(genesWithCancerSpecificEqtl, row.names=F, quote=F) # not much happening here as regards processes/pathways.
write.table(genesWithNormalSpecificEqtl, row.names=F, quote=F)
allSnpGenes <- colnames(theData_conv_Peer[[1]])

## Show a specific eQTL!
# List the top 50 eQTLs in Bulk, and their corresponding P-value in Cancer.
topPvalsSortBulk <- sort(-log10(theData_conv_Peer[[1]][1,])
intmat <- cbind(theData_conv_Peer[[1]][1,eqtlsWithSigInt], theData_intInv_Peer[[1]][1,eqtlsWithSigInt])
intmat["MDGA1~rs6458012", ] # 1.543865e-29 2.635636e-01

#' Lets use MDGA1~rs6458012 as an example.
load(file=theRootDir %&% "Results/rDatas/gwasFilteredInputData_WithQNexprAndCpe.RData") # tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin, geneLocation_filt_qnDat, theCpeProps_matched, covariatesVec_filt, expressionMatNorm_filt, expressionMat_filt, genotypesMat_filt, genotypePcs, snpsLocation_filt, geneLocation_filt, 


svg(file=paste(theRootDir, "paper/figures/figure2/Fig2_MDGA1_rs6458012.svg", sep=""), width=3, height=3)
plot((genotypesMat_filt["rs6458012", ]+rnorm(length(genotypesMat_filt["rs6458012", ]), 0, .1)), tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin["MDGA1",], col="#00000044", pch=20, las=1, ylab="MDGA1 expression", xlab="rs6458012", bty="l", xaxt = "n", cex.axis=.8) 
axis(1, at=0:2, labels=c("A/A","A/G","G/G"), cex.axis=.8)
abline(coef=coef(summary(lm(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin["MDGA1",]~genotypesMat_filt["rs6458012", ])))[,1], col="#e31a1c")
dev.off()


## Make the NewQTL plot.
theCpeProps_matched_INV <- 1 - theCpeProps_matched
d <- sort(theCpeProps_matched_INV)
a <- split(d, ceiling(seq_along(d)/179))

theSes <- numeric()
theBetas <- numeric()
load(file=theRootDir %&% "Results/rDatas/thePeers.RData") # theFactors_NoInteractionAnalysis, theFactors_correctedAnaysis, theFactors_correctedAnaysis_filtered, 

rownames(theFactors_correctedAnaysis) <- colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin)
rownames(theFactors_NoInteractionAnalysis) <- colnames(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin)

# theEqtl <- "ARHGEF5~rs720475"
makeNewQtlsPlot <- function(theEqtl)
{
  theGene <- strsplit(theEqtl, "~")[[1]][1]
  theSnp <- strsplit(theEqtl, "~")[[1]][2]
  theBetas <- numeric()
  theSes <- numeric()
  for(i in 1:length(a))
  {
    theCoefHere <- coef(summary(lm(t(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[theGene, names(a[[i]]),drop=FALSE])~genotypesMat_filt[theSnp, names(a[[i]])]+theCpeProps_matched_INV[names(a[[i]])]+genotypePcs$x[names(a[[i]]),1:3]+theFactors_NoInteractionAnalysis[names(a[[i]]), ])))
    theSes[i] <- theCoefHere[2, 2] *1.96
    theBetas[i] <- theCoefHere[2, 1]
#     print(i)
  }

  seNormal <- theData_int_Peer$errMatAll[1, theEqtl] *1.96
  betaNormal <- theData_int_Peer$betaMatAll[1, theEqtl]
  seCancer <- theData_intInv_Peer$errMatAll[1, theEqtl] *1.96
  betaCancer <- theData_intInv_Peer$betaMatAll[1, theEqtl]
  betaMaxCancer <- betaCancer + seCancer
  betaMinCancer <- betaCancer - seCancer
  betaMaxNormal <- betaNormal + seNormal
  betaMinNormal <- betaNormal - seNormal 
  
  maxY <- max(c(betaMaxCancer, betaMaxNormal, theBetas+theSes)) + .1
  minY <- min(c(betaMinCancer, betaMinNormal,theBetas-theSes)) - .1

  x <- sapply(a, mean)

  plot(x, theBetas, ylim=c(minY, 1.05), xlim=c(-.28, .7), pch=20, xlab="Proportion of Normal Tissue", ylab="eQTL Effect Size", main="", las=1, cex.axis=.8, xaxt = "n", bty="l",panel.first=abline(h=0, col="darkgrey", lty=2))
  axis(1, at=c(-.2, 0, .2, .4, .6), labels=c("Bulk\nTumor", "0%", "20%", "40%", "60%"), cex.axis=.8)
  
  abline(v=-.1, col="grey", lty=1, lwd=1)
  arrows(x, theBetas-theSes, x, theBetas+theSes, length=0.05, angle=90, code=3)
  
  # Plot the effect size and standard error for the bulk tumor
  betaBulk <- theData_conv_Peer$betaMatAll[1, theEqtl]
  seBulk <- theData_conv_Peer$errMatAll[1, theEqtl]*1.96
  points(-.2, theData_conv_Peer$betaMatAll[1, theEqtl], col="#1a9641", pch=20)
  arrows(-.2, betaBulk-seBulk, -.2, betaBulk+seBulk, length=0.05, angle=90, code=3, col="#1a9641")
  
  # add the points for "cancer" and "normal"
  points(0, betaCancer, pch=20, col="#ca0020")
#   points(1, betaNormal, pch=20, col="#2c7bb6")
  
  # add the standard errors for "cancer" and "normal"
  arrows(0, betaMinCancer, 0, betaMaxCancer, length=0.05, angle=90, code=3, col="#ca0020") # cancer
#   arrows(1, betaMinNormal, 1, betaMaxNormal, length=0.05, angle=90, code=3, col="#2c7bb6") # normal
  

#   lines(c( ,), c( , ))
#   points(sapply(a, mean), theBetas+theSes, col="red")
#   points(sapply(a, mean), theBetas-theSes, col="blue")
}

svg(file=paste(theRootDir, "paper/figures/figure2/Fig2_MDGA1_rs6458012_newqtl.svg", sep=""), width=4, height=3.75)
makeNewQtlsPlot("MDGA1~rs6458012")
dev.off()

#' Final question, can I estimate a lower bound on the true FDR in the 57189 things that are purported to be FDR < 0.05 by the conventional model.
#' Under a scenario where the interaction model is doing a better job at identifying borderline significant SNPs this estimate may be conservative, however, under a scenario where the eQTL architecture is reasonably well conserved between tumor infiltraing normal cells and cancer cells, and many of these would inded also be identifed as cancer cells at a larger samples size, this estimate could be too harsh?
sort(adjs_intInv_Peer)[57189] # 0.3820547

57189 * 0.3820547 # [1] 21849.33

#' How many of the top 57189 things are shared between the conventional and interaction models.
sum(names(sort(adjs_intInv_Peer)[1:57189]) %in% names(sort(adjs_conv_Peer)[1:57189])) # [1] 20032

print(sessionInfo())










