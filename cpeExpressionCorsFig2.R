#' This script creates Fig2. a-c. I.e. the associations between gene expression and tumor purity.

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

load(file=theRootDir %&% "Results/rDatas/gwasFilteredInputData_WithQNexprAndCpe.RData") # tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin, geneLocation_filt_qnDat, theCpeProps_matched, covariatesVec_filt, expressionMatNorm_filt, expressionMat_filt, genotypesMat_filt, genotypePcs, snpsLocation_filt, geneLocation_filt, 

#' Fig 2(a): The estimated proportion of cancer cells from CPE.
d <- density(theCpeProps_matched, adjust=2, from=0, to=1)
svg(paste(theRootDir, "paper/figures/figure2/fig2a.svg", sep=""), width=3, height=3)
plot(d, col="#e31a1c", las=1, ylim=c(0,3.5), bty="l", xlab="Proportion of cancer cells (CPE)",  main="", cex.axis=.8, panel.first=hist(theCpeProps_matched, bty="l", col="#1f78b4", border="#1f78b4", breaks=25, freq=F, add=T))
dev.off()

#' Fig 2(b): The number of genes significatly correlated with proportions. 
pVals <- numeric()
theCors <- numeric()
for(i in 1:nrow(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin))
{
  corStuff <- cor.test(tpmDatMat_bc_tpm_logged_tumor_keepPatients_qn_norm_filtFinal_fin[i,], theCpeProps_matched, method="pearson")
  pVals[i] <- corStuff$p.value
  theCors[i] <- corStuff$estimate
  print(i)
}

#' Print the numbers of signfiicant correlations
print(sum(p.adjust(pVals, method="BH") < 0.05))
print(length(pVals))

#' What is the cut-off for something to be considered "significant"
corIsSig <- which(p.adjust(pVals, method="BH") < 0.05)
fdrCutoff <- theCors[corIsSig][which.min(abs(theCors[corIsSig]))]

#' Plot the pearson's correlations, this will be figure 2b.
svg(paste(theRootDir, "paper/figures/figure2/fig2b.svg", sep=""), width=3, height=3)
hist(theCors, col="#1f78b4", border="#1f78b4", las=1, main="", cex.axis=.8, freq=F, xlab="Pearson correlation (genes vs CPE)", breaks=100, bty="l")
abline(v=fdrCutoff, col="#e31a1c", lty=2)
abline(v=-fdrCutoff, col="#e31a1c", lty=2)
dev.off()

#' Plot the distribution of the p-values, can be fig2.(a)
svg(file=paste(theRootDir, "paper/figures/figure2/pValsCpeVsExpression.svg", sep=""), width=3, height=3)
hist(pVals, col="#1f78b4", border="#1f78b4",cex.axis=.8, las=1, main="", xlab="P-value", bty="l")
dev.off()


print(sessionInfo())

