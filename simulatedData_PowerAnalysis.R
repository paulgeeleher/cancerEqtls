#' This script does the power analysis on the simulated data (i.e. it creates the supplementary figure with the power analysis (supplementary figure 1)

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

#' Create actual simulated Expression datasets of 1,000 samples and 600 genes....
cancerExpressionSimMat <- numeric(600*1000)
normalExpressionSimMat <- numeric(600*1000)
bulkExpressionSimMat <- numeric(600*1000)
normalEffectsIndependentOfCancer <- numeric(600*1000)
dim(cancerExpressionSimMat) <- c(600, 1000)
dim(normalExpressionSimMat) <- c(600, 1000)
dim(bulkExpressionSimMat) <- c(600, 1000)
dim(normalEffectsIndependentOfCancer) <- c(600, 1000)
numSamps <- 250

# I will use these as the simulated cell proportions, learned from TCGA data.
set.seed(12345) # Setting seed before running all of the below will ensure consistent result each time this is run

# Load the CPE tumor purity estimates for breast cancer.
nComsProps <- read.csv(paste(theRootDir, "ncomms9971-s2.csv", sep=""), as.is=T)
nComsProps_brca <- nComsProps[nComsProps[,2] == "BRCA" & substring(nComsProps[,1], 14, 16) == "01A", "CPE"] # exctract 01A (primary site) breast cancer samples.
nComsProps_brca_noNa <- nComsProps_brca[!is.na(nComsProps_brca)]

# Proportions that I will use in subsequent analysis
theProp <- sample(nComsProps_brca_noNa, 1000) # select 1000 samples at random (without replacement).
propInv <- (1-theProp)

#' Create 100 genes with a mixture of *different* eQTL in both Cancer and normal component. Increment effect sizes from -1 to +1.
simulatedEffectSizesCancer <- seq(-.5, 0.49, .01) # Create 100 simulated effect sizes....
simulatedEffectSizesNormal <- sample(simulatedEffectSizesCancer) # randomeize the above
for(i in 1:100)
{
  # NOTE: rnorm() will add noise and mean of 0 and sd of 1, because the real data were standardized to a mean of 0 and sd of 1.
  cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) + rnorm(numSamps))
  normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) + rnorm(numSamps))
  bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv) # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
}

#' Create 100 genes with an eQTL in cancer only
simulatedEffectSizesCancer[101:200] <- seq(-.5, 0.49, .01)
simulatedEffectSizesNormal[101:200] <- rep(0, 100) 
for(i in 101:200)
{
  cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) + rnorm(numSamps))
  normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) + rnorm(numSamps))
  bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv) # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
  normalEffectsIndependentOfCancer[i,] <- residuals(lm(normalExpressionSimMat[i,]~cancerExpressionSimMat[i,])) # This will be used later, to show the interaction model values arne't correlated with normal.
}

#' Create 100 genes with an eQTL in normal only
simulatedEffectSizesCancer[201:300] <- rep(0, 100) 
simulatedEffectSizesNormal[201:300] <- seq(-.5, 0.49, .01) 
for(i in 201:300)
{
  cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) + rnorm(numSamps))
  normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) + rnorm(numSamps))
  bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv) # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
}

#' Create 100 genes with an eQTL in neither
simulatedEffectSizesCancer[301:400] <- rep(0, 100) 
simulatedEffectSizesNormal[301:400] <- rep(0, 100) 
for(i in 301:400)
{
  cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) + rnorm(numSamps))
  normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) + rnorm(numSamps))
  bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv) # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
}

#' Create 100 genes with the SAME eQTL in Cancer and normal. 
simulatedEffectSizesCancer[401:500] <- seq(-.5, 0.49, .01)
simulatedEffectSizesNormal[401:500] <- seq(-.5, 0.49, .01)
for(i in 401:500)
{
  cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) + rnorm(numSamps))
  normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) + rnorm(numSamps))
  bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv) # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
}

#' Create 100 genes with SIMILAR eQTL in Cancer and Normal. I.e. add some noise to the "cancer" eQTL effect size to create a "normal" effect size.
simulatedEffectSizesCancer[501:600] <- seq(-.5, 0.49, .01)
a <- seq(-.5, 0.49, .01) + rnorm(100, 0, .1)
aScaled <- (a / max(a)) *.5 # Scale this so the effects won't be bigger than .5, to be consitent with everything else.
simulatedEffectSizesNormal[501:600] <- aScaled
for(i in 501:600)
{
  cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) + rnorm(numSamps))
  normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps), rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) + rnorm(numSamps), rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) + rnorm(numSamps))
  bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv) # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
}

#' Lets Now run an analysis and see how many "cancer" and "normal" eQTLs we can recover...
genotype <- c(rep(0, numSamps), rep(1, (numSamps*2)), rep(2, numSamps)) # genotype is always the same
# genotype <- c(rep(1, numSamps), rep(2, (numSamps*2)), rep(3, numSamps)) # This produces the same results as above
# genotype <- c(rep(2, numSamps), rep(3, (numSamps*2)), rep(4, numSamps)) # However, this produces slightly different results....? Not sure if it can be argued that one makes more sense.? For now, lets stick to convention (encoding genotypes as 0,1,2)
cancerEffectInteractionModel <- numeric()
normalEffectInteractionModel <- numeric()
bulkEffectConventionalModel <- numeric()
pValuesBulkTumor <- numeric()
pValuesCancerInteractionModel <- numeric()
interactionTermPvalue <- numeric()
pValuesCancerInteractionModel_randomCpe <- numeric()

#' Create vector of the estimated proportions, but add measurement noise, to reflect the fact that these will not be estimated exactly in the real data.
thePropNoisier <- theProp + rnorm(length(theProp), 0,0.1)

#' Quantile normalize the noise added proportions so they are on an identical distribution to the original proportions (but have noise added, i.e. they will be reordered to some extent).
thePropSort <- sort(theProp)
thePropNoise <- thePropSort[rank(thePropNoisier)]
propInvNoise <- (1-thePropNoise)

#' Plot this correlation (supplementary figure)
thisCor <- cor(thePropNoise, theProp) # this correlation is reasonable based on the correlations achieved by different genomics methods.
dir.create(theRootDir %&% "paper/figures/", recursive=T, showWarnings=F)
svg(file=paste(theRootDir, "paper/figures/suppFig1.svg", sep=""), width=4, height=4)
plot(theProp, thePropNoise, main=paste("Pearson Correlation = ", format(round(thisCor, 2), nsmall = 2)), las=1, cex.axis=.8, pch=20, col="#00000099", xlab="Simulated known proportion", ylab="Simulated measured (noise added) proportion", bty="l")
dev.off()

#' Recover the eQTLs using the different types of models (i.e. interaction (cancer) and conventional (bulk tumor)
for(i in 1:nrow(bulkExpressionSimMat))
{
  cancerEffectInteractionModel[i] <- coef(summary(lm(bulkExpressionSimMat[i,]~genotype*propInvNoise)))[2, 1] # We use this inverse, as we want the main efffect "genotype" to correspond to 0% normal cells (which is 100% cancer.).
  normalEffectInteractionModel[i] <- coef(summary(lm(bulkExpressionSimMat[i,]~genotype*thePropNoise)))[2, 1]
  bulkEffectConventionalModel[i] <- coef(summary(lm(bulkExpressionSimMat[i,]~genotype)))[2, 1]
  pValuesBulkTumor[i] <- coef(summary(lm(bulkExpressionSimMat[i,]~genotype)))[2, 4]
  pValuesCancerInteractionModel[i] <- coef(summary(lm(bulkExpressionSimMat[i,]~genotype*propInvNoise)))[2, 4]
  pValuesCancerInteractionModel_randomCpe[i] <- coef(summary(lm(bulkExpressionSimMat[i,]~genotype*sample(propInvNoise))))[2, 4]
  interactionTermPvalue[i] <- coef(summary(lm(bulkExpressionSimMat[i,]~genotype*propInvNoise)))[4, 4]
}
fdrCancerInteraction <- p.adjust(pValuesCancerInteractionModel, method="BH")
fdrBulkTumor <- p.adjust(pValuesBulkTumor, method="BH")
fdrSigInCancerInteraction <- fdrCancerInteraction < 0.05
fdrSigInBulkTumor <- fdrBulkTumor < 0.05
sum(p.adjust(pValuesCancerInteractionModel_randomCpe, method="BH") < 0.05) # we get the fewest associations with a randomly arranged tumor purity estimate (same as what happens TCGA breast cancer).



#' QQ-plot of the simulated P-values for bulk, cancer model and random CPE model.
convSort <- sort(-log10(pValuesBulkTumor), decreasing=T)
intSort <- sort(-log10(pValuesCancerInteractionModel), decreasing=T)
intSortRand <- sort(-log10(pValuesCancerInteractionModel_randomCpe), decreasing=T)

#' The interaction term itself is doing a poor job of picking out differential eQTLs. At least at only 1,000 samples. I.e. more samples are likely required for this to be useful.
print(sum(p.adjust(interactionTermPvalue, method="BH") < 0.05))

#' Calculate sensitivy and specificity for interaction term
hasSigIntTerm <- which(p.adjust(interactionTermPvalue, method="BH") < 0.05)
isDifferentCancerNormal <- which(simulatedEffectSizesCancer != simulatedEffectSizesNormal)
sum(hasSigIntTerm %in% isDifferentCancerNormal) / length(isDifferentCancerNormal)


#' Calculate Sensitivity across all data points (TP / P)
#' A function to calculate sensitivty, specificy and FDR from a vectors of estimated beta and p-values. Controlling for a specific FDR.
#' I would like this to also return the intex of the false positives, so I can deduced where they originated!
calcSensitivitySpecificityFdr <- function(estimatedBetas, pValues, realBetas, fdr=0.05)
{
  # Sensitivty (true positive rate) = number of true positives / number of real positives. (TP/P)
  # Calculate the number of "positives", i.e. the number of effects really different from zero.
  positives <- which(realBetas != 0)
  numPositives <- sum(realBetas != 0)
  
  # Calculate the number of True Positives, i.e. the number of effects really different from zero that were identified by our model.
  correctDirectionality <- which((estimatedBetas * realBetas) > 0)
  isSignificant <- which(p.adjust(pValues, method="BH") < fdr)
  numTruePositives <- sum((isSignificant %in% correctDirectionality) & (isSignificant %in% positives)) # things significant with correct directionality.
  
  sensitivity <- numTruePositives/numPositives
  
  # Specificity (true negative rate) = the number of true negatives / number of real negatives.
  # the number of real negatives
  negatives <- which(realBetas == 0)
  numNegatives <- sum(realBetas == 0)
  notSignificant <- which(p.adjust(pValues, method="BH") > fdr)
  numTrueNegatives <- sum(notSignificant %in% negatives) # How many negatives were correctly assigned as "not significant"
  specificity <- numTrueNegatives/numNegatives
  
  # Calculate the type I error rate, i.e. the type I error rate....! I.e. the proportion of things called positive that are actually negative.
  MeasuredFalseDiscoveryRate <- sum(isSignificant %in% negatives) / length(isSignificant)
  FalseDiscoveriesIndex <- isSignificant[isSignificant %in% negatives]
  
  return(list(Sensitivty=sensitivity, Specificty=specificity, MeasuredFalseDiscoveryRate=MeasuredFalseDiscoveryRate, numSigEQtls=length(isSignificant), FalseDiscoveriesIndex=FalseDiscoveriesIndex))
}

#' Print the sensitivity and specificity for the cancer and bulk models.
print(calcSensitivitySpecificityFdr(cancerEffectInteractionModel, pValuesCancerInteractionModel, simulatedEffectSizesCancer))
print(calcSensitivitySpecificityFdr(bulkEffectConventionalModel, pValuesBulkTumor, simulatedEffectSizesCancer))


#' What types of eQTL are driving the false positives?
classesOfEqtlsInSimulation <- c("Different eQTLs", "Cancer Only eQTLs", "Normal Only eQTLs", "No eQTLs", "Same eQTL", "Similar eQTL")
bigClassesVec <- as.character(sapply(classesOfEqtlsInSimulation, rep, 100))
falsePositivesByType_tumor <- table(bigClassesVec[calcSensitivitySpecificityFdr(bulkEffectConventionalModel, pValuesBulkTumor, simulatedEffectSizesCancer)$FalseDiscoveriesIndex])
print(falsePositivesByType_tumor)
falsePositivesByType_cancer <- table(bigClassesVec[calcSensitivitySpecificityFdr(cancerEffectInteractionModel, pValuesCancerInteractionModel, simulatedEffectSizesCancer)$FalseDiscoveriesIndex])
print(falsePositivesByType_cancer)


#' Power analysis is performed here by simulation. To estimate the effect size which we have 80% power using conventional / interaction model,
#' I will simulate the change in probability of an eQTL being correctly identified as nominally significant (P < 0.05) as the effect size changes. I will
#' model this change in power using logistic, then use that model to estimate the effect size that can be recovered with 80% power.
#' First do power calculations for "interaction model"
colors <- c("#e41a1c", "#377eb8", "", "", "#4daf4a", "#984ea3")
effectSize_interaction <- numeric()
svg("/mnt/data_scratch/prediXcanProj/Results/Figures/SuppFig_power_interactionModel.svg", width=5, height=3.5)
for(i in c(1,2,5,6)) # Do analysis for "Different eQTLs", "Cancer Only eQTLs", "Same eQTL", "Similar eQTL".....
{
  identified <- as.numeric(pValuesCancerInteractionModel < 0.05)[(1 + (100*(i-1))): (100+(100*(i-1)))]
  thisEffectSize <- abs(simulatedEffectSizesCancer[(1 + (100*(i-1))): (100+(100*(i-1)))])
  dat <- data.frame(identified, thisEffectSize)
  g=glm(identified~thisEffectSize,family=binomial,dat) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
  if(i == 1)
  {
    curve(predict(g,data.frame(thisEffectSize=x),type="resp"),xlim=c(0,.5), ylim=c(0,1), xlab="Effect Size", ylab="Power", col=colors[i], las=1, cex.axis=.8, bty="l") # draws a curve based on prediction from logistic regression model
  }
  else
  {
    curve(predict(g,data.frame(thisEffectSize=x),type="resp"),add=TRUE, xlim=c(0,.5), ylim=c(0,1), xlab="Effect Size", ylab="Power", col=colors[i], las=1, cex.axis=.8) # draws a curve based on prediction from logistic regression model
  }
  points(thisEffectSize,fitted(g),pch=20, col=paste(colors[i], "66", sep="")) # optional: you could skip this draws an invisible set of points of body size survival based on a 'fit' to glm model. pch= changes type of dots.
  effectSize[i] <- (log(.8/(1-.8)) - coef(g)[1])/coef(g)[2] # The formula below will allow me to recover the X value in a logistic regression (here: Effect Size) from the Y value (Power).
}
abline(h=0.8, col="grey", lty=2)
legend("right",  c("Different eQTLs", "Cancer Only eQTLs", "Same eQTL", "Similar eQTL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"), cex=.7)
dev.off()

#' Now do power calculations for "conventional model"
colors <- c("#e41a1c", "#377eb8", "", "", "#4daf4a", "#984ea3")
effectSize_conventional <- numeric()
svg("/mnt/data_scratch/prediXcanProj/Results/Figures/SuppFig_power_conv.svg", width=5, height=3.5)
for(i in c(1,2,5,6)) # Do analysis for "Different eQTLs", "Cancer Only eQTLs", "Same eQTL", "Similar eQTL".....
{
  identified <- as.numeric(pValuesBulkTumor < 0.05)[(1 + (100*(i-1))): (100+(100*(i-1)))]
  thisEffectSize <- abs(simulatedEffectSizesCancer[(1 + (100*(i-1))): (100+(100*(i-1)))])
  dat <- data.frame(identified, thisEffectSize)
  g=glm(identified~thisEffectSize,family=binomial,dat) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
  if(i == 1)
  {
    curve(predict(g,data.frame(thisEffectSize=x),type="resp"),xlim=c(0,.5), ylim=c(0,1), xlab="Effect Size", ylab="Power", col=colors[i], las=1, cex.axis=.8, bty="l") # draws a curve based on prediction from logistic regression model
  }
  else
  {
    curve(predict(g,data.frame(thisEffectSize=x),type="resp"),add=TRUE, xlim=c(0,.5), ylim=c(0,1), xlab="Effect Size", ylab="Power", col=colors[i], las=1, cex.axis=.8) # draws a curve based on prediction from logistic regression model
  }
  points(thisEffectSize,fitted(g),pch=20, col=paste(colors[i], "66", sep="")) # optional: you could skip this draws an invisible set of points of body size survival based on a 'fit' to glm model. pch= changes type of dots.
  effectSize_conventional[i] <- (log(.8/(1-.8)) - coef(g)[1])/coef(g)[2] # The formula below will allow me to recover the X value in a logistic regression (here: Effect Size) from the Y value (Power).
}
abline(h=0.8, col="grey", lty=2)
legend("right",  c("Different eQTLs", "Cancer Only eQTLs", "Same eQTL", "Similar eQTL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"), cex=.7)
dev.off()


