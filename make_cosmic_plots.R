#' Load cosmis data for RANBP9 and for ERBB2 and make a figure for them.

ranbp <- read.csv("/home/pgeeleher/Downloads/RANBP9_Cnv_detailsTue Dec 19 17_07_00 2017.csv", as.is=T)
erbb2 <- read.csv("/home/pgeeleher/Downloads/ERBB2_Cnv_detailsTue Dec 19 17_06_26 2017.csv", as.is=T)


# test differential expression of ranbp9 in tumor vs matched normal in TCGA BRCA.
theRootDir <- "/mnt/data_scratch/finalData/"
load(file=paste(theRootDir, "dataIn/tenRuvNewStandardApproach.RData", sep="")) # tenRuvNewStandardApproach, cancerTypesVec, 

brcaMat <- tenRuvNewStandardApproach[, cancerTypesVec == "BRCA"]

sampleOriginVecBRCA <- do.call(rbind, strsplit(colnames(brcaMat), ".", fixed=T))[,4]
patientID <- do.call(rbind, strsplit(colnames(brcaMat), ".", fixed=T))[,3]
patientsWithTumAndNormal <- intersect(patientID[sampleOriginVecBRCA == "01A"], patientID[sampleOriginVecBRCA == "11A"])



# filter the lowest 20% of varying genes:
theVars <- numeric()
for(i in 1:nrow(brcaMat))
{
  theVars[i] <- var(brcaMat[i,])
}
toRemove <- order(theVars)[1:(length(theVars)/5)]
brcaMat_filt <- brcaMat[-toRemove,]

# Quantile normalize
library(preprocessCore)
brcaMat_filt_qn <- normalize.quantiles(brcaMat_filt)
rownames(brcaMat_filt_qn) <- rownames(brcaMat_filt)
colnames(brcaMat_filt_qn) <- colnames(brcaMat_filt)

# A not paired t-test
pVals <- numeric()
for(i in 1:nrow(brcaMat_filt_qn))
{
  pVals[i] <- t.test(brcaMat_filt_qn[i, sampleOriginVecBRCA == "01A"], brcaMat_filt_qn[i, sampleOriginVecBRCA == "11A"])$p.value
#   print(i)
}
names(pVals) <- rownames(brcaMat_filt)

sort(pVals)[1:100]
pVals["ERBB2"]
pVals["RANBP9"]



# Find the samples with both tumor and matched normal from the same patient and do paired t-test. This is clearly the best way of doing this as each patient will act as its own control.
brcaMat_filt_qn_01Apats <- brcaMat_filt_qn[, sampleOriginVecBRCA == "01A" & patientID %in% patientsWithTumAndNormal]
brcaMat_filt_qn_11Apats <- brcaMat_filt_qn[, sampleOriginVecBRCA == "11A" & patientID %in% patientsWithTumAndNormal]
colnames(brcaMat_filt_qn_01Apats) <- do.call(rbind, strsplit(colnames(brcaMat_filt_qn_01Apats), ".", fixed=T))[,3]
colnames(brcaMat_filt_qn_11Apats) <- do.call(rbind, strsplit(colnames(brcaMat_filt_qn_11Apats), ".", fixed=T))[,3]
brcaMat_filt_qn_01Apats_ord <- brcaMat_filt_qn_01Apats[, patientsWithTumAndNormal]
brcaMat_filt_qn_11Apats_ord <- brcaMat_filt_qn_11Apats[, patientsWithTumAndNormal]

pValsPaired <- numeric()
meanDiffs <- numeric()
for(i in 1:nrow(brcaMat_filt_qn_01Apats_ord))
{
  pValsPaired[i] <- t.test(brcaMat_filt_qn_01Apats_ord[i, ], brcaMat_filt_qn_11Apats_ord[i, ], paired=T)$p.value
  meanDiffs[i] <- mean(brcaMat_filt_qn_01Apats_ord[i, ]) - mean(brcaMat_filt_qn_11Apats_ord[i, ])
  print(i)
}
names(pValsPaired) <- rownames(brcaMat_filt_qn_01Apats_ord)
pAdjPaired <- p.adjust(pValsPaired, method="BH")
sort(pValsPaired)[1:10]
pValsPaired["ERBB2"]
pValsPaired["RANBP9"]
pAdjPaired["ERBB2"]
pAdjPaired["RANBP9"]


library("ggplot2")
svg("ranbp9_tumor_vs_normal_tcga.svg", width=8, height=8)
d <- data.frame(y = c(brcaMat_filt_qn_01Apats_ord["RANBP9", ], brcaMat_filt_qn_11Apats_ord["RANBP9", ]),
                group = as.factor(rep(c('Tumor', 'Matched Normal'), each = 97)),
                id = rep(1:97, 2))
ggplot(d, aes(y = y)) +
  geom_boxplot(aes(x = rep(c(-3, 3), each = 97), group = group), fill = 'steelblue') +
  geom_point(aes(x = rep(c(-1, 1), each = 97)), size = 3, alpha=.2) +
  geom_line(aes(x = rep(c(-1, 1), each = 97), group = id, alpha=.5)) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),  axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Normalized expression level (TCGA)") +
  xlab("") 
dev.off()

svg("erbb2_tumor_vs_normal_tcga.svg", width=8, height=8)
d <- data.frame(y = c(brcaMat_filt_qn_01Apats_ord["ERBB2", ], brcaMat_filt_qn_11Apats_ord["ERBB2", ]),
                group = as.factor(rep(c('Tumor', 'Matched Normal'), each = 97)),
                id = rep(1:97, 2))
ggplot(d, aes(y = y)) +
  geom_boxplot(aes(x = rep(c(-3, 3), each = 97), group = group), fill = 'steelblue') +
  geom_point(aes(x = rep(c(-1, 1), each = 97)), size = 3, alpha=.2) +
  geom_line(aes(x = rep(c(-1, 1), each = 97), group = id, alpha=.5)) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),  axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Normalized expression level (TCGA)") +
  xlab("") 
dev.off()















