#' What do we see if we compare the results fit on the top 10% of the most pure data, with the results from the conventional or interaction models? We would expect the results from the interaction model to be more similar to those obtained from the conventional model.

#' Set root directiory for this project.
source("theRootDir.R")

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")

#' Load the data for the conventional model fit on the top 10% of the data.
load(file=theRootDir %&% "data/forMatrixEQTL/output/theData_conv_10pct.RData") #theData_conv_10pct

#' Load the data for the analyses on the entire dataset.
load(file=theRootDir %&% "data/forMatrixEQTL/output/allTheDatas.RData") # theData_conv_Peer, theData_int_selectPeers, theData_conv_noPeer, theData_int_noPeer, theData_int_Peer, theData_intInv_noPeer, theData_intInv_Peer, theData_intInv_selectPeers

#' Remove the variables we're not interested in.
rm(list=ls()[which(!ls() %in% c("theData_conv_10pct", "theData_conv_Peer", "theData_intInv_Peer", "theRootDir", "%&%"))])


#' Are the significant p-values recovered by the 10pct model more similar to the interaction / conventional model? Answer: almost all of the eQTLs identified in the 10% most pure samples are identified by both the conventional and the interaction model, so there's little we can learn from that.
sigIntInv <- colnames(theData_intInv_Peer$pMatAll)[p.adjust(theData_intInv_Peer$pMatAll[1,], method="BH") < 0.05]
sigConv <- colnames(theData_conv_Peer$pMatAll)[p.adjust(theData_conv_Peer$pMatAll[1,], method="BH") < 0.05]
sig10pct <- colnames(theData_conv_10pct$pMatAll)[p.adjust(theData_conv_10pct$pMatAll[1,], method="BH") < 0.05]
sum(sigIntInv %in% sig10pct) # [1] 231
length(sig10pct) # [1] 248
sum(sigConv %in% sig10pct) # [1] 233
length(sigConv) # [1] 57189
length(sigIntInv) # [1] 8833

#' Associations tested were the same and are in the smae order? yes.
sum(!colnames(theData_intInv_Peer$betaMatAll) == colnames(theData_conv_10pct$betaMatAll[1,]))

#' Correlations of betas from convetional, interaction models against 10pct model.
print(cor.test(theData_intInv_Peer$betaMatAll[1,], theData_conv_10pct$betaMatAll[1,], method="spearman")) # 0.4265433
print(cor.test(theData_conv_Peer$betaMatAll, theData_conv_10pct$betaMatAll, method="spearman")) # 0.2800685

print(cor.test(theData_intInv_Peer$betaMatAll[1,], theData_conv_10pct$betaMatAll[1,], method="pearson")) # 0.4472866 CI: 0.4464601 0.4481123
print(cor.test(theData_conv_Peer$betaMatAll[1,], theData_conv_10pct$betaMatAll[1,], method="pearson")) # 0.2996332 CI: 0.2986930 0.3005729

summary(lm(theData_conv_10pct$betaMatAll[1,]~theData_intInv_Peer$betaMatAll[1,]))
summary(lm(theData_conv_10pct$betaMatAll[1,]~theData_conv_Peer$betaMatAll[1,]))


#' Plot these associations (use the hexbin feature in ggplot? these don't look great)
# aDf <- data.frame("Interaction Model Betas"=theData_intInv_Peer$betaMatAll[1,], "Top 10% Model Betas"=theData_conv_10pct$betaMatAll[1,])
# library(ggplot2)
# ggplot(aDf, aes(Interaction.Model.Betas, Top.10..Model.Betas)) + geom_hex()  + scale_y_continuous(trans='log10')  + scale_x_continuous(trans='log10')

#' Try plotting randomly selected points in the base plotting package? These look better.
set.seed(12345)
randPoints <- sample(1:length(theData_intInv_Peer$betaMatAll[1,]), 1000)
svg(file="/home/ubuntu/Dropbox/predixcanProj/paper/figures/suppfigs/Suppfig_10pct.svg", width=6, height=3)
par(mfrow=c(1,2))
plot(theData_intInv_Peer$betaMatAll[1,randPoints], theData_conv_10pct$betaMatAll[1,randPoints], col="#00000022", pch=20, bty="l", xlab="Effect size\nInteraction model", ylab="Effect size 10% model", las=1, cex.axis=.8)
int10pctlm <- lm(theData_conv_10pct$betaMatAll[1,randPoints]~theData_intInv_Peer$betaMatAll[1,randPoints])
abline(coef(int10pctlm), col="red")
plot(theData_conv_Peer$betaMatAll[1,randPoints], theData_conv_10pct$betaMatAll[1,randPoints], col="#00000022", pch=20, bty="l", xlab="Effect size\nConventional model", ylab="Effect size 10% model", las=1, cex.axis=.8)
conv10pctlm <- lm(theData_conv_10pct$betaMatAll[1,randPoints]~theData_conv_Peer$betaMatAll[1,randPoints])
abline(coef(conv10pctlm), col="red")
dev.off()
summary(conv10pctlm)
summary(int10pctlm)







