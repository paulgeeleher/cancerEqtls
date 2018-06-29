cpe <- read.csv("/home/pgeeleher/Dropbox/predixcanProj/data/tcgaDeconvolution/ncomms9971-s2.csv")

bymedian <- with(cpe, reorder(Cancer.type, -CPE, median, na.rm = TRUE))

cols <- colorRampPalette(c('#1f78b4', '#e31a1c'))(21)

svg("boxplotCpe.svg", width=6, height=3)
boxplot(CPE~bymedian, cpe, las=2, bty="l", cex.axis=.9, pch=20, col=cols, lty=1, outcol="#00000066") 
dev.off()

