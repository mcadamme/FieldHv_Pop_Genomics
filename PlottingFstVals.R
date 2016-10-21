#R script to plot weir and cockerham fst per loc

setwd("~/temp/bowtie_output/samtoolsANDvcftools_output")

data <- read.table("data.goes.here.weir.fst", header = T)
head(data)

data$CHROM <- gsub('Contig', '', data$CHROM)
data$order <- seq(1:nrow(data))#giving SNPs a unique position in the ordered dataset.

sub_data <- subset(data, WEIR_AND_COCKERHAM_FST > 0.4)#getting contigs with highly divergent SNPs 

png(file = "Fig5_2007and2012WCFst.png", units = "px", height = 800 , width = 1200)
par(mar=c(5,5,4,1))
plot(data$order, data$WEIR_AND_COCKERHAM_FST, ylim = c(-0.05,0.7), ylab = "Weir and Cockerham's Fst", cex.lab = 2, cex.axis = 2, pch = 16, xlab = "SNP position",xaxt='n')
  #points(36453, 0.2042, col = "black", bg = "red", pch = 21, cex = 1.8)#Adding SNP 11655 linked to NaNhp.
  #points(36457, 0.2101, col = "black", bg = "red", pch = 21, cex = 1.8)#Adding SNP 11706 linked to NaNhp.
  abline(0.004, 0, col = "grey", lty = 2, lwd = 2)
dev.off()
