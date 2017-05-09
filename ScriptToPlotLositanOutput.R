#this script plots Lositan output & gets outlier loci sum stats

#setting working directory
setwd("~/temp/bowtie_output/samtoolsANDvcftools_output/LositanReducedDataset_03292017/ConfInt_99/")

#Plotting all 3 by-year comparisons in one panel
plot_ci <- function(f_name, tcolor) {
  cpl <- read.table(f_name, header=FALSE)
  lines(cpl[,1],cpl[,4], lwd=3, type='l', col=tcolor)
}

plot_loci <- function(f_name, color) {
  cpl <- read.table(f_name, header=TRUE, sep='\t')
  points(cpl[,2],cpl[,3], pch=1, col=color)
}

Corner_text <- function(text, location=location){
  legend(location,legend=text, bty ="n", pch=NA, cex=4) 
}

png(file = "SNPoutliersByYearComparison_Lositan.png", units = "px", height = 1200, width = 800)
par(mfrow=c(3,1))
par(mar = c(5,6,4,1))
plot1 <- plot(-10,ylim=c(0,0.4),xlim=c(0.06,0.65), cex.lab=3, cex.axis=2.5, xlab='He', ylab='Fst')
Corner_text(text="A",location= "topleft")
plot_loci('FieldHvirReduced_Fst_1997thru2012.txt', 'black')
plot_ci('FieldHvirReduced_Fst_1997thru2012_confint.txt', 'red')
points(0.522,0.205, pch = 16, cex = 2, col = 'red')
plot2 <- plot(-10,ylim=c(0,0.4),xlim=c(0.06,0.65), cex.lab=3, cex.axis=2.5, xlab='He', ylab='Fst')
Corner_text(text="B",location= "topleft")
plot_loci('FieldHvirReduced_Fst_1997thru2007.txt', 'black')
plot_ci('FieldHvirReduced_Fst_1997thru2007_confint.txt', 'red')
plot3 <- plot(-10,ylim=c(0,0.4),xlim=c(0.06,0.65), cex.lab=3, cex.axis=2.5, xlab='He', ylab='Fst')
Corner_text(text="C",location= "topleft")
plot_loci('FieldHvirReduced_Fst_2007thru2012.txt', 'black')
plot_ci('FieldHvirReduced_Fst_2007thru2012_confint.txt', 'red')

dev.off()

#Getting VennDiagram & list of unique contigs with outliers for all comparisons
library("VennDiagram")

snp_code <- read.table("SNPInfo.txt", head=F)
colnames(snp_code) <- c("Contig","StartPos")

#loading GTF file with structural annotation
GTF <- read.table("H.virescens_OGS1_cleaned.gtf", header = F)
colnames(GTF) <- c("Contig","Prog", "Type", "StartPos", "StopPos", "Score", 
                   "Strand", "Frame", "Attr1", "Attr2", "Attr3", "Attr4", "Attr5", "Attr6", 
                   "Attr7","Attr8","Attr9")

sub_GTF <- GTF[,c(1,4,5,13)]

#Fst Outliers 1997-2012
Hv1997thru2012 <- read.table("FieldHvirReduced_Fst_1997thru2012.txt", header = T)
Hv1997thru2012 <- cbind(Hv1997thru2012,snp_code)
sub_Hv1997thru2012 <- subset(Hv1997thru2012, Pval >= 0.99)
write.table(sub_Hv1997thru2012, file = "SNPOutliers_Hv1997thru2012.txt")
uniqContigs_1997thr2012 <- unique(sub_Hv1997thru2012$Contig)
length(uniqContigs_1997thr2012)

#Fst Outliers 1997-2007
Hv1997thru2007 <- read.table("FieldHvirReduced_Fst_1997thru2007.txt", header = T)
Hv1997thru2007 <- cbind(Hv1997thru2007, snp_code)
sub_Hv1997thru2007 <- subset(Hv1997thru2007, Pval >= 0.99)
write.table(sub_Hv1997thru2007, file = "SNPOutliers_Hv1997thru2007.txt")
uniqContigs_1997thru2007 <- unique(sub_Hv1997thru2007$Contig)
length(uniqContigs_1997thru2007)

#Fst Outliers 2007-2012
Hv2007thru2012 <- read.table("FieldHvirReduced_Fst_2007thru2012.txt", header = T)
Hv2007thru2012 <- cbind(Hv2007thru2012, snp_code)
sub_Hv2007thru2012 <- subset(Hv2007thru2012, Pval >= 0.99)
write.table(sub_Hv2007thru2012, file = "SNPOutliers_Hv2007thru2012.txt")
uniqContigs_2007thru2012 <- unique(sub_Hv2007thru2012$Contig)
length(uniqContigs_2007thru2012)

#Note Heterozygosity values are exceeding 0.5 for some of these SNPs. 
#I understand that this may be due to the way that Lositan calculates heterozygosity.
#In Lositan "Heterozygosity" is the probability of choosing two alleles at random that are different, one from each deme, 
#rather than the prob of picking two alleles at random that are different irrespective of deme. With two demes it can 
#reach 1 when FST is 1. So the reason why it is so high is that we are probably looking at populations that look 
#"inbred" at these loci.

#Getting Venn Diagram of overlapping scaffolds
overlap <- calculate.overlap(x = list( "1997 and 2007" = unique(sub_Hv1997thru2007$Contig), "1997 and 2012" = unique(sub_Hv1997thru2012$Contig),
                                       "2007 and 2012" = unique(sub_Hv2007thru2012$Contig)))

unlist(lapply(overlap, length))

#vector of plotted values comes from above overlap output.  For example, 190 = size of first circle and comes from 120 + 57 + 11 + 2
png(file = "VennDiagram_PlotContigsWithOutlierSNPs.png", units = "px", height = 800, width = 1200)
draw.triple.venn(190,170,33,59,9,13,2, category = c("1997 and 2007", "1997 and 2012", "2007 and 2012"),
                 fill = c("orange", "blue", "green"), cex = 2, cat.cex = 2.5)

dev.off()

#merging significant SNPs with GO terms
merged_sig_1997thru2012 <- merge(sub_Hv1997thru2012, sub_GTF, by = "Contig")
merged_sig_1997thru2007 <- merge(sub_Hv1997thru2007, sub_GTF, by = "Contig")
merged_sig_2007thru2012 <- merge(sub_Hv2007thru2012, sub_GTF, by = "Contig")

#Getting genes within 10kb of SNP outliers
#1997thru2012
subLow_merged_sig_1997thru2012 <- subset(merged_sig_1997thru2012, StartPos.y > merged_sig_1997thru2012$StartPos.x - 10000)
suball_merged_sig_1997thru2012 <- subset(subLow_merged_sig_1997thru2012, StartPos.y < subLow_merged_sig_1997thru2012$StartPos.x + 10000)
write.table(suball_merged_sig_1997thru2012, file = "JamgModels_forSigSNPs_10kb_1997thru2012.txt", row.names = F)

#1997thru2007
subLow_merged_sig_1997thru2007 <- subset(merged_sig_1997thru2007, StartPos.y > merged_sig_1997thru2007$StartPos.x - 10000)
suball_merged_sig_1997thru2007 <- subset(subLow_merged_sig_1997thru2007, StartPos.y < subLow_merged_sig_1997thru2007$StartPos.x + 10000)
write.table(suball_merged_sig_1997thru2007, file = "JamgModels_forSigSNPs_10kb_1997thru2007.txt", row.names = F)

#2007thru2012
subLow_merged_sig_2007thru2012 <- subset(merged_sig_2007thru2012, StartPos.y > merged_sig_2007thru2012$StartPos.x - 10000)
suball_merged_sig_2007thru2012 <- subset(subLow_merged_sig_2007thru2012, StartPos.y < subLow_merged_sig_2007thru2012$StartPos.x + 10000)
write.table(suball_merged_sig_2007thru2012, file = "JamgModels_forSigSNPs_10kb_2007thru2012.txt", row.names = F)

#Getting genes within 20kb of SNP outliers
#1997thru2012
subLow_merged_sig_1997thru2012 <- subset(merged_sig_1997thru2012, StartPos.y > merged_sig_1997thru2012$StartPos.x - 20000)
suball_merged_sig_1997thru2012 <- subset(subLow_merged_sig_1997thru2012, StartPos.y < subLow_merged_sig_1997thru2012$StartPos.x + 20000)
write.table(suball_merged_sig_1997thru2012, file = "JamgModels_forSigSNPs_20kb_1997thru2012.txt", row.names = F)

#1997thru2007
subLow_merged_sig_1997thru2007 <- subset(merged_sig_1997thru2007, StartPos.y > merged_sig_1997thru2007$StartPos.x - 20000)
suball_merged_sig_1997thru2007 <- subset(subLow_merged_sig_1997thru2007, StartPos.y < subLow_merged_sig_1997thru2007$StartPos.x + 20000)
write.table(suball_merged_sig_1997thru2007, file = "JamgModels_forSigSNPs_20kb_1997thru2007.txt", row.names = F)

#2007thru2012
subLow_merged_sig_2007thru2012 <- subset(merged_sig_2007thru2012, StartPos.y > merged_sig_2007thru2012$StartPos.x - 20000)
suball_merged_sig_2007thru2012 <- subset(subLow_merged_sig_2007thru2012, StartPos.y < subLow_merged_sig_2007thru2012$StartPos.x + 20000)
write.table(suball_merged_sig_2007thru2012, file = "JamgModels_forSigSNPs_20kb_2007thru2012.txt", row.names = F)

#Getting genes within 30kb of SNP outliers
#1997thru2012
subLow_merged_sig_1997thru2012 <- subset(merged_sig_1997thru2012, StartPos.y > merged_sig_1997thru2012$StartPos.x - 30000)
suball_merged_sig_1997thru2012 <- subset(subLow_merged_sig_1997thru2012, StartPos.y < subLow_merged_sig_1997thru2012$StartPos.x + 30000)
write.table(suball_merged_sig_1997thru2012, file = "JamgModels_forSigSNPs_30kb_1997thru2012.txt", row.names = F)

#1997thru2007
subLow_merged_sig_1997thru2007 <- subset(merged_sig_1997thru2007, StartPos.y > merged_sig_1997thru2007$StartPos.x - 30000)
suball_merged_sig_1997thru2007 <- subset(subLow_merged_sig_1997thru2007, StartPos.y < subLow_merged_sig_1997thru2007$StartPos.x + 30000)
write.table(suball_merged_sig_1997thru2007, file = "JamgModels_forSigSNPs_30kb_1997thru2007.txt", row.names = F)

#2007thru2012
subLow_merged_sig_2007thru2012 <- subset(merged_sig_2007thru2012, StartPos.y > merged_sig_2007thru2012$StartPos.x - 30000)
suball_merged_sig_2007thru2012 <- subset(subLow_merged_sig_2007thru2012, StartPos.y < subLow_merged_sig_2007thru2012$StartPos.x + 30000)
write.table(suball_merged_sig_2007thru2012, file = "JamgModels_forSigSNPs_30kb_2007thru2012.txt", row.names = F)

#Getting genes within 36kb of SNP outliers
#1997thru2012
subLow_merged_sig_1997thru2012 <- subset(merged_sig_1997thru2012, StartPos.y > merged_sig_1997thru2012$StartPos.x - 36000)
suball_merged_sig_1997thru2012 <- subset(subLow_merged_sig_1997thru2012, StartPos.y < subLow_merged_sig_1997thru2012$StartPos.x + 36000)
write.table(suball_merged_sig_1997thru2012, file = "JamgModels_forSigSNPs_36kb_1997thru2012.txt", row.names = F)

#1997thru2007
subLow_merged_sig_1997thru2007 <- subset(merged_sig_1997thru2007, StartPos.y > merged_sig_1997thru2007$StartPos.x - 36000)
suball_merged_sig_1997thru2007 <- subset(subLow_merged_sig_1997thru2007, StartPos.y < subLow_merged_sig_1997thru2007$StartPos.x + 36000)
write.table(suball_merged_sig_1997thru2007, file = "JamgModels_forSigSNPs_36kb_1997thru2007.txt", row.names = F)

#2007thru2012
subLow_merged_sig_2007thru2012 <- subset(merged_sig_2007thru2012, StartPos.y > merged_sig_2007thru2012$StartPos.x - 36000)
suball_merged_sig_2007thru2012 <- subset(subLow_merged_sig_2007thru2012, StartPos.y < subLow_merged_sig_2007thru2012$StartPos.x + 36000)
write.table(suball_merged_sig_2007thru2012, file = "JamgModels_forSigSNPs_36kb_2007thru2012.txt", row.names = F)
