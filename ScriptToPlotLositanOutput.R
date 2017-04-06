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

png('output.png', units = "px", height = 1200, width = 800)
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

#for all comparisons

snp_code <- read.table("SNPInfo.txt", head=F)

#Fst Outliers 1997-2012
Hv1997thru2012 <- read.table("FieldHvirReduced_Fst_1997thru2012.txt", header = T)
Hv1997thru2012 <- cbind(Hv1997thru2012,snp_code)
sub_Hv1997thru2012 <- subset(Hv1997thru2012, Pval >= 0.99)
uniqContigs_1997thr2012 <- unique(sub_Hv1997thru2012$V1)
uniqContigs_1997thr2012

#Fst Outliers 1997-2007
Hv1997thru2007 <- read.table("FieldHvirReduced_Fst_1997thru2007.txt", header = T)
Hv1997thru2007 <- cbind(Hv1997thru2007, snp_code)
sub_Hv1997thru2007 <- subset(Hv1997thru2007, Pval >= 0.99)
uniqContigs_1997thru2007 <- unique(sub_Hv1997thru2007$V1)
uniqContigs_1997thru2007

#Fst Outliers 2007-2012
Hv2007thru2012 <- read.table("FieldHvirReduced_Fst_2007thru2012.txt", header = T)
Hv2007thru2012 <- cbind(Hv2007thru2012, snp_code)
sub_Hv2007thru2012 <- subset(Hv2007thru2012, Pval >= 0.99)
uniqContigs_2007thru2012 <- unique(sub_Hv2007thru2012$V1)
uniqContigs_2007thru2012

#Getting Venn Diagram of overlapping scaffolds
overlap <- calculate.overlap(x = list( "1997 and 2007" = sub_Hv1997thru2007$V1, "1997 and 2012" = sub_Hv1997thru2012$V1,
                                       "2007 and 2012" = sub_Hv2007thru2012$V1))

unlist(lapply(overlap, length))

draw.triple.venn(201,184,35,75,15,15,6, category = c("1997 and 2007", "1997 and 2012", "2007 and 2012"))

