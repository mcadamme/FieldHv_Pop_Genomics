#This is the script that I used to examine changes in allele frequencies for significantly diverged SNPs.

#used equations on pg 28 of Falconer and Mackay and solved for s.
library(scatterplot3d)

#function for selection coefficient assuming dominance of p
Dom.Sel.Coef = function (q.1, q.2, g) {
  delta.q = (q.1-q.2)/g
  s = (delta.q)/((q.1^2)*(q.1+delta.q-1))
  dom.plot <- scatterplot3d(delta.q,q.1,abs(s),highlight.3d = TRUE, col.axis = "blue",
                            col.grid = "lightblue", pch = 20)
 print(dom.plot)
}

#function for selection coefficint assuming no dominance of p
NoDom.Sel.Coef = function (q.1, q.2, g) {
  delta.q = (q.1-q.2)/g
  s = (delta.q)/(q.1 *(0.5*q.1+delta.q-0.5))
}

setwd("~/temp/bowtie_output/samtoolsANDvcftools_output/LositanReducedDataset_03292017/Allele_Frequencies_by_year")

#reading in allele frequency datasets
Freqs1997 <- read.csv("FieldHvir1997_freq.csv", header = T)
Freqs2012 <- read.csv("FieldHvir2012_freq.csv", header = T)
Freqs2007 <- read.csv("FieldHvir2007_freq.csv", header = T)

#reading in subsets of signficant outliers
sub_Hv1997thru2012 <- read.table("SNPOutliers_Hv1997thru2012.txt", header = T)
sub_Hv1997thru2007 <- read.table("SNPOutliers_Hv1997thru2007.txt", header = T)
sub_Hv2007thru2012 <- read.table("SNPOutliers_Hv2007thru2012.txt", header = T)

#Pulling out frequencies for divergent loci
#1997thru2012
merged_Hv1997 <- merge(sub_Hv1997thru2012, Freqs1997, by = c("Contig","StartPos"))
merged_Hv1997$Year <- rep("1997", times = nrow(merged_Hv1997))
merged_Hv2012 <- merge(sub_Hv1997thru2012, Freqs2012, by = c("Contig","StartPos"))
merged_Hv2012$Year <- rep("2012", times = nrow(merged_Hv2012))

Full_Hv1997thru2012 <- rbind(merged_Hv1997,merged_Hv2012)

#Reshaping from long to wide to get delta q
reshaped_Hv1997thru2012 <- reshape(Full_Hv1997thru2012, v.names = c("Ancest_Freq","Alt_Freq",
                          "N_CHR"), idvar = "Locus", timevar = "Year", direction = "wide")

#getting decreasing allele for q
reshaped_Hv1997thru2012$Diff_Ancest_Freq <- reshaped_Hv1997thru2012$Ancest_Freq.2012-reshaped_Hv1997thru2012$Ancest_Freq.1997
reshaped_Hv1997thru2012$Diff_Alt_Freq <- reshaped_Hv1997thru2012$Alt_Freq.2012-reshaped_Hv1997thru2012$Alt_Freq.1997

sub1 <- subset(reshaped_Hv1997thru2012, Diff_Ancest_Freq > 0)
sub2 <- subset(reshaped_Hv1997thru2012, Diff_Alt_Freq > 0)

g = (2012-1997)*4 #60 generations

#removing the few situations where starting allele frequencies are zero, 
#which is unrealistic and inflates s.
sub1_nozeroes <- subset(sub1, Ancest_Freq.1997 > 0)#removes 2 samples
sub2_nozeroes <- subset(sub2, Alt_Freq.1997 > 0)#removes 12 samples

#starting with assumption of dominance
sub1_nozeroes$dom_sel_coef_1997thru2012 <- Dom.Sel.Coef(sub1_nozeroes$Alt_Freq.1997,sub1_nozeroes$Alt_Freq.2012,g)
sub2_nozeroes$dom_sel_coef_1997thru2012 <- Dom.Sel.Coef(sub2_nozeroes$Ancest_Freq.1997,sub2_nozeroes$Ancest_Freq.2012,g)

#no dominance
sub1_nozeroes$nd_sel_coef_1997thru2012 <- NoDom.Sel.Coef(sub1_nozeroes$Alt_Freq.1997,sub1_nozeroes$Alt_Freq.2012,g)
sub2_nozeroes$nd_sel_coef_1997thru2012 <- NoDom.Sel.Coef(sub2_nozeroes$Ancest_Freq.1997,sub2_nozeroes$Ancest_Freq.2012,g)


Full_SelCoef_Hv1997thru2012 <- rbind(sub1_nozeroes, sub2_nozeroes)

#getting a vector of q.1
Full_SelCoef_Hv1997thru2012$q.1_1997thru2012 <- c(sub1_nozeroes$Alt_Freq.1997,sub2_nozeroes$Ancest_Freq.1997)


#1997thru2007
merged_Hv1997 <- merge(sub_Hv1997thru2007, Freqs1997, by = c("Contig","StartPos"))
merged_Hv1997$Year <- rep("1997", times = nrow(merged_Hv1997))
merged_Hv2007 <- merge(sub_Hv1997thru2007, Freqs2007, by = c("Contig","StartPos"))
merged_Hv2007$Year <- rep("2007", times = nrow(merged_Hv2007))

Full_Hv1997thru2007 <- rbind(merged_Hv1997,merged_Hv2007)

#Reshaping from long to wide to get delta q
reshaped_Hv1997thru2007 <- reshape(Full_Hv1997thru2007, v.names = c("Ancest_Freq","Alt_Freq",
                                                                    "N_CHR"), idvar = "Locus", timevar = "Year", direction = "wide")

#getting decreasing allele for q
reshaped_Hv1997thru2007$Diff_Ancest_Freq <- reshaped_Hv1997thru2007$Ancest_Freq.2007-reshaped_Hv1997thru2007$Ancest_Freq.1997
reshaped_Hv1997thru2007$Diff_Alt_Freq <- reshaped_Hv1997thru2007$Alt_Freq.2007-reshaped_Hv1997thru2007$Alt_Freq.1997

sub1 <- subset(reshaped_Hv1997thru2007, Diff_Ancest_Freq > 0)
sub2 <- subset(reshaped_Hv1997thru2007, Diff_Alt_Freq > 0)

g = (2007-1997)*4 #40 generations

#starting with assumption of dominance
sub1$dom_sel_coef_1997thru2007 <- Dom.Sel.Coef(sub1$Alt_Freq.1997,sub1$Alt_Freq.2007,g)
sub2$dom_sel_coef_1997thru2007 <- Dom.Sel.Coef(sub2$Ancest_Freq.1997,sub2$Ancest_Freq.2007,g)

#no dominance
sub1$nd_sel_coef_1997thru2007 <- NoDom.Sel.Coef(sub1$Alt_Freq.1997,sub1$Alt_Freq.2007,g)
sub2$nd_sel_coef_1997thru2007 <- NoDom.Sel.Coef(sub2$Ancest_Freq.1997,sub2$Ancest_Freq.2007,g)

#removing the few situations where starting allele frequencies are zero, 
#which is unrealistic and inflates s.
sub1_nozeroes <- subset(sub1, Ancest_Freq.1997 > 0)#removes 1 sample
sub2_nozeroes <- subset(sub2, Alt_Freq.1997 > 0)#removes 11 samples

Full_SelCoef_Hv1997thru2007 <- rbind(sub1_nozeroes, sub2_nozeroes)

#2007thru2012
merged_Hv2007 <- merge(sub_Hv2007thru2012, Freqs2007, by = c("Contig","StartPos"))
merged_Hv2007$Year <- rep("2007", times = nrow(merged_Hv2007))
merged_Hv2012 <- merge(sub_Hv2007thru2012, Freqs2012, by = c("Contig","StartPos"))
merged_Hv2012$Year <- rep("2012", times = nrow(merged_Hv2012))

Full_Hv2007thru2012 <- rbind(merged_Hv2007,merged_Hv2012)

#Reshaping from long to wide to get delta q
reshaped_Hv2007thru2012 <- reshape(Full_Hv2007thru2012, v.names = c("Ancest_Freq","Alt_Freq",
                                                                    "N_CHR"), idvar = "Locus", timevar = "Year", direction = "wide")

#getting decreasing allele for q
reshaped_Hv2007thru2012$Diff_Ancest_Freq <- reshaped_Hv2007thru2012$Ancest_Freq.2012-reshaped_Hv2007thru2012$Ancest_Freq.2007
reshaped_Hv2007thru2012$Diff_Alt_Freq <- reshaped_Hv2007thru2012$Alt_Freq.2012-reshaped_Hv2007thru2012$Alt_Freq.2007

sub1 <- subset(reshaped_Hv2007thru2012, Diff_Ancest_Freq > 0)
sub2 <- subset(reshaped_Hv2007thru2012, Diff_Alt_Freq > 0)

g = (2012-2007)*4 #20 generations

#starting with assumption of dominance
sub1$dom_sel_coef_2007thru2012 <- Dom.Sel.Coef(sub1$Alt_Freq.2007,sub1$Alt_Freq.2012,g)
sub2$dom_sel_coef_2007thru2012 <- Dom.Sel.Coef(sub2$Ancest_Freq.2007,sub2$Ancest_Freq.2012,g)

#no dominance
sub1$nd_sel_coef_2007thru2012 <- NoDom.Sel.Coef(sub1$Alt_Freq.2007,sub1$Alt_Freq.2012,g)
sub2$nd_sel_coef_2007thru2012 <- NoDom.Sel.Coef(sub2$Ancest_Freq.2007,sub2$Ancest_Freq.2012,g)

#removing the few situations where starting allele frequencies are zero, 
#which is unrealistic and inflates s.
sub1_nozeroes <- subset(sub1, Ancest_Freq.2007 > 0)#removes 0 samples
sub2_nozeroes <- subset(sub2, Alt_Freq.2007 > 0)#removes 1 sample

Full_SelCoef_Hv2007thru2012 <- rbind(sub1_nozeroes, sub2_nozeroes)



