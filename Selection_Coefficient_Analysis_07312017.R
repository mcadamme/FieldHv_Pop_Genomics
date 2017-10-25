#This is the script that I used to examine changes in allele frequencies for significantly diverged SNPs.

#used equations on pg 28 of Falconer and Mackay and solved for s.
library(scatterplot3d)

#function for selection coefficient assuming dominance of p
Dom.Sel.Coef = function (q.1, q.2, g) {
  delta.q = (q.2-q.1)/g
  s = (delta.q)/((q.1^2)*(q.1+delta.q-1))
}

#To plot selection coefficient assuming dominance of p
plot.Dom.Sel.Coef = function(q.1, q.2, g, s, dataframe) {
  delta.q = (q.2-q.1)/g
  df.name <- deparse(substitute(dataframe))
  png(file = print(paste0(df.name,"Dom.png")), units = "px", height = 600, width = 900)
  scatterplot3d(abs(delta.q), q.1, s, highlight.3d = TRUE, col.axis = "blue", 
                cex = 2.5, cex.axis = 1.5, cex.lab = 2, cex.main = 2, col.grid = "lightblue", 
                main = "Dominance of p", angle = 15, xlab = "Delta q", ylab = "", 
                zlab = "Selection Coefficient", pch = 20, zlim = c(0,.8),xlim = c(0,0.025))
 dev.off()
}

#function for selection coefficient assuming incomplete dominance (or no dominance) of p
NoDom.Sel.Coef = function (q.1, q.2, g) {
  delta.q = (q.2-q.1)/g
  s = (delta.q)/(q.1 *(0.5*q.1+delta.q-0.5))
}

#To plot selection coefficient assuming incomplete dominance of p
plot.NoDom.Sel.Coef = function(q.1, q.2, g, s, dataframe) {
  delta.q = (q.2-q.1)/g
  df.name <- deparse(substitute(dataframe))
  png(file = print(paste0(df.name,"NoDom.png")), units = "px", height = 600, width = 900)
  scatterplot3d(abs(delta.q), q.1, abs(s),highlight.3d = TRUE, col.axis = "blue", 
                cex = 2.5, cex.axis = 1.5, angle = 15, cex.lab = 2, cex.main = 2, col.grid = "lightblue", 
                main = "Incomplete Dominance of p", xlab = "Delta q", ylab = "",
                zlab = "Selection Coefficient", pch = 20, zlim = c(0,.8), xlim = c(0,0.025))
  dev.off()
}

#function for selection coefficient assuming recessiveness of p
Rec.Sel.Coef = function (q.1, q.2, g) {
  delta.q = (q.2-q.1)/g
  s = (delta.q)/(((-(1-q.1)^2)*q.1)+(2*q.1*(1-q.1)*delta.q)+(q.1^2*delta.q))
}

#To plot selection coefficient assuming Recessiveness of p
plot.Rec.Sel.Coef = function(q.1, q.2, g, s, dataframe) {
  delta.q = (q.2-q.1)/g
  df.name <- deparse(substitute(dataframe))
  png(file = print(paste0(df.name,"Rec.png")), units = "px", height = 600, width = 900)
  scatterplot3d(abs(delta.q), q.1, abs(s),highlight.3d = TRUE, col.axis = "blue", 
                cex = 2.5, cex.axis = 1.5, angle = 15, cex.lab = 2, cex.main = 2, col.grid = "lightblue", 
                main = "Recessiveness of p", xlab = "Delta q", ylab = "",
                zlab = "Selection Coefficient", pch = 20, zlim = c(0,.8),xlim = c(0,0.025))
  dev.off()
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

########1997thru2012
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
sub1_nozeroes <- data.frame(subset(sub1, Ancest_Freq.1997 > 0.05))
sub2_nozeroes <- data.frame(subset(sub2, Alt_Freq.1997 > 0.05))

Full_SelCoef_Hv1997thru2012 <- rbind(sub1_nozeroes, sub2_nozeroes)
Full_SelCoef_Hv1997thru2012$q.1 <- c(sub1_nozeroes$Alt_Freq.1997, sub2_nozeroes$Ancest_Freq.1997)
Full_SelCoef_Hv1997thru2012$q.2 <- c(sub1_nozeroes$Alt_Freq.2012, sub2_nozeroes$Ancest_Freq.2012)

#starting with assumption of dominance
Full_SelCoef_Hv1997thru2012$dom_sel_coef <- Dom.Sel.Coef(Full_SelCoef_Hv1997thru2012$q.1,Full_SelCoef_Hv1997thru2012$q.2,g)
plot.Dom.Sel.Coef(Full_SelCoef_Hv1997thru2012$q.1,Full_SelCoef_Hv1997thru2012$q.2,g,Full_SelCoef_Hv1997thru2012$dom_sel_coef,Full_SelCoef_Hv1997thru2012)

#Calculating average selection coefficient and sd
mean(Full_SelCoef_Hv1997thru2012$dom_sel_coef)
sd(Full_SelCoef_Hv1997thru2012$dom_sel_coef)
min(Full_SelCoef_Hv1997thru2012$dom_sel_coef)
max(Full_SelCoef_Hv1997thru2012$dom_sel_coef)

#No dominance
Full_SelCoef_Hv1997thru2012$nodom_sel_coef <- NoDom.Sel.Coef(Full_SelCoef_Hv1997thru2012$q.1,Full_SelCoef_Hv1997thru2012$q.2,g)
plot.NoDom.Sel.Coef(Full_SelCoef_Hv1997thru2012$q.1,Full_SelCoef_Hv1997thru2012$q.2,g,Full_SelCoef_Hv1997thru2012$nodom_sel_coef,Full_SelCoef_Hv1997thru2012)

#Calculating average selection coefficient and sd
mean(Full_SelCoef_Hv1997thru2012$nodom_sel_coef)
sd(Full_SelCoef_Hv1997thru2012$nodom_sel_coef)
min(Full_SelCoef_Hv1997thru2012$nodom_sel_coef)
max(Full_SelCoef_Hv1997thru2012$nodom_sel_coef)

#Recessiveness of p
Full_SelCoef_Hv1997thru2012$rec_sel_coef <- Rec.Sel.Coef(Full_SelCoef_Hv1997thru2012$q.1,Full_SelCoef_Hv1997thru2012$q.2,g)
plot.Rec.Sel.Coef(Full_SelCoef_Hv1997thru2012$q.1,Full_SelCoef_Hv1997thru2012$q.2,g,Full_SelCoef_Hv1997thru2012$rec_sel_coef,Full_SelCoef_Hv1997thru2012)

#Calculating average selection coefficient and sd
mean(Full_SelCoef_Hv1997thru2012$rec_sel_coef)
sd(Full_SelCoef_Hv1997thru2012$rec_sel_coef)
min(Full_SelCoef_Hv1997thru2012$rec_sel_coef)
max(Full_SelCoef_Hv1997thru2012$rec_sel_coef)



########1997thru2007
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

sub1_nozeroes <- data.frame(subset(sub1, Ancest_Freq.1997 > 0.05))
sub2_nozeroes <- data.frame(subset(sub2, Alt_Freq.1997 > 0.05))

Full_SelCoef_Hv1997thru2007 <- rbind(sub1_nozeroes, sub2_nozeroes)
Full_SelCoef_Hv1997thru2007$q.1 <- c(sub1_nozeroes$Alt_Freq.1997, sub2_nozeroes$Ancest_Freq.1997)
Full_SelCoef_Hv1997thru2007$q.2 <- c(sub1_nozeroes$Alt_Freq.2007, sub2_nozeroes$Ancest_Freq.2007)

#starting with assumption of dominance
Full_SelCoef_Hv1997thru2007$dom_sel_coef <- Dom.Sel.Coef(Full_SelCoef_Hv1997thru2007$q.1,Full_SelCoef_Hv1997thru2007$q.2,g)
plot.Dom.Sel.Coef(Full_SelCoef_Hv1997thru2007$q.1,Full_SelCoef_Hv1997thru2007$q.2,g,Full_SelCoef_Hv1997thru2007$dom_sel_coef,Full_SelCoef_Hv1997thru2007)

#Calculating average selection coefficient and sd
mean(Full_SelCoef_Hv1997thru2007$dom_sel_coef)
sd(Full_SelCoef_Hv1997thru2007$dom_sel_coef)
min(Full_SelCoef_Hv1997thru2007$dom_sel_coef)
max(Full_SelCoef_Hv1997thru2007$dom_sel_coef)


#No dominance
Full_SelCoef_Hv1997thru2007$nodom_sel_coef <- NoDom.Sel.Coef(Full_SelCoef_Hv1997thru2007$q.1,Full_SelCoef_Hv1997thru2007$q.2,g)
plot.NoDom.Sel.Coef(Full_SelCoef_Hv1997thru2007$q.1,Full_SelCoef_Hv1997thru2007$q.2,g,Full_SelCoef_Hv1997thru2007$nodom_sel_coef,Full_SelCoef_Hv1997thru2007)

##Calculating average selection coefficient and sd
mean(Full_SelCoef_Hv1997thru2007$nodom_sel_coef)
sd(Full_SelCoef_Hv1997thru2007$nodom_sel_coef)
min(Full_SelCoef_Hv1997thru2007$nodom_sel_coef)
max(Full_SelCoef_Hv1997thru2007$nodom_sel_coef)

#Recessiveness of p
Full_SelCoef_Hv1997thru2007$rec_sel_coef <- Rec.Sel.Coef(Full_SelCoef_Hv1997thru2007$q.1,Full_SelCoef_Hv1997thru2007$q.2,g)
plot.Rec.Sel.Coef(Full_SelCoef_Hv1997thru2007$q.1,Full_SelCoef_Hv1997thru2007$q.2,g,Full_SelCoef_Hv1997thru2007$rec_sel_coef,Full_SelCoef_Hv1997thru2007)

#Calculating average selection coefficient and sd
mean(Full_SelCoef_Hv1997thru2007$rec_sel_coef)
sd(Full_SelCoef_Hv1997thru2007$rec_sel_coef)
min(Full_SelCoef_Hv1997thru2007$rec_sel_coef)
max(Full_SelCoef_Hv1997thru2007$rec_sel_coef)


########2007thru2012
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

sub1_nozeroes <- data.frame(subset(sub1, Ancest_Freq.2007 > 0.05))
sub2_nozeroes <- data.frame(subset(sub2, Alt_Freq.2007 > 0.05))

Full_SelCoef_Hv2007thru2012 <- rbind(sub1_nozeroes, sub2_nozeroes)
Full_SelCoef_Hv2007thru2012$q.1 <- c(sub1_nozeroes$Alt_Freq.2007, sub2_nozeroes$Ancest_Freq.2007)
Full_SelCoef_Hv2007thru2012$q.2 <- c(sub1_nozeroes$Alt_Freq.2012, sub2_nozeroes$Ancest_Freq.2012)

#starting with assumption of dominance
Full_SelCoef_Hv2007thru2012$dom_sel_coef <- Dom.Sel.Coef(Full_SelCoef_Hv2007thru2012$q.1,Full_SelCoef_Hv2007thru2012$q.2,g)
plot.Dom.Sel.Coef(Full_SelCoef_Hv2007thru2012$q.1,Full_SelCoef_Hv2007thru2012$q.2,g,Full_SelCoef_Hv2007thru2012$dom_sel_coef,Full_SelCoef_Hv2007thru2012)

#Calculating average selection coefficient and sd
mean(Full_SelCoef_Hv2007thru2012$dom_sel_coef)
sd(Full_SelCoef_Hv2007thru2012$dom_sel_coef)
min(Full_SelCoef_Hv2007thru2012$dom_sel_coef)
max(Full_SelCoef_Hv2007thru2012$dom_sel_coef)

#No dominance
Full_SelCoef_Hv2007thru2012$nodom_sel_coef <- NoDom.Sel.Coef(Full_SelCoef_Hv2007thru2012$q.1,Full_SelCoef_Hv2007thru2012$q.2,g)
plot.NoDom.Sel.Coef(Full_SelCoef_Hv2007thru2012$q.1,Full_SelCoef_Hv2007thru2012$q.2,g,Full_SelCoef_Hv2007thru2012$nodom_sel_coef,Full_SelCoef_Hv2007thru2012)

#Calculating average selection coefficient and sd
mean(Full_SelCoef_Hv2007thru2012$nodom_sel_coef)
sd(Full_SelCoef_Hv2007thru2012$nodom_sel_coef)
min(Full_SelCoef_Hv2007thru2012$nodom_sel_coef)
max(Full_SelCoef_Hv2007thru2012$nodom_sel_coef)

#Recessiveness of p
Full_SelCoef_Hv2007thru2012$rec_sel_coef <- Rec.Sel.Coef(Full_SelCoef_Hv2007thru2012$q.1,Full_SelCoef_Hv2007thru2012$q.2,g)
plot.Rec.Sel.Coef(Full_SelCoef_Hv2007thru2012$q.1,Full_SelCoef_Hv2007thru2012$q.2,g,Full_SelCoef_Hv2007thru2012$rec_sel_coef,Full_SelCoef_Hv2007thru2012)

#Calculating average selection coefficient and sd
mean(Full_SelCoef_Hv2007thru2012$rec_sel_coef)
sd(Full_SelCoef_Hv2007thru2012$rec_sel_coef)
min(Full_SelCoef_Hv2007thru2012$rec_sel_coef)
max(Full_SelCoef_Hv2007thru2012$rec_sel_coef)