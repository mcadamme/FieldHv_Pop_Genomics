#This is the script I used to examine differences in GO distributions.
library(fdrtool)

setwd("~/temp/bowtie_output/samtoolsANDvcftools_output/LositanReducedDataset_03292017/ConfInt_99/For_GO_Analysis/GOAnalysisForR")

#Loading datasets
All_GOs <- read.table("FullDatasetGOAnalysis.txt", header = T)

Go1997thru2007_10kb <- read.table("GODistByLevel_10kb_1997thru2007.txt", header = T)
Go1997thru2012_10kb <- read.table("GODistByLevel_10kb_1997thru2012.txt", header = T)
Go2007thru2012_10kb <- read.table("GODistByLevel_10kb_2007thru2012.txt", header = T)

Go1997thru2007_36kb <- read.table("GODistByLevel_36kb_1997thru2007.txt", header = T)
Go1997thru2012_36kb <- read.table("GODistByLevel_36kb_1997thru2012.txt", header = T)
Go2007thru2012_36kb <- read.table("GODistByLevel_36kb_2007thru2012.txt", header = T)

#merging & formatting datasets for analysis
merged_10kb <- Reduce(function(x, y) merge(x, y, by = "GO_id", all=TRUE), list(All_GOs, Go1997thru2007_10kb, Go1997thru2012_10kb, Go2007thru2012_10kb))
merged_10kb <- merged_10kb[,c(1:4,6,8,10)]
names(merged_10kb) <- c("GO_id", "GO_category", "GO_term", "NumSeqs_All", "NumSeqs_1997thru2007", "NumSeqs_1997thru2012", "NumSeqs_2007thru2012")
merged_10kb[is.na(merged_10kb)] <- 0

merged_36kb <- Reduce(function(x, y) merge(x, y, by = "GO_id", all=TRUE), list(All_GOs, Go1997thru2007_36kb, Go1997thru2012_36kb, Go2007thru2012_36kb))
merged_36kb <- merged_36kb[,c(1:4,6,8,10)]
names(merged_36kb) <- c("GO_id", "GO_category", "GO_term", "NumSeqs_All", "NumSeqs_1997thru2007", "NumSeqs_1997thru2012", "NumSeqs_2007thru2012")
merged_36kb[is.na(merged_36kb)] <- 0


#reshaping to long to get diff counts
reshaped_10kb <- reshape(merged_10kb, direction = "long", varying = list(4:7), v.names = "Counts")
reshaped_10kb <- reshaped_10kb[,c(1:5)]
colnames(reshaped_10kb)[colnames(reshaped_10kb)=="time"] <- "Comparison" #1=All,2=1997thru2007,3=1997thru2012,4=2007thru2012
reshaped_10kb <- reshaped_10kb[with(reshaped_10kb, order(reshaped_10kb[,4], reshaped_10kb[,2])), ]
reshaped_10kb$CatComp <- paste(reshaped_10kb$Comparison, reshaped_10kb$GO_category, sep="_") 
reshaped_10kb$Tot_Counts <- rep(tapply(reshaped_10kb$Counts, reshaped_10kb$CatComp, sum), 
                                times = unlist(tapply(reshaped_10kb$Counts, reshaped_10kb$CatComp, length)), each = 1)
reshaped_10kb$Diff_counts <- reshaped_10kb$Tot_Counts-reshaped_10kb$Counts


reshaped_36kb <- reshape(merged_36kb, direction = "long", varying = list(4:7), v.names = "Counts")
reshaped_36kb <- reshaped_36kb[,c(1:5)]
colnames(reshaped_36kb)[colnames(reshaped_36kb)=="time"] <- "Comparison" #1=All,2=1997thru2007,3=1997thru2012,4=2007thru2012
reshaped_36kb <- reshaped_36kb[with(reshaped_36kb, order(reshaped_36kb[,4], reshaped_36kb[,2])), ]
reshaped_36kb$CatComp <- paste(reshaped_36kb$Comparison, reshaped_36kb$GO_category, sep="_") 
reshaped_36kb$Tot_Counts <- rep(tapply(reshaped_36kb$Counts, reshaped_36kb$CatComp, sum), 
                                times = unlist(tapply(reshaped_36kb$Counts, reshaped_36kb$CatComp, length)), each = 1)
reshaped_36kb$Diff_counts <- reshaped_36kb$Tot_Counts-reshaped_36kb$Counts


#reshaping to wide for easier fisher exact tests.
rereshaped_10kb <- reshape(reshaped_10kb, v.names = c("Counts","Diff_counts"), idvar = "GO_term",
                timevar = "Comparison", direction = "wide")
rereshaped_10kb <- rereshaped_10kb[,c(1:3,6:13)]

rereshaped_36kb <- reshape(reshaped_36kb, v.names = c("Counts","Diff_counts"), idvar = "GO_term",
                           timevar = "Comparison", direction = "wide")
rereshaped_36kb <- rereshaped_36kb[,c(1:3,6:13)]

#Pairwise comparisons
data1 <- rereshaped_10kb[1:nrow(rereshaped_10kb),c(4:7)]#Allwith1997&2007Comparison
data2 <- rereshaped_10kb[1:nrow(rereshaped_10kb),c(4:5,8:9)]#Allwith1997&2012Comparison
data3 <- rereshaped_10kb[1:nrow(rereshaped_10kb),c(4:5,10:11)]#Allwith1997&2012Comparison

data4 <- rereshaped_36kb[1:nrow(rereshaped_36kb),c(4:7)]#Allwith1997&2007Comparison
data5 <- rereshaped_36kb[1:nrow(rereshaped_36kb),c(4:5,8:9)]#Allwith1997&2012Comparison
data6 <- rereshaped_36kb[1:nrow(rereshaped_36kb),c(4:5,10:11)]#Allwith1997&2012Comparison


#looping over pairs of columns for fisher's test, bonferroni-corrected alpha = 0.00093
##inputs and output file name modified for each pairwise comparison.

fish.results <-c()

for (i in 1:nrow(data6)) {
  results <-fisher.test(matrix(unlist(data6[i,c(1:4)]),nrow=2))
  fish.results <- append(fish.results, unlist(results[1]))
  write.table(fish.results, file = "fish.test.2007thru2012_36kb")
}