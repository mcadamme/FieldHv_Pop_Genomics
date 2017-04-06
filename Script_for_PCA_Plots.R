#This is the script I used to produce my PCA to summarize genetic diversity among my populations.

library(adegenet)

setwd('~/temp/bowtie_output/samtoolsANDvcftools_output/LositanReducedDataset_03292017')

FieldHvirData <- read.genepop("Filtered_FieldHvir03292017_ALL.gen", missing="mean")

data_sum <- summary(FieldHvirData)

#Replacing missing data and tranforming allele freqs for PCA

#Trans_data1 <- genind2genpop(FieldHvirData)

#Trans_data2 <- scaleGen(Trans_data1)
Trans_data1 <- scaleGen(FieldHvirData)

#computing PCA
pca1 <- dudi.pca(Trans_data1, scale=FALSE, scannf=TRUE)

#getting pca plots for all 4 axes
s.class(pca1$li, pop(FieldHvirData), xax=1, yax=2, sub="PCA 1-2") #These demonstrate that populations are not well differentiated.
s.class(pca1$li, pop(FieldHvirData), xax=1, yax=3, sub="PCA 1-3")
s.class(pca1$li, pop(FieldHvirData), xax=1, yax=4, sub="PCA 1-4")

#A publishable version
col <- funky(6)
s.class(pca1$li, pop(FieldHvirData), xax=1, yax=2, col=transp(col, .6), label= NULL, axesell=F, cstar=0, cpoint=3, grid=F)

s.class(pca1$li, pop(FieldHvirData), xax=1, yax=3, col=transp(col, .6), label= NULL, axesell=F, cstar=0, cpoint=3, grid=F)

