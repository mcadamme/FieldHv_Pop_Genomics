#This is my code to analyze my pyrethroid resistance allele frequencies over time.

setwd("~/Downloads") #tells R where to find the file - you may have to change this.

data <- read.csv("NaNhp_data_2.2.16.csv", header = T)
head(data) #Is R interpreting my datafile like I think it should?
str(data)

data$Year <- as.factor(data$Year)

y <- cbind(data$X.A, data$Total-data$X.A) #making a matrix of successes and failures.
head(y)

model_full <- glm(y ~ data$Year * data$EOL + data$Location, binomial)
summary(model_full)

model_red1 <- glm(y ~ data$Year + data$EOL + data$Location, binomial)

anova(model_full, model_red1, test ="Chisq")


model_red2 <- glm(y ~ data$EOL + data$Location, binomial) #without year

anova(model_red1, model_red2, test = "Chisq")#Year is very important to predicting allele frequency.
#To look at how important, or get predicted allele frequencies, use plogis.



model_red3 <- glm(y ~ data$Year + data$Location, binomial)#without early or late

anova(model_red1, model_red3, test = "Chisq")



model_red4 <- glm(y ~ data$Year + data$EOL, binomial)#without location

anova(model_red1, model_red4, test = "Chisq")


#A better data format for plotting - LA first
Counts <- c(rep(1, times = (94+35)), rep(0, times = (40+25)), rep(1, times = (58+56)), rep(0, times = (38+52)),
  rep(1, times = (38+68)), rep(0, times = (66+96)), rep(1, times = (25+60)), rep(0, times = (25+84)))
Year <- c(rep(1997,times = (94+35+40+25)), rep(2002, times = (58+56+38+52)), rep(2007, times = (38+68+66+96)), 
  rep(2012, times = (25+25+60+84)))
Loc <- rep("LA", times = length(Counts))

LA_data <- data.frame(cbind(Counts, Year, Loc))
LA_data$Counts <- as.numeric(as.character(LA_data$Counts))

#now for TX data
Counts <- c(rep(1, times = (7+83)), rep(0, times = (1+51)), rep(1, times = (43)), rep(0, times = (77)), rep(1, times = (7+63)), 
  rep(0, times = (15+111)))
Year <- c(rep(1997,times = (7+83+1+51)), rep(2007, times = (43+77)), rep(2012, times = (7+63+15+111)))
Loc <- rep("TX", times = length(Counts))

TX_data <- data.frame(cbind(Counts, Year, Loc))
TX_data$Counts <- as.numeric(as.character(TX_data$Counts))

#getting total counts for pyrethroid resistance allele
Tot_data <- data.frame(rbind(LA_data, TX_data))


#Finally RADtag data 
Counts <- c(rep(1, times = 61), rep(0, times = 23), rep(1, times = 62), rep(0, times = 31), rep(1, times = 41),
  rep(0, times = 41))
Year <- c(rep(1997, times = (61+23)), rep(2007, times = (62+31)), rep(2012, times = (41+41)))
Loc <- rep("RADtag", times = length(Counts))

RADtag_data <- data.frame(cbind(Counts, Year, Loc))
RADtag_data$Counts <- as.numeric(as.character(RADtag_data$Counts))



#getting means and confidence intervals

boot.fn.lower <- function(x=mean, N=5000) {
  lower.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  lower.CI <- quantile(lower.1, probs=c(0.025))
  lower.CI
}

boot.fn.upper <- function(x=mean, N=5000) {
  upper.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  upper.CI <- quantile(upper.1, probs=c(0.975))
  upper.CI
}

#getting Na channel gene res allele freq
Tot_mean_freq <- tapply(Tot_data$Counts, Tot_data$Year, mean)
Tot_LCI.obs <- tapply(Tot_data$Counts, Tot_data$Year, boot.fn.lower)
Tot_UCI.obs <- tapply(Tot_data$Counts, Tot_data$Year, boot.fn.upper)


#getting ddRAD-seq data
RADtag_mean_freq <- tapply(RADtag_data$Counts, RADtag_data$Year, mean)
LCI.obs <- tapply(RADtag_data$Counts, RADtag_data$Year, boot.fn.lower)
LCI.obs

UCI.obs <- tapply(RADtag_data$Counts, RADtag_data$Year, boot.fn.upper)
UCI.obs


Year <- c(1997, 2002, 2007, 2012, 1997,2007,2012)
All_means <- c(Tot_mean_freq, RADtag_mean_freq)
Marker <- c("SNP","SNP","SNP","SNP","RADtag","RADtag","RADtag")
UCI <- c(Tot_UCI.obs,UCI.obs)
LCI <- c(Tot_LCI.obs,LCI.obs)

All_data <- data.frame(cbind(Year,All_means,Marker,UCI,LCI))
All_data$Year <- as.numeric(as.character(All_data$Year))
All_data$All_means <- as.numeric(as.character(All_data$All_means))
All_data$UCI <- as.numeric(as.character(All_data$UCI))
All_data$LCI <- as.numeric(as.character(All_data$LCI))


#Figure2 for pub
library(ggplot2)

png(file = "Fig2_PyRes_Allele_Freq.png", units = "px", height = 800, width = 1200)
par(mar = c(5,5,4,1))
ggplot(data = All_data, aes(x = Year, y = All_means, ymin = LCI, ymax = UCI, colour = Marker)) +
    geom_line(aes(linetype=Marker), cex=2) +
    geom_errorbar(position = position_dodge(width = 0.4), width = 0.3) +
    scale_colour_manual(values = c("black", "red")) + ylab("Allele Frequency") +
    theme_bw(base_size = 18) +
    theme(panel.grid.major.y = element_line(colour = "white", linetype = "dashed"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(), legend.position="none", legend.key = element_blank(), axis.text=element_text(size=24, face= "bold"), axis.title = element_text(size=28, face= "bold"), axis.title.y = element_text(vjust = 0.3)) + scale_x_continuous(limit=c(1996,2014), breaks = seq(1996, 2014, 3)) + ylim(0.2,0.9)
dev.off()

#Test for differences between slopes - resistance allele and hv_11322_hap1

#first combining datasets
Tot_data$MarkerType <- rep("PYR", times = nrow(Tot_data))
RADtag_data$MarkerType <- rep("RAD", times = nrow(RADtag_data))

RADpPYR <- rbind(Tot_data, RADtag_data)
RADpPYR$Year <- as.numeric(as.character(RADpPYR$Year))

str(RADpPYR)

#model 1 with interaction term
mod1 <- glm(RADpPYR$Counts ~ RADpPYR$Year + RADpPYR$MarkerType + RADpPYR$Year*RADpPYR$MarkerType, binomial)

#model 2 without interaction term
mod2 <- glm(RADpPYR$Counts ~ RADpPYR$Year + RADpPYR$MarkerType, binomial)

#test for significance
anova(mod1, mod2, test = "Chisq")


