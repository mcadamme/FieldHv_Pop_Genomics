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

#getting total counts
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
LCI.obs <- tapply(Tot_data$Counts, Tot_data$Year, boot.fn.lower)
LCI.obs

UCI.obs <- tapply(Tot_data$Counts, Tot_data$Year, boot.fn.upper)
UCI.obs

#getting ddRAD-seq data
RADtag_mean_freq <- tapply(RADtag_data$Counts, RADtag_data$Year, mean)
LCI.obs <- tapply(RADtag_data$Counts, RADtag_data$Year, boot.fn.lower)
LCI.obs

UCI.obs <- tapply(RADtag_data$Counts, RADtag_data$Year, boot.fn.upper)
UCI.obs

Year <- c(1997, 2002, 2007, 2012)

#Figure for pub
library(gplots)

plotCI(x = Year, y = means, uiw = SE, add=T)


#Figure for publication
LA_Freqs <- c(0.66, 0.56, 0.40, 0.44)
TX_Freqs <- c(0.63, NA, 0.36, 0.36)
RAD_Tag <- c(0.73, NA, 0.66, 0.5)
Year <- c(1997, 2002, 2007, 2012)

png(file = "Fig2_PyRes_Allele_Freq.png", units = "px", height = 800, width = 800)
par(mar = c(5,5,4,1))
plot(Year,LA_Freqs, type = "l", ylim = c(0.1, 0.9), xaxt='n', lty = 1, lwd = 2, ylab = "Pyrethroid Resistance Allele Frequency", xlab = "Collection Year", cex.axis = 2, cex.lab = 2.5)
axis(side = 1, at = c("1997", "2002", "2007", "2012"), cex.axis = 2)

lines(Year, approx(TX_Freqs, xout=seq_along(TX_Freqs))$y, lty = 2, lwd = 2)#approximates missing 2002 value
lines(Year, approx(RAD_Tag, xout=seq_along(RAD_Tag))$y, lty = 3, lwd = 2)#approximates missing 2002 value

legend("topright", legend = c("LA", "TX", "RAD-seq Allele"), lty = c(1,2,3), cex = 2)

dev.off()

