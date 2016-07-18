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


#Figure for publication
LA_Freqs <- c(0.66, 0.56, 0.40, 0.44)
TX_Freqs <- c(0.63, NA, 0.36, 0.36)
RAD_Tag <- c(0.71, NA, 0.65, 0.5)
Year <- c(1997, 2002, 2007, 2012)

png(file = "Fig2_PyRes_Allele_Freq.png", units = "px", height = 800, width = 800)
par(mar = c(5,5,4,1))
plot(Year,LA_Freqs, type = "l", ylim = c(0.1, 0.9), xaxt='n', lty = 1, lwd = 2, ylab = "Pyrethroid Resistance Allele Frequency", xlab = "Collection Year", cex.axis = 2, cex.lab = 2.5)
axis(side = 1, at = c("1997", "2002", "2007", "2012"), cex.axis = 2)

lines(Year, approx(TX_Freqs, xout=seq_along(TX_Freqs))$y, lty = 2, lwd = 2)#approximates missing 2002 value
lines(Year, approx(RAD_Tag, xout=seq_along(RAD_Tag))$y, lty = 3, lwd = 2)#approximates missing 2002 value

legend("topright", legend = c("LA", "TX", "RAD-seq Allele"), lty = c(1,2,3), cex = 2)

dev.off()

