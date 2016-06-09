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

anova(model_full, model_red1, test ="Chisq") #Interation term is not important



model_red2 <- glm(y ~ data$EOL + data$Location, binomial) #without year

anova(model_red1, model_red2, test = "Chisq")#Year is very important to predicting allele frequency.
#To look at how important, or get predicted allele frequencies, use plogis.



model_red3 <- glm(y ~ data$Year + data$Location, binomial)#without early or late

anova(model_red1, model_red3, test = "Chisq")#Time of year does not appear to be that important



model_red4 <- glm(y ~ data$Year + data$EOL, binomial)#without location

anova(model_red1, model_red4, test = "Chisq")#Location is also not important.

