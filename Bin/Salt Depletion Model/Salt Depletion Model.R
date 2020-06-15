#Load necessary libraries
library(zoo)
library(lmtest)
library(MASS)
library(carData)
library(car)
library(lattice)
library(ggplot2)
library(caret)
library(pls)
library(caret)
library(dplyr)

#Setting a variable for the data
#Model - all predictors are standardized
#Model2 - response variable is standardized
#Model3 - zero depletion rates are omitted, i.e. zero kerflow, deprate is not standardized
Model = read.csv("SaltDepRateData3.csv")
DepFit=lm((DepRate)~.,data=Model)
summary(DepFit)
par (mfrow=c(1,2))
plot(fitted(DepFit), resid(DepFit), col = "grey", pch = 20,
     xlab = "Fitted", ylab = "Residuals", main = "Data from Model 1")
abline(h = 0, col = "darkorange", lwd = 2)

par(mfrow=c(1,2))
x=log1p((Model3$DepRate))
y=(Model3$DepRate)
hist(x)
hist(y)

qqnorm(Model3$DepRate)
qqline(Model3$DepRate, col="dodgerblue", lwd="3")
qqnorm(log1p(Model3$DepRate))
qqline(log1p(Model3$DepRate), col="dodgerblue", lwd="3")


