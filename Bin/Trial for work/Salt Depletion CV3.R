#Load necessary libraries
library(zoo)
library(lmtest)
library(MASS)
library(carData)
library(car)
library(lattice)
library(ggplot2)
library(caret)
library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)

#OLS
Model = read.csv("DepModelCV6Final.csv")
pairs(Model, col="dodgerblue")
round(cor(Model),2)
DepFit=lm(DepRate~.,data=Model)
summary(DepFit)
par (mfrow=c(2,2),mar=c(2,2,2,2))
plot(DepFit)


#outliers
par(mfrow = c(1,2))
plot(hatvalues(DepFit))
abline(h=15/1596, col="red")
plot(rownames(Model),rstudent(DepFit))
abline(h=3,col="red")


#Stepwise OLS
ModelStep=step(DepFit)
DepFitReduced = lm(DepRate~KerFeed + KerTemp + WaterLev + WaterpH + ECPLev + ReactorIP + 
                     ReactordP + ClaydP + ClayOP, data=Model) 
summary(DepFitReduced)
par (mfrow=c(2,2),mar=c(2,2,2,2))
plot(DepFitReduced)

#outliers
par(mfrow = c(1,2))
plot(hatvalues(DepFitReduced))
abline(h=15/1596, col="red")
plot(rownames(Model),rstudent(DepFitReduced))
abline(h=3,col="red")

#Log Transformation of Response Variable
DepFitTrans=lm(log1p(DepRate)~., data=Model)
summary(DepFitTrans)
par (mfrow=c(2,2),mar=c(2,2,2,2))
plot(DepFitTrans)

#Stepwise Log Transformation of Response Variable
ModelStep2=step(DepFitTrans)
DepFitReduced2=lm(log1p(DepRate) ~ (KerFeed) + (KerTemp) + (WaterLev) +
                    (BrineLev) + (WaterpH) + (WCBootLev) + (ECPdP) + 
                    (ECPLev) + (ReactorIP) + (ReactordP) + (ReactorCauLev) + 
                    (ClaydP) + (ClayOP), data = Model)
summary(DepFitReduced2)
DepFitReduced2$coefficients
par (mfrow=c(2,2),mar=c(2,2,2,2))
plot(DepFitReduced2)

#outliers
par(mfrow = c(1,2))
plot(hatvalues(DepFitReduced2))
abline(h=15/1596, col="red")
plot(rownames(Model),rstudent(DepFitReduced2))
abline(h=3,col="red")

#BoxCox Transformation
boxcox(DepFit, lambda = seq(0.5, 1, by=0.025), plotit=TRUE)
#lambda = 2/3
y=2/3
DepFitBoxCox = lm((((DepRate^y)-1)/y)~., data=Model)
summary(DepFitBoxCox)
par (mfrow=c(2,2),mar=c(2,2,2,2))
plot(DepFitBoxCox)

#Stepwise BoxCox Transformation
ModelStep3=step(DepFitBoxCox)
DepFitReduced3 = lm((((DepRate^y)-1)/y)~KerFeed + KerTemp + WaterLev + BrineLev + 
                      WaterpH + WCBootLev + ECPdP + ECPLev + ReactorIP + ClaydP + 
                      ClayOP, data=Model)
summary(DepFitReduced3)
par (mfrow=c(2,2),mar=c(2,2,2,2))
plot(DepFitReduced3)

#outliers
par(mfrow = c(1,2))
plot(hatvalues(DepFitReduced3))
abline(h=15/1596, col="red")
plot(rownames(Model),rstudent(DepFitReduced3))
abline(h=3,col="red")


#Predictor Transformations
a=Model$DepRate
b=Model$KerFeed
par(mfrow = c(1,2))
plot(a~b, data=Model, col="dodgerblue", pch=20, cex=1.5)
PreModel=lm(a~b,data=Model)
abline(PreModel, col="darkorange", lwd=2)
plot(fitted(PreModel),col="dodgerblue", pch=20, cex=1.5, xlab="Fitted", ylab="Residuals", ylim=c(-5,18))
abline(h=0,lty=1,col="darkorange", lwd=2)
#attempts
a=Model$DepRate
b=Model$ECPdP
par(mfrow = c(1,2))
plot(a~b, data=Model, col="dodgerblue", pch=20, cex=1.5)
PreModel=lm((a)~((b)),data=Model)
abline(PreModel, col="darkorange", lwd=2)
plot(fitted(PreModel),col="dodgerblue", pch=20, cex=1.5, xlab="Fitted", ylab="Residuals",)
abline(h=0,lty=1,col="darkorange", lwd=2)


#Test For Multicollinearity
#Heuristics, VIF>5 is significant. Try omitting these and see the effect on the moddle
#Recall Stepwise Log Transformation of Response Variable with highest Rsquared
#omit Reactor IP due to high multicollinearity
vif(DepFitReduced2)
DepFitTrans2=lm(log1p(DepRate) ~ KerTemp + WaterLev + BrineLev + WaterpH + 
                  WCBootLev + ECPdP + ECPLev +  ReactordP + ReactorCauLev + 
                  ClaydP , data = Model)
summary(DepFitTrans2)
par (mfrow=c(2,2),mar=c(2,2,2,2))
plot(DepFitTrans2)
#no more omission
ModelStep5=step(DepFitTrans2)

#PCA and PLS for Multicollinearity
library(pls)
library(caret)
library(dplyr)
set.seed(12345)
trainIndex = createDataPartition(Model$DepRate, p=0.9, list=FALSE)
DepRateTrain = Model[trainIndex,]
DepRateTest = Model[-trainIndex,]

PCAfit=pcr(log1p(DepRate) ~ (KerFeed) + (KerTemp) + (WaterLev) + (BrineLev) + (WaterpH) + 
             (WCBootLev) + (ECPdP) + (ECPLev) + (ReactorIP) + (ReactordP) + (ReactorCauLev) + 
             (ClaydP) + (ClayOP), data=Model, scale=TRUE, validation="CV")
summary(PCAfit)
par (mfrow=c(1,2))
plot(PCAfit$residuals, col="grey",pch=20)
abline(h=0,col="darkorange", lwd=2)
qqnorm(PCAfit$residuals, col="grey")
qqline(PCAfit$residuals,col="darkorange")
PCAfit$coefficients


# Plot the root mean squared error
par(mfrow = c(1,1))
validationplot(PCAfit, val.type = "R2")

predplot(PCAfit)
PCApred=predict(DepFitReduced2,Model)
Prediction=exp(PCApred)
View(Prediction)
write.table(Prediction,"C:/Users/rpdlpascual/Desktop/prediction3.txt",sep="\t")
