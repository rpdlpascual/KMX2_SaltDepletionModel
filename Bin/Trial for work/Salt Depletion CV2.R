#Ordinary Least Squares Regression
#Ordinary Least Squares (OLS) regression is a linear model that seeks to find a set of coefficients for a line/hyper-plane that minimise the sum of the squared errors.
Model=read.csv("DepModelCV2.csv")
DepFit=lm(DepRate~.,data=Model)
summary(DepFit)
plot(DepFit)
attach(Model)
par(mfrow=c(1,1)
plot(WaterLev,DepRate, xlab="KerFeed", ylab="DepRate")

#Log Transformation
DepFitTrans=lm(log1p(DepRate) ~ (KerFeed) + (KerTemp) + (WaterLev) + (BrineLev) + (WaterpH), data=Model)
summary(DepFitTrans)
plot(DepFitTrans)

plot(fitted(DepFitTrans), resid(DepFitTrans), col="grey", pch=20, xlab="Fitted", ylab="Residuals", main="Fitted versus Residuals")
shapiro.test(DepFitTrans$residuals)
library(zoo)
library(lmtest)
bptest(DepFitTrans)


#BoxCox Transformation
library(MASS)
DepFitAbs=lm(abs(DepRate)~., data=Model)
boxcox(DepFitAbs, lambda = seq(-0.5, 1, by=0.05), plotit=TRUE)

#Test For Multicollinearity
library(car)
vif(DepFit)

#Stepwise Regression - finding the best combination of variables 
#Stepwise Linear Regression is a method that makes use of linear regression to discover which subset of attributes in the dataset result in the best performing model. It is step-wise because each iteration of the method makes a change to the set of attributes and creates a model to evaluate the performance of the set.
ModelStep=step(DepFit)
ModelStep2=step(DepFit, scope=list(lower=~1,upper=~KerFeed+KerTemp+WaterLev+BrineLev+WaterpH, direction="both", trace=1))
DepFitReduced = lm(DepRate~KerFeed + KerTemp + WaterLev + WaterpH, data=Model) 
summary(DepFitReduced)
shapiro.test(DepFitReduced$residuals)
bptest(DepFitReduced)

#Model Adequacy Checking - mathematical assumption verification
#Check Normality
par (mfrow=c(2,2),mar=c(2,2,2,2))
plot (DepFitReduced)

plot(hatvalues(DepFitReduced))
abline(h=7/1254, col="red")

plot(rownames(Model),rstudent(DepFitReduced))
abline(h=3,col="red")
rstudent(DepFitReduced)

#K-folds Cross Validating
library(lattice)
library(ggplot2)
library(caret)
set.seed(123)
trctrl = trainControl(method="repeatedcv", number=10,repeats="10")
ModelCV = train(DepRate~.,data=Model,method="lm",trControl=trctrl)
print(ModelCV)
ModelCV$finalModel

ModelCVReduced = train(DepRate~KerFeed + KerTemp + WaterLev + WaterpH, data=Model, method="lm", trControl=trctrl)
summary(ModelCVReduced)
ModelCVReduced$finalModel


