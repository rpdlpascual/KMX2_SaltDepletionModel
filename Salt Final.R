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
library(olsrr)
library(Hmisc)
library(corrplot)
library(rcompanion)
library(DescTools)
library(effects)
library(psych)
library(relaimpo)

Population = read.csv("SaltDepRateData.csv")
PopulationWithoutShutdowns = read.csv("SaltDepRateData3.csv")
set.seed(123456)
trainIndex = createDataPartition(Population$DepRate, p=0.80, list=FALSE)
TrainSet = Population[trainIndex,]
TestSet = Population[-trainIndex,]


trainIndex2 = createDataPartition(PopulationWithoutShutdowns$DepRate, p=0.90, list=FALSE)
TrainSet2 = PopulationWithoutShutdowns[trainIndex2,]
TestSet2 = PopulationWithoutShutdowns[-trainIndex2,]


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#DATA TRANSFORMATIONS
#Simple Linear Regression Model = SLRM
SLRM=lm(DepRate~.,data=TrainSet)
summary(SLRM)
#1. Check for Normality of Residuals
qqnorm(SLRM$residuals, col="dodgerblue", pch = 20, xlab="Theoretical Quantiles", ylab="Standardized Residuals", main="Model 1 - SLRM Normal Q-Q Plot")
qqline(SLRM$residuals,col="darkorange", lwd=2)
hist(SLRM$residuals, main = "Histogram of SLRM Residuals", xlab="SLRM Residuals", breaks=75, col="dodgerblue", border="black")


#2. Check for constant variance
plot(fitted(SLRM), resid(SLRM), col = "dodgerblue", pch = 20,
     xlab = "Predicted Values", ylab = "Residuals", main = "Model 1 - SLRM")
abline(h = 0, col = "darkorange", lwd = 2)

#3. Transformation of Response Variable
#We used ln(DepRate+C) where C = 15 since there are zero values for deprate
SLRMTrans=lm(log1p(DepRate+16)~., data=TrainSet)
summary(SLRMTrans)

qqnorm(SLRMTrans$residuals, col="dodgerblue", pch = 20, xlab="Theoretical Quantiles", ylab="Standardized Residuals", main="Model 2 - SLRM Trans Normal Q-Q Plot")
qqline(SLRMTrans$residuals,col="darkorange", lwd=2)

plot(fitted(SLRMTrans), resid(SLRMTrans), col = "dodgerblue", pch = 20,
     xlab = "Fitted", ylab = "Residuals", main = "Model 2 - SLRM Trans")
abline(h = 0, col = "darkorange", lwd = 2)

#Plot 2 - No Shutdowns
#No Additional constant since all values of DepRate is greater than zero
SLRMTransNoSD = lm(log1p(DepRate)~., data=TrainSet2)
summary(SLRMTransNoSD)

qqnorm(SLRMTransNoSD$residuals, col="dodgerblue", pch = 20, xlab="Theoretical Quantiles", ylab="Standardized Residuals", main="Model 2 - SLRM Trans Normal Q-Q Plot")
qqline(SLRMTransNoSD$residuals,col="darkorange", lwd=2)

plot(fitted(SLRMTransNoSD), resid(SLRMTransNoSD), col = "dodgerblue", pch = 20,
     xlab = "Fitted", ylab = "Residuals", main = "Model 3 - SLRM Trans No SD")
abline(h = 0, col = "darkorange", lwd = 2)



#4. Check for Influential Points
LeveragePoint=abs(hatvalues(SLRMTransNoSD))
plot(LeveragePoint, col="dodgerblue", ylab="Leverage Statistic", xlab="DepRate Observations", main = "Model 3 - Leverage Points")
abline(h=(2*(21+1)/1440), col="darkorange", lwd=2)
which(LeveragePoint>(2*(21+1)/1440),)

Outliers=abs(rstudent(SLRMTransNoSD))
plot(Outliers, col="dodgerblue", ylab="Studentized Residuals", xlab="DepRate Observations", main = "Model 3 - Outliers")
abline(h=3, col="darkorange", lwd=2)
which(Outliers>3,)


ols_plot_resid_lev(SLRMTransNoSD)
InfluentialPlot=ols_plot_resid_lev(SLRMTransNoSD)
Influential=InfluentialPlot$data
write.table(Influential,"C:/Users/rpdlpascual/Desktop/Influential.txt",sep="\t" )
write.table(TrainSet2,"C:/Users/rpdlpascual/Desktop/TrainSet2.txt",sep="\t" )

TrainSet2NoInfluential = read.csv("SaltDepRateData5.csv")           
SLRMTransNoSDNoInfluential = lm(log1p(DepRate)~., data=TrainSet2NoInfluential)
summary(SLRMTransNoSDNoInfluential)

qqnorm(SLRMTransNoSDNoInfluential$residuals, col="dodgerblue", pch = 20, xlab="Theoretical Quantiles", ylab="Standardized Residuals")
qqline(SLRMTransNoSDNoInfluential$residuals,col="darkorange", lwd=2)

plot(fitted(SLRMTransNoSDNoInfluential), resid(SLRMTransNoSDNoInfluential), col = "dodgerblue", pch = 20,
     xlab = "Fitted", ylab = "Residuals", main = "Model 4 - SLRM No Influential", xlim=c(1.65,3), ylim=c(-0.2,0.2))
abline(h = 0, col = "darkorange", lwd = 2)

#5. Check for Multicollinearity
flattenCorrMatrix = function(cormat,pmat){
  ut = upper.tri(cormat)
  data.frame(
    row=rownames(cormat)[row(cormat)[ut]],
    column=rownames(cormat)[col(cormat)[ut]],
    cor=(cormat)[ut],
    p = pmat[ut]
  )
}
res=rcorr(as.matrix(TrainSet2[,10:22]))
Correlation=flattenCorrMatrix(res$r,res$P)
write.table(Correlation,"C:/Users/rpdlpascual/Desktop/Correlation.txt",sep="\t")
#Correlation Graph, insignificant correlations are crossed
corrplot(res$r, type="upper", order="hclust", p.mat=res$P, sig.level=0.05, insig="blank")

#Include interactions with significant p-value (which conceptually makes sense)
SLRMTransNoSDWithInt=lm(log1p(DepRate)~ D1+D2+D3+D4+D5+D6+D7+D8+
                               KerFeed+ KerTemp + WaterLev+BrineLev + WaterpH+WCBootLev+
                               ECPdP + ECPLev+ReactorIP+ReactordP + ReactorCauLev+
                               ClaydP+ClayOP+ReactorIP:ClayOP+KerFeed:ReactorIP+
                               ReactorIP:ClaydP+KerFeed:ClaydP+KerFeed:ClayOP+
                               ReactorIP:ReactordP+ReactordP:ClayOP+KerFeed:KerTemp+
                               KerFeed:ReactordP+WaterLev:WaterpH+KerTemp:ReactorIP+
                               ReactordP:ClaydP+KerFeed:BrineLev+KerTemp:ClayOP+
                               ClaydP:ClayOP, data=TrainSet2)
summary(SLRMTransNoSDWithInt)

SLRMTransNoSDNoInfWithInt=lm(log1p(DepRate)~ D1+D2+D3+D4+D5+D6+D7+D8+
                          KerFeed+KerTemp+WaterLev+BrineLev + WaterpH+WCBootLev+
                          ECPdP + ECPLev+ReactorIP+ReactordP + ReactorCauLev+
                          ClaydP+ClayOP+ReactorIP:ClayOP+KerFeed:ReactorIP+
                          ReactorIP:ClaydP+KerFeed:ClaydP+KerFeed:ClayOP+
                          ReactorIP:ReactordP+ReactordP:ClayOP+KerFeed:KerTemp+
                          KerFeed:ReactordP+WaterLev:WaterpH+KerTemp:ReactorIP+
                          ReactordP:ClaydP+KerFeed:BrineLev+KerTemp:ClayOP+
                          ClaydP:ClayOP, data=TrainSet2NoInfluential)
summary(SLRMTransNoSDNoInfWithInt)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Model Accuracy and Stepwise Linear Regression
#Model 3
summary(SLRMTransNoSD)
PredictModel3 = exp(predict(SLRMTransNoSD,TestSet2))
MAE(PredictModel3, TestSet2$DepRate)
MSE(PredictModel3, TestSet2$DepRate)
RMSE(PredictModel3, TestSet2$DepRate)

#Model4
summary(SLRMTransNoSDNoInfluential)
PredictModel4 = exp(predict(SLRMTransNoSDNoInfluential,TestSet2))
MAE(PredictModel4, TestSet2$DepRate)
MSE(PredictModel4, TestSet2$DepRate)
RMSE(PredictModel4, TestSet2$DepRate)

#Model5
summary(SLRMTransNoSDWithInt)
PredictModel5 = exp(predict(SLRMTransNoSDWithInt,TestSet2))
MAE(PredictModel5, TestSet2$DepRate)
MSE(PredictModel5, TestSet2$DepRate)
RMSE(PredictModel5, TestSet2$DepRate)

#Model5a
summary(SLRMTransNoSDNoInfWithInt)
PredictModel5a = exp(predict(SLRMTransNoSDNoInfWithInt,TestSet2))
MAE(PredictModel5a, TestSet2$DepRate)
MSE(PredictModel5a, TestSet2$DepRate)
RMSE(PredictModel5a, TestSet2$DepRate)


#Model6 - Stepwise Model 5a
ModelStep = step(SLRMTransNoSDNoInfWithInt)
SLRMTransNoSDNoInfWithIntStepReduced = lm(log1p(DepRate) ~ D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + KerFeed + 
                                            KerTemp + WaterLev + BrineLev + WaterpH + WCBootLev + ECPdP + 
                                            ReactorIP + ReactordP + ClaydP + ClayOP + KerFeed:ReactorIP + 
                                            ReactorIP:ClaydP + KerFeed:ClaydP + KerFeed:ClayOP + ReactorIP:ReactordP + 
                                            ReactordP:ClayOP + KerFeed:KerTemp + KerFeed:ReactordP + 
                                            KerTemp:ReactorIP + KerFeed:BrineLev + KerTemp:ClayOP + ClaydP:ClayOP, data = TrainSet2NoInfluential)
summary(SLRMTransNoSDNoInfWithIntStepReduced)
Model6=SLRMTransNoSDNoInfWithIntStepReduced$coefficients
write.table(Model6,"C:/Users/rpdlpascual/Desktop/Model6.txt",sep="\t")
PredictModel6 = exp(predict(SLRMTransNoSDNoInfWithIntStepReduced,TestSet2))
MAE(PredictModel6, TestSet2$DepRate)
MSE(PredictModel6, TestSet2$DepRate)
RMSE(PredictModel6, TestSet2$DepRate)


#Check VIF of Model 6 - Model 7
VIF(SLRMTransNoSDNoInfWithIntStepReduced)
#Remove ReactorIP with VIF of 19.311
SLRMTransNoSDNoInfWithIntStepReduced2 = lm(log1p(DepRate)~D1+D2+D3+D4+D5+D6+D7+D8+KerFeed+KerTemp+
                                           WaterLev+BrineLev+WaterpH+WCBootLev+ECPdP+ReactordP+ClaydP+ClayOP+
                                             KerFeed:ClaydP + KerFeed:ClayOP + ReactordP:ClayOP + KerFeed:KerTemp + KerFeed:ReactordP +  
                                             KerFeed:BrineLev + ClayOP:KerTemp + ClaydP:ClayOP, data = TrainSet2NoInfluential )
VIF(SLRMTransNoSDNoInfWithIntStepReduced2)



#Remove ClaydP with VIF of 13.56
SLRMTransNoSDNoInfWithIntStepReduced3 = lm(log1p(DepRate) ~ D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + KerFeed + KerTemp +
                                             WaterLev + BrineLev + WaterpH + WCBootLev + ECPdP+ ReactordP + ClayOP + 
                                             KerFeed:ClayOP + ReactordP:ClayOP + KerFeed:KerTemp + KerFeed:ReactordP +  
                                             KerFeed:BrineLev + ClayOP:KerTemp, data = TrainSet2NoInfluential)
summary(SLRMTransNoSDNoInfWithIntStepReduced3)
VIF(SLRMTransNoSDNoInfWithIntStepReduced3)


ModelStep2=step(SLRMTransNoSDNoInfWithIntStepReduced3)

#Model 7
SLRMTransFinalModel = lm(log1p(DepRate)~D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + KerFeed + 
                           KerTemp + BrineLev + WCBootLev + ECPdP + ReactordP + ClayOP + 
                           KerFeed:ClayOP + ReactordP:ClayOP + KerFeed:KerTemp + KerFeed:ReactordP + 
                           KerFeed:BrineLev + KerTemp:ClayOP, data = TrainSet2NoInfluential)
summary(SLRMTransFinalModel)
VIF(SLRMTransFinalModel)


Model7=SLRMTransFinalModel$coefficients
write.table(Model7,"C:/Users/rpdlpascual/Desktop/Model7Coeffs.txt",sep="\t")


PredictModel7 = exp(predict(SLRMTransFinalModel,TestSet2))
MAE(PredictModel7, TestSet2$DepRate)
MSE(PredictModel7, TestSet2$DepRate)
RMSE(PredictModel7, TestSet2$DepRate)


#Model 8 - Simplified Model 7
SLRMTransFinalModelSimplified = lm(log1p(DepRate)~D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + KerFeed + 
                           ReactordP + ClayOP + KerFeed:ClayOP + ReactordP:ClayOP + KerFeed:ReactordP
                           , data = TrainSet2NoInfluential)
summary(SLRMTransFinalModelSimplified)
Model8=SLRMTransFinalModelSimplified$coefficients
write.table(Model8,"C:/Users/rpdlpascual/Desktop/Model8Coeffs.txt",sep="\t")

VIF(SLRMTransFinalModelSimplified)
PredictModel8 = exp(predict(SLRMTransFinalModelSimplified,TestSet2))
MAE(PredictModel8, TestSet2$DepRate)
MSE(PredictModel8, TestSet2$DepRate)
RMSE(PredictModel8, TestSet2$DepRate)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Model 8A
DummyVarInt=read.csv("SaltDepRateData7.csv")
DumVarMod = lm(log1p(DepRate)~D1+D2+D3+D4+D5+D6+D7+D8+Theta+
                              D1:Theta+D2:Theta+D3:Theta+
                              D4:Theta+D5:Theta+D6:Theta+
                              D7:Theta+D8:Theta, data=DummyVarInt)
summary(DumVarMod)
Model8a=DumVarMod$coefficients
write.table(Model8a,"C:/Users/rpdlpascual/Desktop/Model8aCoeffs.txt",sep="\t")
attach(DummyVarInt)

#Edit SaltDepRateData7 and Make zeroes in D8 to 2 to be able to detect D9
plot(Theta[D1=="1"], log1p(DepRate)[D1=="1"], col=rgb(1,0.4,0.1,0.5), xlab="Phi", ylab="ln(DepRate)", 
     main="DepRate vs P,TimeFrame", pch=16, ylim=c(1.6,3.1), xlim=c(-0.4,0.7))
points(Theta[D2=="1"], log1p(DepRate)[D2=="1"],col="dodgerblue", pch=16)
points(Theta[D3=="1"], log1p(DepRate)[D3=="1"],col=rgb(0.8,0,0.9), pch=16)
points(Theta[D4=="1"], log1p(DepRate)[D4=="1"],col=rgb(0,0,0,0.75), pch=16)
points(Theta[D5=="1"], log1p(DepRate)[D5=="1"],col=rgb(0.5,1,0.1,0.55),pch=16)
points(Theta[D6=="1"], log1p(DepRate)[D6=="1"],col=rgb(1,0,0), pch=16)
points(Theta[D7=="1"], log1p(DepRate)[D7=="1"],col=rgb(0.7,0.3,0.1,0.65), pch=16)
points(Theta[D8=="1"], log1p(DepRate)[D8=="1"],col=rgb(0.5,0.5,0.6), pch=16)
points(Theta[D8=="2"], log1p(DepRate)[D8=="2"],col=rgb(1,0.8,0,0.65), pch=16)
legend(-0.42,3.0,legend=c("D1","D2","D3","D4","D5","D6","D7","D8","D9"), col=c(rgb(1,0.4,0.1,0.5),"dodgerblue",rgb(0.8,0,0.9),
                          rgb(0,0,0,0.75),rgb(0.5,1,0.1,0.55),
                          rgb(1,0,0),rgb(0.7,0.3,0.1,0.65),
                          rgb(0.5,0.5,0.6),rgb(1,0.8,0,0.65)),pch=16, cex=1.12)


#substituted equation lines
plot(Theta[D1=="1"], log1p(DepRate)[D1=="1"], col=rgb(0,0,0,0), xlab="Phi", ylab="ln(DepRate)", 
     main="DepRate vs P,TimeFrame", pch=16, ylim=c(1.6,3.1), xlim=c(-0.8,0.7))

abline(a=2.221, b=1.105,col=rgb(1,0.4,0.1), lwd=2)
abline(a=1.855, b=0.986,col="dodgerblue", lwd=2)
abline(a=2.475, b=0.972,col=rgb(0.8,0,0.9), lwd=2)
abline(a=2.262, b=0.973,col=rgb(0,0,0), lwd=2)
abline(a=2.240, b=1.010,col=rgb(0.5,1,0.1), lwd=2)
abline(a=2.157, b=0.904,col=rgb(1,0,0), lwd=2)
abline(a=2.223, b=0.945,col=rgb(0.7,0.3,0.1), lwd=2)
abline(a=2.764, b=1.091,col=rgb(0.5,0.5,0.6), lwd=2)
abline(a=2.439, b=0.991,col=rgb(1,0.8,0), lwd=2)
legend(-0.8,3.0,legend=c("D1","D2","D3","D4","D5","D6","D7","D8","D9"), col=c(rgb(1,0.4,0.1),"dodgerblue",rgb(0.8,0,0.9),
                            rgb(0,0,0),rgb(0.5,1,0.1),
                            rgb(1,0,0),rgb(0.7,0.3,0.1),
                            rgb(0.5,0.5,0.6),rgb(1,0.8,0)), lwd = 3, cex=0.8)

#DepRate vs Kerfeed,TimeFrame
DummyVarInt2=read.csv("SaltDepRateData9.csv")
attach(DummyVarInt2)
plot(KerFeed[D1=="1"], log1p(DepRate)[D1=="1"], col=rgb(1,0.4,0.1,0.5), xlab="KerFeed", ylab="ln(DepRate)", 
     main="DepRate vs KerFeed,TimeFrame", pch=16, ylim=c(1.5,3.1), xlim=c(-1,1.7))
points(KerFeed[D2=="1"], log1p(DepRate)[D2=="1"],col="dodgerblue", pch=16)
points(KerFeed[D3=="1"], log1p(DepRate)[D3=="1"],col=rgb(0.8,0,0.9), pch=16)
points(KerFeed[D4=="1"], log1p(DepRate)[D4=="1"],col=rgb(0,0,0,0.75), pch=16)
points(KerFeed[D5=="1"], log1p(DepRate)[D5=="1"],col=rgb(0.5,1,0.1,0.55),pch=16)
points(KerFeed[D6=="1"], log1p(DepRate)[D6=="1"],col=rgb(1,0,0), pch=16)
points(KerFeed[D7=="1"], log1p(DepRate)[D7=="1"],col=rgb(0.7,0.3,0.1,0.65), pch=16)
points(KerFeed[D8=="1"], log1p(DepRate)[D8=="1"],col=rgb(0.5,0.5,0.6), pch=16)
points(KerFeed[D8=="2"], log1p(DepRate)[D8=="2"],col=rgb(1,0.8,0,0.65), pch=16)
legend(-1,3.0,legend=c("D1","D2","D3","D4","D5","D6","D7","D8","D9"), col=c(rgb(1,0.4,0.1,0.5),"dodgerblue",rgb(0.8,0,0.9),
                                                                               rgb(0,0,0,0.75),rgb(0.5,1,0.1,0.55),
                                                                               rgb(1,0,0),rgb(0.7,0.3,0.1,0.65),
                                                                               rgb(0.5,0.5,0.6),rgb(1,0.8,0,0.65)),pch=16, cex=1.12)
#DepRate vs ReactordP,TimeFrame
plot(ReactordP[D1=="1"], log1p(DepRate)[D1=="1"], col=rgb(1,0.4,0.1,0.5), xlab="ReactordP", ylab="ln(DepRate)", 
     main="DepRate vs ReactordP,TimeFrame", pch=16, xlim=c(-0.7,0.1), ylim=c(1.5,3.1))
points(ReactordP[D2=="1"], log1p(DepRate)[D2=="1"],col="dodgerblue", pch=16)
points(ReactordP[D3=="1"], log1p(DepRate)[D3=="1"],col=rgb(0.8,0,0.9), pch=16)
points(ReactordP[D4=="1"], log1p(DepRate)[D4=="1"],col=rgb(0,0,0,0.75), pch=16)
points(ReactordP[D5=="1"], log1p(DepRate)[D5=="1"],col=rgb(0.5,1,0.1,0.55),pch=16)
points(ReactordP[D6=="1"], log1p(DepRate)[D6=="1"],col=rgb(1,0,0), pch=16)
points(ReactordP[D7=="1"], log1p(DepRate)[D7=="1"],col=rgb(0.7,0.3,0.1,0.65), pch=16)
points(ReactordP[D8=="1"], log1p(DepRate)[D8=="1"],col=rgb(0.5,0.5,0.6), pch=16)
points(ReactordP[D8=="2"], log1p(DepRate)[D8=="2"],col=rgb(1,0.8,0,0.65), pch=16)
legend(-0.7,3.0,legend=c("D1","D2","D3","D4","D5","D6","D7","D8","D9"), col=c(rgb(1,0.4,0.1,0.5),"dodgerblue",rgb(0.8,0,0.9),
                                                                               rgb(0,0,0,0.75),rgb(0.5,1,0.1,0.55),
                                                                               rgb(1,0,0),rgb(0.7,0.3,0.1,0.65),
                                                                               rgb(0.5,0.5,0.6),rgb(1,0.8,0,0.65)),pch=16, cex=1.12)
#DepRate vs ClayOP,TimeFrame
plot(ClayOP[D1=="1"], log1p(DepRate)[D1=="1"], col=rgb(1,0.4,0.1,0.5), xlab="ClayOP", ylab="ln(DepRate)", 
     main="DepRate vs ClayOP,TimeFrame", pch=16, ylim=c(1.5,3.2), xlim=c(-0.2,1))
points(ClayOP[D2=="1"], log1p(DepRate)[D2=="1"],col="dodgerblue", pch=16)
points(ClayOP[D3=="1"], log1p(DepRate)[D3=="1"],col=rgb(0.8,0,0.9), pch=16)
points(ClayOP[D4=="1"], log1p(DepRate)[D4=="1"],col=rgb(0,0,0,0.75), pch=16)
points(ClayOP[D5=="1"], log1p(DepRate)[D5=="1"],col=rgb(0.5,1,0.1,0.55),pch=16)
points(ClayOP[D6=="1"], log1p(DepRate)[D6=="1"],col=rgb(1,0,0), pch=16)
points(ClayOP[D7=="1"], log1p(DepRate)[D7=="1"],col=rgb(0.7,0.3,0.1,0.65), pch=16)
points(ClayOP[D8=="1"], log1p(DepRate)[D8=="1"],col=rgb(0.5,0.5,0.6), pch=16)
points(ClayOP[D8=="2"], log1p(DepRate)[D8=="2"],col=rgb(1,0.8,0,0.65), pch=16)
legend(-0.2,3.2,legend=c("D1","D2","D3","D4","D5","D6","D7","D8","D9"), col=c(rgb(1,0.4,0.1,0.5),"dodgerblue",rgb(0.8,0,0.9),
                                                                              rgb(0,0,0,0.75),rgb(0.5,1,0.1,0.55),
                                                                              rgb(1,0,0),rgb(0.7,0.3,0.1,0.65),
                                                                              rgb(0.5,0.5,0.6),rgb(1,0.8,0,0.65)),pch=16, cex=1.12)

attach(TrainSet2NoInfluential)
plot(KerFeed, ReactordP, col=rgb(1,0.4,0.1,0.5), xlab="KerFeed", ylab="ReactordP", 
     main="KerFeed vs ReactordP", pch=16, xlim=c(-0.7,1.7), ylim=c(-0.7,0.2))


plot(KerFeed, ClayOP, col="dodgerblue", xlab="KerFeed", ylab="ClayOP", 
     main="KerFeed vs ClayOP", pch=16, xlim=c(-2,5), ylim=c(-5,5))

plot(ClayOP, ReactordP, col=rgb(0.8,0,0.9,0.5), xlab="ClayOP", ylab="ReactordP", 
     main="ClayOP vs ReactordP", pch=16, xlim=c(-0.1,0.9), ylim=c(-0.6,0.2))


plot((0.407*KerFeed+0.073*ClayOP)[D3=="1"], log1p(DepRate)[D3=="1"], col=rgb(1,0.4,0.1,0.5), xlab="Z(KerFeed,ClayOP)", ylab="ln(DepRate)", 
     main="Z vs ln(DepRate)", pch=16, xlim=c(-0.2,0.7))
points((0.407*KerFeed+0.073*ClayOP-(0.097*KerFeed*ClayOP))[D3=="1"], log1p(DepRate)[D3=="1"],col="dodgerblue", pch=16)
legend(-0.2,3,legend=c("w/o Interactions","w/ Interactions"), 
       col=c(rgb(1,0.4,0.1,0.5),"dodgerblue"),pch=16, cex=1.12)

plot((-0.021*ReactordP+0.073*ClayOP)[D3=="1"], log1p(DepRate)[D3=="1"], col=rgb(1,0.4,0.1,0.5), xlab="Z(ReactordP,ClayOP)", ylab="ln(DepRate)", 
     main="Z vs ln(DepRate)", pch=16, xlim=c(0.039,0.048))
points((-0.021*ReactordP+0.073*ClayOP+(0.035*ReactordP*ClayOP))[D3=="1"], log1p(DepRate)[D3=="1"],col="dodgerblue", pch=16)
legend(0.0447,3,legend=c("w/o Interactions","w/ Interactions"), 
       col=c(rgb(1,0.4,0.1,0.5),"dodgerblue"),pch=16, cex=1.12)

plot((0.407*KerFeed-0.021*ReactordP)[D3=="1"], log1p(DepRate)[D3=="1"], col=rgb(1,0.4,0.1,1), xlab="Z(KerFeed,ReactordP)", ylab="ln(DepRate)", 
     main="Z vs ln(DepRate)", pch=16, xlim=c(-0.2,0.7))
points((0.407*KerFeed-0.021*ReactordP-(0.03*KerFeed*ReactordP))[D3=="1"], log1p(DepRate)[D3=="1"],col=rgb(0.12,0.565,0.9,0.65), pch=16)
legend(-0.2,3,legend=c("w/o Interactions","w/ Interactions"), 
       col=c(rgb(1,0.4,0.1,1),rgb(0.12,0.565,0.9,1)),pch=16, cex=1.12)

plot((0.407*KerFeed-0.021*ReactordP+0.073*ClayOP), log1p(DepRate), 
     col=rgb(1,0.4,0.1), xlab="Total Z", ylab="ln(DepRate)", ylim=c(1.25,3.1), 
     main="Total Z vs ln(DepRate)", pch=16, xlim=c(-0.2,0.7))
points((0.407*KerFeed-0.021*ReactordP+0.073*ClayOP+(-0.03*KerFeed*ReactordP)+
     (0.035*ReactordP*ClayOP)-(0.097*KerFeed*ClayOP)),
      log1p(DepRate),col=rgb(0.12,0.565,0.9,0.4), pch=16)
legend(0.35,1.85,legend=c("w/o Interactions","w/ Interactions"), 
       col=c(rgb(1,0.4,0.1,0.5),rgb(0.12,0.565,0.9,0.4)),pch=16, cex=1.12)


#Relative Importance

#Model 9
NoDummy = read.csv("SaltDepRateData11.csv")
SLRMTransNoDummy = lm(log1p(DepRate)~., data = NoDummy)
summary(SLRMTransNoDummy)
step(SLRMTransNoDummy)
SLRMTransNoDummyReduced = lm(log1p(DepRate) ~ KerFeed + KerTemp + WaterLev + BrineLev+
                          WaterpH + ECPLev + ReactorIP + ReactordP + ClaydP + 
                          ClayOP + KerFeed:ReactorIP + KerFeed:KerTemp+
                          KerFeed:ReactordP + KerFeed:ClaydP + KerFeed:ClayOP+
                          KerFeed:BrineLev+WaterpH:WaterLev+ReactorIP:ReactordP+
                          ReactorIP:ClaydP + ReactorIP:ClayOP + ReactordP:ClaydP+
                          ReactordP:ClayOP + ClaydP:ClayOP, data = NoDummy)
summary(SLRMTransNoDummyReduced)
step(SLRMTransNoDummyReduced)
SLRMTransNoDummyReduced2 = lm(log1p(DepRate) ~ KerFeed + KerTemp + WaterLev + BrineLev + 
                           WaterpH + ReactorIP + ReactordP + ClaydP + ClayOP + 
                           KerFeed:ReactorIP + KerFeed:ClaydP + KerFeed:ClayOP + 
                           KerFeed:BrineLev + ReactorIP:ClaydP + ReactordP:ClaydP + 
                           ReactordP:ClayOP + ClaydP:ClayOP, data = NoDummy)
summary(SLRMTransNoDummyReduced2)
vif(SLRMTransNoDummyReduced2)
#Remove ReactorIP
SLRMTransNoDummyReduced3 = lm(log1p(DepRate) ~ KerFeed + KerTemp + WaterLev + BrineLev + 
                                WaterpH + ReactordP + ClaydP + ClayOP + 
                                KerFeed:ClaydP + KerFeed:ClayOP + 
                                KerFeed:BrineLev + ReactordP:ClaydP + 
                                ReactordP:ClayOP + ClaydP:ClayOP, data = NoDummy)
summary(SLRMTransNoDummyReduced3)
vif(SLRMTransNoDummyReduced3)
step(SLRMTransNoDummyReduced3)

SLRMTransNoDummyReduced4 = lm(log1p(DepRate) ~ KerFeed + KerTemp + BrineLev + 
                                WaterpH + ReactordP + ClaydP + ClayOP + KerFeed:ClaydP + 
                                KerFeed:ClayOP + KerFeed:BrineLev+ ReactordP:ClayOP + 
                                ClaydP:ClayOP, data = NoDummy)
summary(SLRMTransNoDummyReduced4)
PredictModel9 = exp(predict(SLRMTransNoDummyReduced4,TestSet2))
MAE(PredictModel9, TestSet2$DepRate)
MSE(PredictModel9, TestSet2$DepRate)
RMSE(PredictModel9, TestSet2$DepRate)

#Relative Importance
calc.relimp(SLRMTransFinalModelSimplified, rela=TRUE)
calc.relimp(SLRMTransFinalModelSimplified, rela=FALSE)

#Accuracy Comparison
#Model0
Base = read.csv("SaltDepRateData13.csv")
SLRMNoStandardizationNoDummy = lm(DepRate~., data = Base)
summary(SLRMNoStandardizationNoDummy)
Pred0 = predict(SLRMNoStandardizationNoDummy,Base)
write.table(Pred0,"C:/Users/rpdlpascual/Desktop/Pred0.txt",sep="\t")

#Model3
Pred3 = exp(predict(SLRMTransNoSD,Population))
write.table(Pred3,"C:/Users/rpdlpascual/Desktop/Pred3.txt",sep="\t")

#Model4
Pred4 = exp(predict(SLRMTransNoSDNoInfluential,Population))
write.table(Pred4,"C:/Users/rpdlpascual/Desktop/Pred4.txt",sep="\t")

#Model5
Pred5 = exp(predict(SLRMTransNoSDWithInt,Population))
write.table(Pred5,"C:/Users/rpdlpascual/Desktop/Pred5.txt",sep="\t")

#Model5a
Pred5a = exp(predict(SLRMTransNoSDNoInfWithInt,Population))
write.table(Pred5a,"C:/Users/rpdlpascual/Desktop/Pred5a.txt",sep="\t")

#Model6
Pred6 = exp(predict(SLRMTransNoSDNoInfWithIntStepReduced,Population))
write.table(Pred6,"C:/Users/rpdlpascual/Desktop/Pred6.txt",sep="\t")

#Model7
Pred7=exp(predict(SLRMTransFinalModel,Population))
write.table(Pred7,"C:/Users/rpdlpascual/Desktop/Pred7.txt",sep="\t")

#Model8
Pred8=exp(predict(SLRMTransFinalModelSimplified,Population))
write.table(Pred8,"C:/Users/rpdlpascual/Desktop/Pred8.txt",sep="\t")

#Model9
Pred9=exp(predict(SLRMTransNoDummyReduced4,Population))
write.table(Pred9,"C:/Users/rpdlpascual/Desktop/Pred9.txt",sep="\t")