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
library(PerformanceAnalytics)

RealModel = read.csv ("SaltDepRateData.csv")
Model = read.csv("SaltDepRateData3.csv")
DepFit=lm((DepRate)~.,data=Model)
summary(DepFit)
par (mfrow=c(1,1))
plot(fitted(DepFit), resid(DepFit), col = "dodgerblue", pch = 20,
     xlab = "Fitted", ylab = "Residuals", main = "Simple Linear Model")
abline(h = 0, col = "darkorange", lwd = 2)

qqnorm(DepFit$residuals, col="dodgerblue", pch = 20, xlab="Theoretical Quantiles", ylab="Standardized Residuals")
qqline(DepFit$residuals,col="darkorange", lwd=2)

#Transformation
DepFitTrans=lm(log1p(DepRate)~., data=Model)
summary(DepFitTrans)
par (mfrow=c(1,1))
plot(fitted(DepFitTrans), resid(DepFitTrans), col = "dodgerblue", pch = 20,
     xlab = "Fitted", ylab = "Residuals", main = "Transformed SLRM")
abline(h = 0, col = "darkorange", lwd = 2)

qqnorm(DepFitTrans$residuals, col="dodgerblue", pch = 20, xlab="Theoretical Quantiles", ylab="Standardized Residuals")
qqline(DepFitTrans$residuals,col="darkorange", lwd=2)

#Leverages and Outliers
par(mfrow = c(1,2))
#Plot for detecting outliers. Studentized deleted residuals 
#(or externally studentized residuals) is the deleted residual divided
#by its estimated standard deviation. Studentized residuals are going to 
#be more effective for detecting outlying Y observations than standardized 
#residuals. If an observation has an externally studentized residual that is 
#larger than 3 (in absolute value) we can call it an outlier.
ols_plot_resid_stud(DepFitTrans)
#outliers - R-student (y-values)
x=abs(rstudent(DepFitTrans))
which(x>3,)

#leverage (x-values), criteria (2*(k+1)/n) where k = # of predictors and n = #values
z=abs(hatvalues(DepFitTrans))
which(z>(44/1596),)
#zvsx plot
ols_plot_resid_lev(DepFitTrans)

#----------------------------------------------------------------
#----------------------------------------------------------------
#New Data without Influential Points
Model2 = read.csv("SaltDepRateData4.csv")
DepFit2=lm((DepRate)~.,data=Model2)
summary(DepFit2)
par (mfrow=c(1,1))
plot(fitted(DepFit2), resid(DepFit2), col = "dodgerblue", pch = 20,
     xlab = "Fitted", ylab = "Residuals", main = "Simple Linear Model")
abline(h = 0, col = "darkorange", lwd = 2)

qqnorm(DepFit2$residuals, col="dodgerblue", pch = 20, xlab="Theoretical Quantiles", ylab="Standardized Residuals")
qqline(DepFit2$residuals,col="darkorange", lwd=2)

#Transformation
DepFitTrans2=lm(log1p(DepRate)~., data=Model2)
summary(DepFitTrans2)
par (mfrow=c(1,1))
plot(fitted(DepFitTrans2), resid(DepFitTrans2), col = "dodgerblue", pch = 20,
     xlab = "Fitted", ylab = "Residuals", main = "Transformed SLRM")
abline(h = 0, col = "darkorange", lwd = 2)

qqnorm(DepFitTrans2$residuals, col="dodgerblue", pch = 20, xlab="Theoretical Quantiles", ylab="Standardized Residuals")
qqline(DepFitTrans2$residuals,col="darkorange", lwd=2)
ols_plot_resid_lev(DepFitTrans2)

ModelStep3=step(DepFitTrans2)
DepFitTrans2Reduced = lm(log1p(DepRate) ~ D1 + D2 + D4 + D5 + D6 + D7 + D8 + KerFeed + 
                           KerTemp + WaterLev + BrineLev + WCBootLev + ECPdP + ECPLev + 
                           ReactorIP + ClaydP + ClayOP, data=Model)
summary(DepFitTrans2Reduced)
Predict3=exp(predict(DepFitTrans2Reduced,RealModel))
write.table(Predict3,"C:/Users/rpdlpascual/Desktop/Predict3.txt",sep="\t" )

#----------------------------------------------------------------
#----------------------------------------------------------------
#Test For Multicollinearity
#Heuristics, VIF>5 is significant. Try omitting these and see the effect on the moddle
#Recall Stepwise Log Transformation of Response Variable with highest Rsquared
#omit Reactor IP due to high multicollinearity
w=vif(DepFitTrans)
which(w>5,)
#list of all interactions
flattenCorrMatrix = function(cormat,pmat){
  ut = upper.tri(cormat)
  data.frame(
    row=rownames(cormat)[row(cormat)[ut]],
    column=rownames(cormat)[col(cormat)[ut]],
    cor=(cormat)[ut],
    p = pmat[ut]
  )
}
attach(Model)
res=rcorr(as.matrix(Model[,2:22]))
k=flattenCorrMatrix(res$r,res$P)
write.table(k,"C:/Users/rpdlpascual/Desktop/correlation.txt",sep="\t")
#Correlation Graph, insignificant correlations are crossed
corrplot(res$r, type="upper", order="hclust", p.mat=res$P, sig.level=0.05, insig="blank")
help(corrplot)

#NewModel Accounting for Multicollinearity
DepFitTransInt=lm(log1p(DepRate)~ D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8
                  + (KerFeed) + (KerTemp) + (WaterLev) + (BrineLev) + (WaterpH) + 
                  (WCBootLev) + (ECPdP) + (ECPLev) + (ReactorIP) + (ReactordP) + 
                  (ReactorCauLev) + (ClaydP) + (ClayOP) + KerFeed:ReactorIP +
                   ReactorIP:ClayOP + KerFeed:ClayOP + KerFeed:KerTemp +
                   KerFeed:WaterpH + KerTemp:ReactorIP + WaterLev:WaterpH, data=Model)
summary(DepFitTransInt)
par (mfrow=c(1,1))
plot(fitted(DepFitTransInt), resid(DepFitTransInt), col = "dodgerblue", pch = 20,
     xlab = "Fitted", ylab = "Residuals", main = "Transformed SLRM")
abline(h = 0, col = "darkorange", lwd = 2)

qqnorm(DepFitTransInt$residuals, col="dodgerblue", pch = 20, xlab="Theoretical Quantiles", ylab="Standardized Residuals")
qqline(DepFitTransInt$residuals,col="darkorange", lwd=2)

ModelStep = step(DepFitTransInt)

DepFitTransIntReduced = lm(log1p(DepRate) ~ D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + KerFeed + 
                             KerTemp + WaterLev + BrineLev + WaterpH + ECPdP + ECPLev + 
                             ReactorIP + ReactorCauLev + ClaydP + ClayOP + KerFeed:ReactorIP + 
                             ReactorIP:ClayOP + KerFeed:ClayOP + KerFeed:KerTemp + KerFeed:WaterpH + 
                             KerTemp:ReactorIP + WaterLev:WaterpH, data=Model)
summary(DepFitTransIntReduced)
g=DepFitTransIntReduced$coefficients
Predict1 = exp(predict(DepFitTransIntReduced,RealModel))
write.table(Predict1,"C:/Users/rpdlpascual/Desktop/Predict1.txt",sep="\t" )
write.table(g,"C:/Users/rpdlpascual/Desktop/coefficients.txt",sep="\t")

#NO INTERACTIONS
ModelStep2 = step(DepFitTrans)
DepFitTransReduced = lm(log1p(DepRate) ~ D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + KerFeed + 
                          KerTemp + WaterLev + WaterpH + ECPdP + ECPLev + ReactorIP + 
                          ReactorCauLev + ClaydP + ClayOP , data=Model)
summary(DepFitTransReduced)
j=DepFitTransReduced$coefficients
Predict2 = exp(predict(DepFitTransReduced, RealModel))
write.table(Predict2,"C:/Users/rpdlpascual/Desktop/Predict2.txt",sep="\t" )
write.table(j,"C:/Users/rpdlpascual/Desktop/coefficients2.txt",sep="\t")

#PCA
PCAfit=pcr(log1p(DepRate) ~., data=Model, scale=TRUE, validation="CV")
PCApred=predict(PCAfit,RealModel)
Prediction4=exp(PCApred)
write.table(Prediction4,"C:/Users/rpdlpascual/Desktop/prediction5.txt",sep="\t")

