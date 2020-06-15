Model = read.csv ("DepModel.csv")
DepFit = lm(DepRate ~ KerFeed + KerTemp + WaterpH + WaterLev + BrineLev, data = Model)
summary (DepFit)

#model adequacy checking
#check for normality
par (mfrow=c(2,2),mar=c(2,2,2,2))
plot (DepFit)

plot(hatvalues(DepFit))
abline(h=7/1254, col="red")

plot(rownames(Model),rstudent(DepFit))
abline(h=3,col="red")
rstudent(DepFit)

DepFit2 = lm(DepRate ~ KerFeed + KerTemp + WaterLev + BrineLev, data = Model)
summary(DepFit2)
ResM = mean(DepFit2[["residuals"]])
print (ResM)

Stepwise = step (DepFit, scope = list(lower=~1, upper=~KerFeed + KerTemp + WaterpH + WaterLev + BrineLev, direction="both",trace=1))
