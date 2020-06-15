#Installation of necessary libraries
install.packages("ggplot2")
install.packages("lattice")
install.packages("psych")
library(ggplot2)
library(lattice)
library(caret)
library(psych)

#read the data
data = sat.act
head(data)

#set the number of splits, generally, the fewer the number of splits, the lower the variance and the higher the bias/error and vice versa
datactrl = trainControl(method="cv",number=5)
model = train(ACT~gender+age+SATV+SATQ, data=data, trControl=datactrl, method="lm",na.action = na.pass)
model

model$finalModel
model$resample

data2 = swiss
set.seed(123)
training = swiss$Fertility, createDataPartition=(p=0.8, list=FALSE)

model3=lm(Fertility~., data=data2)
predictions=model3,predict(data2)
