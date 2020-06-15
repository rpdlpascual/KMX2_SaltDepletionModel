#Loading necessary libraries
library(lattice)
library(ggplot2)
library(caret)

DepModelCV=read.csv("DepModelCV.csv")
set.seed(123)
trctrl = trainControl(method="repeatedcv", number=10,repeats="10")
Model = train(DepRate~.,data=DepModelCV,method="lm",trControl=trctrl)
print(Model)
Model$finalModel

ModelZ=lm(DepRate~.,data=DepModelCV)
summary(ModelZ)
