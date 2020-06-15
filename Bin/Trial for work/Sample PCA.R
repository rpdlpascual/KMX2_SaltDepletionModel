#Ron Exploring PCA with so much stress lol
train = read.csv("trainRon.csv")
test = read.csv("testRon.csv")
test$Item_Outlet_Sales = 1
combi = rbind(train,test)
combi$Item_Weight[is.na(combi$Item_Weight)] = median(combi$Item_Weight, na.rm = TRUE)
combi$Item_Visibility = ifelse(combi$Item_Visibility == 0, median(combi$Item_Visibility), combi$Item_Visibility)                                
table(combi$Outlet_Size, combi$Outlet_Type)
levels(combi$Outlet_Size)[1] <- "Other"
my_data = subset(combi, select = -c(Item_Outlet_Sales, Item_Identifier, Outlet_Identifier))
colnames(my_data)
str(my_data)

library(dummies)                                 
new_my_data <- dummy.data.frame(my_data, names = c("Item_Fat_Content","Item_Type",
                                                   "Outlet_Establishment_Year","Outlet_Size",
                                                   "Outlet_Location_Type","Outlet_Type"))
str(new_my_data)
pca.train <- new_my_data[1:nrow(train),]
pca.test <- new_my_data[-(1:nrow(train)),]

prin_comp <- prcomp(pca.train, scale. = T)
names(prin_comp)

#outputs the mean of variables
prin_comp$center

#outputs the standard deviation of variables
prin_comp$scale

prin_comp$rotation
prin_comp$rotation[1:5,1:4]

dim(prin_comp$x)

biplot(prin_comp, scale = 0)

#compute standard deviation of each principal component
std_dev <- prin_comp$sdev

#compute variance
pr_var <- std_dev^2

#check variance of first 10 components
pr_var[1:10]

#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]

#scree plot
plot(prop_varex, xlab = "Principal Component",
       ylab = "Proportion of Variance Explained",
       type = "b")

#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
       ylab = "Cumulative Proportion of Variance Explained",
       type = "b")

#add a training set with principal components
train.data <- data.frame(Item_Outlet_Sales = train$Item_Outlet_Sales, prin_comp$x)

#we are interested in first 30 PCAs
train.data <- train.data[,1:31]

#run a decision tree
library(rpart)
rpart.model <- rpart(Item_Outlet_Sales ~ .,data = train.data, method = "anova")
rpart.model
rpart.model$

#transform test into PCA
test.data <- predict(prin_comp, newdata = pca.test)
test.data <- as.data.frame(test.data)

#select the first 30 components
test.data <- test.data[,1:30]

#make prediction on test data
rpart.prediction <- predict(rpart.model, test.data)

library(pls)
modell = pcr(Item_Outlet_Sales~., data = train.data, validation="CV")
modell$coefficients
modell.predict = predict(modell,test.data)
