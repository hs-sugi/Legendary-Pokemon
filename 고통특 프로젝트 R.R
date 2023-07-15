### 고통특 프로젝트
# 222STG10 김희숙

### 데이터 출처 : https://www.kaggle.com/datasets/alopez247/pokemon
library(tidyverse)

df <- read.csv(file="/Users/ssugi/Downloads/pokemon_alopez247.csv", header=TRUE)
dim(df) # 721  23
head(df)
str(df)
df$Name <- as.factor(df$Name)
df$Type_1 <- as.factor(df$Type_1)
df$Type_2 <- as.factor(df$Type_2)
df$isLegendary <- as.factor(df$isLegendary)
df$Color <- as.factor(df$Color)
df$hasGender <- as.factor(df$hasGender)
df$Egg_Group_1 <- as.factor(df$Egg_Group_1)
df$Egg_Group_2 <- as.factor(df$Egg_Group_2)
df$hasMegaEvolution <- as.factor(df$hasMegaEvolution)
df$Body_Style <- as.factor(df$Body_Style)

summary(df) # 675 46 # isLegendary가 15:1 정도
df[is.na(df$Pr_Male),]

library(plyr)
library(psych)
multi.hist(df) #error, not numeric
multi.hist(df[,sapply(df, is.numeric)], global = FALSE)

# Pr_Male -> NA: gender가 없는 경우
df[is.na(df$Pr_Male),]
which(is.na(df$Pr_Male))
sum(is.na(df$Pr_Male)) #77개
df[is.na(df)]
df$Pr_Male[is.na(df$Pr_Male)] <- 1.5 # 결측값 채우기
df[is.na(df)] #없음

which(df$Height_m>summary(df$Height_m)[5] + 1.5*IQR(df$Height_m)) # Height_m에서 이상값
which(df$Weight_kg>summary(df$Weight_kg)[5] + 1.5*IQR(df$Weight_kg)) # Weight_kg에서 이상값
intersect(which(df$Height_m>summary(df$Height_m)[5] + 1.5*IQR(df$Height_m)),which(df$Weight_kg>summary(df$Weight_kg)[5] + 1.5*IQR(df$Weight_kg))) #공통 값 추출
# 95 130 208 249 250 321 350 382 383 384 483 484 486 487 493 614 623 643 644 646 699 716 717 718
which(df$Catch_Rate>summary(df$Catch_Rate)[5] + 1.5*IQR(df$Catch_Rate)) # 이상치 없음
outInd <- intersect(which(df$Height_m>summary(df$Height_m)[5] + 1.5*IQR(df$Height_m)),which(df$Weight_kg>summary(df$Weight_kg)[5] + 1.5*IQR(df$Weight_kg)))

length(outInd) # 24

df2 <- df[-outInd,]
dim(df2) # 697  23


# 범주형 drop, classification을 진행할 column은 label encoding
pokemon <- df2[, -c(1,2,3,4,12,14,15,17,18,19,23)] # categorical drop

head(pokemon)
dim(pokemon) # 693  12

summary(pokemon) # 667 30 # isLegendary가 22:1 정도
pokemon$isLegendary <- as.numeric(factor(pokemon$isLegendary))

pairs(pokemon[,-8], main = "scatter plot",pch = 21, bg = c("red", "blue")[unclass(pokemon$isLegendary)])



#install.packages("corrplot")
library(corrplot) 
pokemon_cor <- cor(pokemon, use = "complete.obs")
corrplot(pokemon_cor, method="circle")



# 오버샘플링
# https://rpubs.com/yoompubs/467234
# Majority Weighted Minority Oversampling TEchnique
# MWMOTE is an extension of the original SMOTE algorithm. It assigns higher weight to borderline instances, undersized minority clusters and examples near the borderline of the two classes.
library(imbalance)
newMWMOTE <- mwmote(pokemon, numInstances = 600, classAttr = "isLegendary")
dim(newMWMOTE) # 600  12
newMWMOTE$isLegendary <- as.factor(newMWMOTE$isLegendary)
pokemon$isLegendary <- as.factor(pokemon$isLegendary)
plotComparison(pokemon, rbind(pokemon, newMWMOTE), attrs = names(pokemon)[1:3], classAttr = "isLegendary")
plotComparison(pokemon, rbind(pokemon, newMWMOTE), attrs = names(pokemon)[4:6], classAttr = "isLegendary")
newdf <- rbind(pokemon, newMWMOTE)
summary(newdf) # 667 vs 630
newdf <- data.frame(newdf)

# x, y 정의
y <- newdf[, 8]
x <- newdf[,-8]
pairs(x,main = "scatter plot - original")

# boxcox
n <- dim(x)[1] #1293
p <- dim(x)[2] #11
lam <- rep(1,p)
eps <- 0.5
lamopt <- gauss(x,lam,.5,20)
xlam <- bocotranmat(x,lamopt,eps)
x_xlam <- data.frame(xlam)
names(x_xlam) <- c("Total",  "HP", "Attack", "Defense", "Sp_Atk", "Sp_Def", "Speed", "Pr_Male", "Height_m", "Weight_kg", "Catch_Rate")
pairs(x_xlam,main = "scatter plot - boxcox") # elliptical shape

newPokemon <- cbind(y,x_xlam) # 최종 사용 데이터


# ladle
par(mfrow=c(1,3)) 
y <- as.numeric(y)

# sir plot
n <- dim(x_xlam)[1]; p=11; nboot=200; method="sir";ytype="categorical"; h=2
out=ladle(x_xlam, as.numeric(y),h,nboot,method,ytype); 
kset=out$kset;gn=out$gn 
plot(kset,gn,type="l",xlab="k",main="Ladle - SIR")


# save plot
n <- dim(x_xlam)[1] ; p=11; nboot=200; method="save";ytype="categorical"; h=2
out=ladle(x_xlam,y,h,nboot,method,ytype);
kset=out$kset;gn=out$gn 
plot(kset,gn,type="l",xlab="k",main="Ladle - SAVE")


# dr plot
n <- dim(x_xlam)[1] ; p=11; nboot=200; method="dr"; ytype="categorical" 
out=ladle(x_xlam,y,h,nboot,method,ytype);
kset=out$kset;gn=out$gn 
plot(kset,gn,type="l",xlab="k",main="Ladle - DR")





# model&accuracy

# train_test_split
library(creditmodel)
library(e1071)
train_test = train_test_split(newPokemon, prop = 0.7,
                              seed = 1004, save_data = FALSE)

xy.tra = train_test$train
dim(train_test$train) # 908  12
xy.tra.12 = rbind(xy.tra[xy.tra[,1]==1,], # 473
                  xy.tra[xy.tra[,1]==2,]) # 435
x.tra=xy.tra.12[,2:12];
y.tra=as.matrix(xy.tra.12[,1])

xy.tes = train_test$test
dim(train_test$test) # 389  12
xy.tes.12 = rbind(xy.tes[xy.tes[,1]==1,], # 194
                  xy.tes[xy.tes[,1]==2,]) # 195
x.tes=xy.tes.12[,2:12];
y.tes=xy.tes.12[,1]



#################################
#            Pure               #
#################################
# SVM
### pure
svm_pure<- svm(factor(y) ~ ., data=xy.tra ) 
yhat_test <- predict(svm_pure, x.tes)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9922879


# random forest
library(randomForest)
### pure
rf_pure<- randomForest(factor(y) ~ ., data=xy.tra ) 
yhat_test <- predict(rf_pure, x.tes)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9922879



# xgboost
library(xgboost)
### pure
dtrain <- xgb.DMatrix(data=as.matrix(x.tra), label=as.matrix(y.tra))
dtest <- xgb.DMatrix(data=as.matrix(x.tes), label=as.matrix(y.tes))
xgb_pure <- xgb.train(data=dtrain,nrounds=50)
yhat_test <- round(predict(xgb_pure, as.matrix(x.tes))) 
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9820051


# naive bayes
library(naivebayes)
### pure
nb_pure<- naive_bayes(factor(y) ~ ., data=xy.tra ) 
yhat_test <- predict(nb_pure, x.tes)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9794344
##################################################################




#################################
#   Dimension reduction        #
#################################

# ladle을 통해 찾은 값 넣기 #################################
# SIR #d=1
# 전체 데이터
h=2;r=1;ytype="categorical";sir_beta=sir(x_xlam,y,h,r,ytype) 
x_xlam <- as.matrix(x_xlam)
sir_pred=center(x_xlam)%*%sir_beta 
par(mfrow=c(1,1)) 
plot(sir_pred,y,xlab="first SIR predictor",ylab="y", main="SIR Plot")
points(sir_pred[y==1,1],y[y==1],col="red")
points(sir_pred[y==2,1],y[y==2],col="blue")


#################################

# svm
# train
sir_beta = sir(x.tra,y.tra,h,r,ytype)
sir_pred <-  as.matrix(center(x.tra))%*%sir_beta # pred=center(x)%*%beta
input.data <- cbind(y.tra, sir_pred)
colnames(input.data) <- c('S1','S2')
input.data <- data.frame(input.data)
input.data$S1 <- as.integer(input.data$S1)
input.data$S2 <- as.integer(input.data$S2)

svm_sir<- svm(factor(S1) ~ ., data=input.data)

# test
pred_test <- as.matrix(center(x.tes))%*%sir_beta 
colnames(pred_test) <- c('S2')
yhat_test <- predict(svm_sir,pred_test)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.966581


# random forest
rf_sir<- randomForest(factor(S1) ~ ., data=input.data)
yhat_test <- predict(rf_sir,pred_test)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9640103



# xgboost
dtrain <- xgb.DMatrix(data=as.matrix(sir_pred), label=as.matrix(y.tra))
dtest <- xgb.DMatrix(data=as.matrix(pred_test), label=as.matrix(y.tes))
xgb_sir <- xgb.train(data=dtrain,nrounds=50)
yhat_test <- round(predict(xgb_sir, as.matrix(pred_test))) 
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9614396


# naive bayes
nb_sir<- naive_bayes(factor(S1) ~ ., data=input.data)
yhat_test <- predict(nb_sir, pred_test)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9717224





# SAVE #d=3
# 전체 데이터
h=2;r=3;ytype="categorical";save_beta=save(x_xlam,y,h,r,ytype) 
save_pred=center(x_xlam)%*%save_beta 
plot_data_save <- data.frame(cbind(save_pred[,1],save_pred[,2],save_pred[,3],y))
colnames(plot_data_save) <- c('1PC', '2PC', '3PC', 'Y')
pairs(plot_data_save, pch = 21, bg = c("red", "blue")[unclass(plot_data_save$Y)], main="SAVE Plot")

#################################

# svm
# train
save_beta=save(x.tra,y.tra,h,r,ytype) 
save_pred <-  as.matrix(center(x.tra))%*%save_beta # pred=center(x)%*%beta
input.data <- cbind(y.tra, save_pred)
colnames(input.data) <- c('S1','S2','S3','S4')
input.data <- data.frame(input.data)
input.data$S1 <- as.integer(input.data$S1)
input.data$S2 <- as.integer(input.data$S2)
input.data$S3 <- as.integer(input.data$S3)
input.data$S4 <- as.integer(input.data$S4)

svm_save<- svm(factor(S1) ~ ., data=input.data)

# test
pred_test <- as.matrix(center(x.tes))%*%save_beta 
colnames(pred_test) <- c('S2','S3','S4')
yhat_test <- predict(svm_save,pred_test)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) # 0.9794344


# random forest
rf_save<- randomForest(factor(S1) ~ ., data=input.data)
yhat_test <- predict(rf_save,pred_test)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9717224



# xgboost
dtrain <- xgb.DMatrix(data=as.matrix(save_pred), label=as.matrix(y.tra))
dtest <- xgb.DMatrix(data=as.matrix(pred_test), label=as.matrix(y.tes))
xgb_save <- xgb.train(data=dtrain,nrounds=50)
yhat_test <- round(predict(xgb_save, as.matrix(pred_test))) 
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9794344


# naive bayes
nb_save<- naive_bayes(factor(S1) ~ ., data=input.data)
yhat_test <- predict(nb_save, pred_test)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.6632391




# DR #d=3
# 전체 데이터
h=2;r=3;ytype="categorical"; dr_beta=dr(x_xlam,y,h,r,ytype) 
dr_pred=center(x_xlam)%*%dr_beta 
par(mfrow=c(1,1)) 
plot_data_dr <- data.frame(cbind(dr_pred[,1],dr_pred[,2],dr_pred[,3],y))
colnames(plot_data_dr) <- c('1PC', '2PC', '3PC', 'Y')
pairs(plot_data_dr, pch = 21, bg = c("red", "blue")[unclass(plot_data_save$Y)], main="DR Plot")

#################################

# svm
# train
dr_beta=dr(x.tra,y.tra,h,r,ytype)
dr_pred <-  as.matrix(center(x.tra))%*%dr_beta # pred=center(x)%*%beta
input.data <- cbind(y.tra, dr_pred)
colnames(input.data) <- c('S1','S2','S3','S4')
input.data <- data.frame(input.data)
input.data$S1 <- as.integer(input.data$S1)
input.data$S2 <- as.integer(input.data$S2)
input.data$S3 <- as.integer(input.data$S3)
input.data$S4 <- as.integer(input.data$S4)

svm_dr<- svm(factor(S1) ~ ., data=input.data)

# test
pred_test <- as.matrix(center(x.tes))%*%dr_beta 
colnames(pred_test) <- c('S2','S3','S4')
yhat_test <- predict(svm_dr,pred_test)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) # 0.907455


# random forest
rf_dr<- randomForest(factor(S1) ~ ., data=input.data)
yhat_test <- predict(rf_dr,pred_test)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9717224



# xgboost
dtrain <- xgb.DMatrix(data=as.matrix(dr_pred), label=as.matrix(y.tra))
dtest <- xgb.DMatrix(data=as.matrix(pred_test), label=as.matrix(y.tes))
xgb_dr <- xgb.train(data=dtrain,nrounds=50)
yhat_test <- round(predict(xgb_dr, as.matrix(pred_test))) 
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.9691517


# naive bayes
nb_dr<- naive_bayes(factor(S1) ~ ., data=input.data)
yhat_test <- predict(nb_dr, pred_test)
table <- table(true=y.tes, predict=yhat_test)
table
(table[1,1] + table[2,2]) / sum(table) #정확도: 0.7789203

###################################

###### 번외. total 변수 값이 큰 포켓몬 확인 ######
pokemon_t <- df2[, -c(1,3,4,12,14,15,17,18,19,23)]

df3 <- pokemon_t[order(pokemon_t$Total), ]
View(df3)
totalBig <- df3[649:693,c(1,2,9)] # 45
nrow(totalBig[totalBig$isLegendary == 1,]) # 15
nrow(totalBig[totalBig$isLegendary == 2,]) # 30
plot()
totalBig_1 <- totalBig[totalBig$isLegendary == 1,]
totalBig_2 <- totalBig[totalBig$isLegendary == 2,]
par(mfrow=c(1,2)) 
boxplot(totalBig_1$Total)
boxplot(totalBig_2$Total) #비슷비슷~
