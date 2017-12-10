# title: "IODS Final Assignment: Boston Dataset"
# author: "Tuure Parviainen"
# date: "December 5, 2017"

# As there is not that much manipulation to be done, for such used test dataset as Boston,
# I'll look into normality of assumptions in the data as often at least prices are log-normal,
# so perhaps some transformations are in order.
# I'll use shapiro test and fitdistrplus package to analyse the distribution of the each variable.

# LOG
#Try reloading and running again? MGCV might do something odd?

# https://m-clark.github.io/docs/GAM.html#generalized_additive_models
# PCA GAM
# http://www.win-vector.com/blog/2016/05/pcr_part2_yaware/
# https://github.com/WinVector/PreparingDataWorkshop/blob/master/YAwarePCA/YAwarePCA.Rmd
# Loading packages

require("easypackages") # for loading and installing many packages at one time
packages_to_load <- c("dplyr", "tidyverse", "ggplot2", "GGally", "MASS", "reshape2", "FactoMineR", "factoextra", "pls", "visreg", "mgcv", "ggpubr")
packages(packages_to_load, prompt = TRUE) # lock'n'load install/load

# rm(list = ls())
######################################################################################
#DataSet Loading and division to train and test sets
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # to set

# set the seed to make your partition reproductible
set.seed(777)

# 1) Normal unchanged dataset
boston <- read.table("data/boston.txt")

# 70% of the sample size
sample_size <- floor(0.70 * nrow(boston))

train_ind <- sample(seq_len(nrow(boston)), size = sample_size)
# splitting data by the vector
trainb <- boston[train_ind, ]
testb <- boston[-train_ind, ]

# 2) Scaled dataset
bostons <- read.table("data/bostons.txt")
# 70% of the sample size

train_inds <- sample(seq_len(nrow(bostons)), size = sample_size)
# splitting data by the vector
trainbs <- bostons[train_inds, ]
testbs <- bostons[-train_inds, ]


# 3) Scaled dataset with the crim transformed with the square transformation (log can't deal with zeros).
bostonssq <- read.table("data/bostonssq.txt")
# General note: It seems that some of the data has been censored for the high tax category
# 70% of the sample size

train_indssq <- sample(seq_len(nrow(bostonssq)), size = sample_size)
# splitting data by the vector
trainbssq <- bostonssq[train_indssq, ]
testbssq <- bostonssq[-train_indssq, ]
######################################################################################
# We are going to look at the nitrous oxide with the data available.

# Data exploration

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) +
    geom_point() +
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}
 g = ggpairs(boston, lower = list(continuous = my_fn))
 g
######################################################################################
# Normal dataset
######################################################################################

m <- lm(nox ~ ., data = trainb)

# Model is able to explain 78 %
stepAIC(m, direction = "both")

m2 <- lm(nox ~ crim + indus + age + dis + rad + ptratio + black + medv, data = trainb)

# Which is better m or m2 by AIC?

AIC(m,m2)
anova(m2, m)

# By AIC more parsimonious model is better.
plot(m2)

# Let's look at the data again with the same variables as in M2
g = ggpairs(trainb,columns = c(1,3,7:8,10,14), lower = list(continuous = my_fn))
g
summary(m2)
summary(m2)$r.squared
# We still have quite many dimensions in the data. Maybe we could lower the dimensions using PCA?
# R2 = 78 %

# Testing how the predicted values fit with the actual values
pred.test <- predict(m2, newdata = testb)

#simple
#plot(pred,testb$nox, xlab="predicted",ylab="actual")
#abline(a=0,b=1)

#Neat
lm.scatter <- data.frame(pred.test, testb$nox)
colnames(lm.scatter) <- c("pred.test","nox")
ggscatterhist(lm.scatter, x = "pred.test", y = "nox",
              color = "#00AFBB",
              margin.params = list(fill = "lightgray"))

#Calculating mean squared error
#pred.test <- predict(m2, newdata = testb)
mse.train <- summary(m2)$sigma^2
mse.test  <- sum((pred.test - testb$nox)^2)/(nrow(testb)-length(trainb)-2)

mse.train
mse.test
(1-mse.train/mse.test) *100
# Mean squared error ratio is 7.7 %
# This suggests that the new test data can be predicted with similar accuracy than the train data.



######################################################################################
# PCA

# Note that, by default, the function PCA() [in FactoMineR], standardizes the data automatically during the PCA; so you don’t need do this transformation before the PCA.

res.pca <- PCA(boston, graph = FALSE)
#summary(res.pca)

# How many dimensions we need
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
# Individuals and dimensions
fviz_pca_biplot(res.pca)
# Contributions of variables to PC1, PC2, PC3
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
#fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
#fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
#fviz_contrib(res.pca, choice = "var", axes = 1:3, top = 10)

#eig.val <- get_eigenvalue(res.pca)
#eig.val

# With eigenvalue higher than 1 are included. dim 1:3.

# Most influential variables.
fviz_pca_var(res.pca, select.var= list(cos2 = 5))

# With LM we had: crim + indus + age + dis + rad + ptratio + medv
# With PCA:  medv + indus + dis + lstat and nox
# The trouble is our data is still highly dimensional. What to do?


######################################################################################
# LM on PCA so PCR

pcr_model <- pcr(nox~., data = boston, scale = TRUE, validation = "CV")
#summary(pcr_model)

#validationplot(pcr_model)
#validationplot(pcr_model, val.type="MSEP")
validationplot(pcr_model, val.type = "R2")

# What you would like to see is a low cross validation error with a lower number of components than the number of variables in your dataset. If this is not the case or if the smalles cross validation error occurs with a number of components close to the number of variables in the original data, then no dimensionality reduction occurs. In the example above, it looks like 3 components are enough to explain more than 90% of the variability in the data although the CV score is a little higher than with 4 or 5 components. Finally, note that 6 components explain all the variability as expected.
predplot(pcr_model)
#coefplot(pcr_model)

#Test train
pcr_model <- pcr(nox~., data = trainb ,scale =TRUE, validation = "CV")

pcr_pred <- predict(pcr_model, testb, ncomp = 2)
mean((pcr_pred - testb$nox)^2)

# 0.2822
# PCR seems worse than linear regression.

plot(pcr_model, ncomp = 2) # Plot of cross-validated predictions
#plot(pcr_model, "scores") # Score plot
#plot(pcr_model, "loadings", comps = 1:3) # The three first loadings
#plot(pcr_model, "coef", ncomp = 5) # Coefficients
#plot(pcr_model, "val") # RMSEP curves
plot(pcr_model, "val", val.type = "MSEP", estimate = "CV") # CV MSEP

#library(ggpubr)
pcr_pred <- data.frame(pcr_pred)
colnames(pcr_pred) <- "pcr_pred"
scatter <- cbind(testb$nox, pcr_pred)
ggscatterhist(scatter, x = "pcr_pred", y = "testb$nox",
                             color = "#00AFBB",
                             margin.params = list(fill = "lightgray"))
ggscatterhist(
  scatter, x = "pcr_pred", y = "testb$nox",
  color = "blue", size = 3, alpha = 0.6,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.plot = "boxplot",
  ggtheme = theme_bw()
)


######################################################################################
# scaled dataset
######################################################################################

ms <- lm(nox ~ ., data = trainbs)

# Model is able to explain 78 %
stepAIC(m, direction = "both")

ms2 <- lm(nox ~ crim + indus + zn + age + dis + rad + rm + ptratio + black + medv, data = trainbs)

# Which is better m or m2 by AIC?

AIC(ms,ms2)
anova(ms2, ms)

# By AIC more parsimonious model is better.
par(mfrow=c(2,2))
plot(ms2)
dev.off()
# Let's look at the data again with the same variables as in M2
# g = ggpairs(trainb,columns = c(1,3,7:8,10,14), lower = list(continuous = my_fn))
# g
summary(ms2)
summary(ms2)$r.squared
# We still have quite many dimensions in the data. Maybe we could lower the dimensions using PCA?
# R2 = 76.5 % R-squared is worse.

# Testing how the predicted values fit with the actual values
pred.s <- predict(m2, newdata = trainbs)

plot(pred.s,trainbs$nox, xlab="predicted",ylab="actual")
abline(a=0,b=1)

#Calculating mean squared error
pred.test.s <- predict(m2, newdata = testbs)
mse.train.s <- summary(m2)$sigma^2
mse.test.s  <- sum((pred.test.s - testbs$nox)^2)/(nrow(testbs)-length(trainbs)-2)

mse.train.s
mse.test.s
(1-mse.train.s/mse.test.s)*100
# Mean squared error ratio is 11.3 % Poorer prediction
# This suggests that the new data can be predicted with good accuracy, scaling decreased peformance by 11 % on new data. Perhaps overfit?

######################################################################################
# Transformed dataset
######################################################################################

mssq <- lm(nox ~ ., data = trainbssq)

# Model is able to explain 78 %
stepAIC(mssq, direction = "both")

mssq2 <- lm(nox ~ indus + chas + age + dis + rad + ptratio + medv, data = trainbssq)

# Which is better m or m2 by AIC?

AIC(mssq,mssq2)
anova(mssq2, mssq , test="Chisq")

# By AIC and Chi's square test model 2 is better.
par(mfrow=c(2,2))
plot(mssq2)
dev.off()
# Let's look at the data again with the same variables as in M2
# g = ggpairs(trainb,columns = c(1,3,7:8,10,14), lower = list(continuous = my_fn))
# g
summary(mssq2)
summary(mssq2)$r.squared
# We still have quite many dimensions in the data. Maybe we could lower the dimensions using PCA?
# R2 = 78.9 %

# Testing how the predicted values fit with the actual values
pred.ssq <- predict(mssq2, newdata = trainbssq)

plot(pred.ssq,trainbssq$nox, xlab="predicted",ylab="actual")
abline(a=0,b=1)

#Calculating mean squared error
pred.test.ssq <- predict(mssq2, newdata = testbssq)
mse.train.ssq <- summary(mssq2)$sigma^2
mse.test.ssq  <- sum((pred.test.ssq - testbssq$nox)^2)/(nrow(testbssq)-length(trainbssq)-2)

mse.train.ssq
mse.test.ssq
(1 - mse.train.ssq/mse.test.ssq) *100
# Mean squared error ratio is 12.57 %
# This suggests that the new data can be predicted with even better accuracy compared to scaled linear regression.


######################################################################################
# LM on PCA so PCR - transformed dataset

pcr_model.ssq <- pcr(nox~., data = bostonssq, scale = TRUE, validation = "CV")
#summary(pcr_model.ssq)

#validationplot(pcr_model.ssq)
#validationplot(pcr_model.ssq, val.type="MSEP")
validationplot(pcr_model.ssq, val.type = "R2")

# What you would like to see is a low cross validation error with a lower number of components than the number of variables in your dataset.
# If this is not the case or if the smalles cross validation error occurs with a number of components close to the number of variables in the original data,
# then no dimensionality reduction occurs. In the example above, it looks like 3 components are enough to explain more than 90% of the variability in the
# data although the CV score is a little higher than with 4 or 5 components. Finally, note that 6 components explain all the variability as expected.
predplot(pcr_model.ssq)
#coefplot(pcr_model.ssq)

#Test train
pcr_model.ssq <- pcr(nox~., data = trainbssq ,scale =TRUE, validation = "CV")

pcr_pred.ssq <- predict(pcr_model.ssq, testbssq, ncomp = 2)
mean((pcr_pred.ssq - testbssq$nox)^2)

# normal PCR MSE      0.2822
# transformed PCR MSE 0.2514871
# transformed PCR gives better MSE.


plot(pcr_model.ssq, ncomp = 2) # Plot of cross-validated predictions
#plot(pcr_model, "scores") # Score plot
#plot(pcr_model, "loadings", comps = 1:3) # The three first loadings
#plot(pcr_model, "coef", ncomp = 5) # Coefficients
#plot(pcr_model, "val") # RMSEP curves
plot(pcr_model.ssq, "val", val.type = "MSEP", estimate = "CV") # CV MSEP

library(ggpubr)
pcr_pred.ssq <- data.frame(pcr_pred.ssq)
colnames(pcr_pred.ssq) <- "pcr_pred.ssq"
scatter <- cbind(testbssq$nox, pcr_pred.ssq)
ggscatterhist(scatter, x = "pcr_pred.ssq", y = "testbssq$nox",
              color = "#00AFBB",
              margin.params = list(fill = "lightgray"))
# ggscatterhist(
#   scatter, x = "pcr_pred", y = "testb$nox",
#   color = "blue", size = 3, alpha = 0.6,
#   palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#   margin.plot = "boxplot",
#   ggtheme = theme_bw()
# )

######################################################################################
# Y-aware PCA
######################################################################################
# x’ such that a unit change in x’ corresponds to a unit change in y. Under this rescaling, all the independent variables are in the same units,
# which are indeed the natural units for the problem at hand: characterizing their effect on y. (We also center the transformed variables x’ to be zero mean,
# as is done with standard centering and scaling).

# Using guide from http://www.win-vector.com/blog/2016/05/pcr_part2_yaware/
# https://github.com/WinVector/PreparingDataWorkshop/blob/master/YAwarePCA/YAwarePCA.Rmd

# vtreat package
library(vtreat)
library('WVPlots') # devtools::install_github('WinVector/WVPlots',build_vignettes=TRUE)

# design treatment plan
treatmentsN <- designTreatmentsN(trainb,setdiff(colnames(trainb),'nox'),'nox',
                                 verbose=FALSE)

scoreFrame = treatmentsN$scoreFrame
scoreFrame$vartype = ifelse(grepl("noise", scoreFrame$varName), "noise", "signal")

######################################################################################
#Defining functions
dotplot_identity = function(frame, xvar, yvar, colorvar=NULL) {
  if(is.null(colorvar)) {
    gplot = ggplot(frame, aes_string(x=xvar, y=yvar, ymax=yvar))
  } else {
    gplot = ggplot(frame, aes_string(x=xvar, y=yvar, ymax=yvar, color=colorvar))
  }

  gplot + geom_point() + geom_linerange(aes(ymin=0))
}

barbell_plot = function(frame, xvar, ymin, ymax, colorvar=NULL) {
  if(is.null(colorvar)) {
    gplot = ggplot(frame, aes_string(x=xvar))
  } else {
    gplot = ggplot(frame, aes_string(x=xvar, color=colorvar))
  }

  gplot + geom_point(aes_string(y=ymin)) +
    geom_point(aes_string(y=ymax)) +
    geom_linerange(aes_string(ymin=ymin, ymax=ymax))
}

extractProjection <- function(ndim,princ) {
  # pull off the rotation.
  proj <- princ$rotation[,1:ndim]
  # sign was arbitrary, so flip in convenient form
  for(i in seq_len(ndim)) {
    si <- sign(mean(proj[,i]))
    if(si!=0) {
      proj[,i] <- proj[,i]*si
    }
  }
  proj
}

rsq <- function(x,y) {
  1 - sum((y-x)^2)/sum((y-mean(y))^2)
}

######################################################################################
# The high sig means insignificant variables, remove chas?
# dotplot_identity(scoreFrame, "varName", "sig", "vartype") +
#   coord_flip()  + ggtitle("vtreat variable significance estimates")+
#   scale_color_manual(values = c("noise" = "#d95f02", "signal" = "#1b9e77"))

# prepare the treated frames, with y-aware scaling
examplePruneSig = 1
trainbNTreatedYScaled <- prepare(treatmentsN,trainb,pruneSig=examplePruneSig,scale=TRUE)
testbNTreatedYScaled <- prepare(treatmentsN,trainb,pruneSig=examplePruneSig,scale=TRUE)

# get the variable ranges
ranges = vapply(trainbNTreatedYScaled, FUN=function(col) c(min(col), max(col)), numeric(2))
rownames(ranges) = c("vmin", "vmax")
rframe = as.data.frame(t(ranges))  # make ymin/ymax the columns
rframe$varName = rownames(rframe)
varnames = setdiff(rownames(rframe), "nox")
rframe = rframe[varnames,]
rframe$vartype = ifelse(grepl("noise", rframe$varName), "noise", "signal")

names(trainbNTreatedYScaled )
# show a few columns
#summary(trainbNTreatedYScaled)
summary(trainbNTreatedYScaled[, c("nox", "rad_clean", "indus_clean", "ptratio_clean", "medv_clean")])

# barbell_plot(rframe, "varName", "vmin", "vmax", "vartype") +
#   coord_flip() + ggtitle("y-scaled variables: ranges") +
#   scale_color_manual(values = c("noise" = "#d95f02", "signal" = "#1b9e77"))

vars <- setdiff(colnames(trainbNTreatedYScaled),'nox')
# prcomp defaults to scale. = FALSE, but we already scaled/centered in vtreat- which we don't want to lose.
dmTrain <- as.matrix(trainbNTreatedYScaled[,vars])
dmTest <- as.matrix(trainbNTreatedYScaled[,vars])
princ <- prcomp(dmTrain, center = FALSE, scale. = FALSE)
dotplot_identity(frame = data.frame(pc=1:length(princ$sdev),
                                    magnitude=princ$sdev),
                 xvar="pc",yvar="magnitude") +
  ggtitle("Y-Scaled variables: Magnitudes of singular values")

proj <- extractProjection(3,princ) # taking the first 3 principal components
rot5 <- extractProjection(5,princ) # taking the first 3 principal components
rotf = as.data.frame(rot5)
rotf$varName = rownames(rotf)
rotflong = gather(rotf, "PC", "loading", starts_with("PC"))
rotflong$vartype = ifelse(grepl("noise", rotflong$varName), "noise", "signal")

dotplot_identity(rotflong, "varName", "loading", "vartype") +
  facet_wrap(~PC,nrow=1) + coord_flip() +
  ggtitle("Y-Scaled Variable loadings, first five principal components") +
  scale_color_manual(values = c("noise" = "#d95f02", "signal" = "#1b9e77"))

# apply projection
projectedTrain <- as.data.frame(dmTrain %*% proj,
                                stringsAsFactors = FALSE)
# plot data sorted by principal components
projectedTrain$nox <- dTrainNTreatedYScaled$nox
ScatterHistN(projectedTrain,'PC1','PC2','nox',
             "Y-Scaled Training Data projected to first two principal components")

model <- lm(nox~PC1+PC2+PC3,data=projectedTrain)
model2 <- lm(nox~PC1+PC2,data=projectedTrain)
AIC(model,model2)
summary(model2)

#WHY IS THE MODEL SUDDENLY SO ODD? TRACEBACK!

dev.off()
plot(model, 2)


projectedTrain$estimate <- predict(model,newdata=projectedTrain)
trainrsq = rsq(projectedTrain$estimate,projectedTrain$nox)

ScatterHist(projectedTrain,'estimate','nox','Recovered model versus truth (y aware PCA train)',
            smoothmethod='identity',annot_size=3)

# apply projection
projectedTest <- as.data.frame(dmTest %*% proj,
                               stringsAsFactors = FALSE)
# plot data sorted by principal components
projectedTest$nox <- dTestNTreatedYScaled$nox
ScatterHistN(projectedTest,'PC1','PC2','nox',
             "Y-Scaled Test Data projected to first two principal components")

projectedTest$estimate <- predict(model,newdata=projectedTest)
testrsq = rsq(projectedTest$estimate,projectedTest$nox)
testrsq

ScatterHist(projectedTest,'estimate','nox','Recovered model versus truth (y aware PCA test)',
            smoothmethod='identity',annot_size=3)

#gam
#library(mgcv)


gmodel <- gam(nox ~ s(PC1) + PC2 + PC3, family  = gaussian, data=projectedTrain,
    trace=FALSE)

gmodel2 <- gam(nox ~ s(PC1) + PC2 + PC3, family  = gaussian, data=projectedTrain,
              trace=FALSE)
gmodel3 <- gam(nox ~ s(PC1, bs="cr") + PC2, data=projectedTrain)
gmodel4 <- gam(nox ~ s(PC1, bs="cr") + PC2 + PC3, data=projectedTrain)
gmodel5 <- gam(nox ~ te(PC1, PC2),method="REML", data=projectedTrain)
gmodel6 <- gam(nox ~ te(PC1, PC2) + PC3, data=projectedTrain)

# te =  tensor product smooth, and by smoothing the marginal smooths of PC1 and PC2

AIC(gmodel,gmodel2,gmodel3,gmodel4,gmodel5,gmodel6)

#par(mfrow=c(1,3)) #to partition the Plotting Window
plot(gmodel,se = TRUE)
plot(gmodel4,se = TRUE)
#se stands for standard error Bands

summary(gmodel2)$r.sq
summary(gmodel)$sp.criterion
anova(gmodel,gmodel2,gmodel3,gmodel4,gmodel5,gmodel6, test="Chisq")

vis.gam(gmodel6, type='response', plot.type='contour')
visreg2d(gmodel6, xvar='PC1', yvar='nox', scale='response')
visreg2d(gmodel6, xvar='PC2', yvar='nox', scale='response')

# Lets confirm out model with replications of K (the smoother term)
par(mfrow=c(1,3))
gam.check(gmodel5, k.rep=1000, type=c("pearson"))
qq.gam(gmodel5)
plot(gmodel5,pages=1)

gam.check(gmodel6, k.rep=1000, type=c("pearson"))

#visreg2d(gmodel6, xvar='PC3', yvar='nox', scale='response')

#We have perfectly divided are variance to two dimensions, but we can only make use of this with similar data.
# How about predictions on our test data?

projectedTest$g.estimate <- predict(gmodel5, newdata = projectedTest)

trainrsq_gam = rsq(projectedTest$g.estimate ,projectedTrain$nox)

ScatterHist(projectedTest,'g.estimate','nox','Recovered model versus truth (y aware PCA test) GAM ',
            smoothmethod='identity',annot_size=3)

#### WITH FACTOMINER

res.pcaY <- PCA(dmTrain,scale.unit = FALSE, graph = FALSE)
# How many dimensions we need
fviz_eig(res.pcaY, addlabels = TRUE, ylim = c(0, 50))
# Individuals and dimensions
fviz_pca_biplot(res.pcaY)
# Contributions of variables to PC1, PC2, PC3
fviz_contrib(res.pcaY, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pcaY, choice = "var", axes = 2, top = 10)
fviz_contrib(res.pcaY, choice = "var", axes = 3, top = 10)
fviz_contrib(res.pcaY, choice = "var", axes = 1:2, top = 10)

eig.val <- get_eigenvalue(res.pcaY)
eig.val


