# title: "IODS Final Assignment: Boston Dataset"
# author: "Tuure Parviainen"
# date: "December 5, 2017"

# As there is not that much manipulation to be done, for such used test dataset as Boston, I'll look into normality of assumptions in the data as often at least prices are log-normal, so perhaps some transformations are in order. I'll use shapiro test and fitdistrplus package to analyse the distribution of the each variable.

# Loading packages

require("easypackages") # for loading and installing many packages at one time
packages_to_load <- c("dplyr", "tidyverse", "ggplot2", "GGally", "MASS", "reshape2", "FactoMineR", "factoextra", "pls")
packages(packages_to_load, prompt = TRUE) # lock'n'load install/load

# rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # to set


# 1) Normal unchanged dataset
boston <- read.table("data/boston.txt")

# 70% of the sample size
sample_size <- floor(0.70 * nrow(boston))

# set the seed to make your partition reproductible
set.seed(777)
train_ind <- sample(seq_len(nrow(boston)), size = sample_size)
# splitting data by the vector
trainb <- boston[train_ind, ]
testb <- boston[-train_ind, ]

# 2) Scaled dataset
bostons <- read.table("data/bostons.txt")
# 3) Scaled dataset with the crim transformed with the square transformation (log can't deal with zeros).
bostonssq <- read.table("data/bostonssq.txt")
# General note: It seems that some of the data has been censored for the high tax category

# We are going to look at the nitrous oxide with the data available.

# Data exploration

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) +
    geom_point() +
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}
# g = ggpairs(boston, lower = list(continuous = my_fn))
# g
# columns = c(2:7)
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
pred <- predict(m2, newdata = trainb)

plot(pred,trainb$nox, xlab="predicted",ylab="actual")
abline(a=0,b=1)

#Calculating mean squared error
pred.test <- predict(m2, newdata = testb)
mse.train <- summary(m2)$sigma^2
mse.test  <- sum((pred.test - testb$nox)^2)/(nrow(testb)-length(trainb)-2)

mse.train
mse.test
mse.train/mse.test
# Mean squared error ratio is 92 %
# This suggests that the new data can be predicted with good accuracy.

######################################################################################
# PCA

# Note that, by default, the function PCA() [in FactoMineR], standardizes the data automatically during the PCA; so you don’t need do this transformation before the PCA.

res.pca <- PCA(boston, graph = FALSE)
# How many dimensions we need
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
# Individuals and dimensions
fviz_pca_biplot(res.pca)
# Contributions of variables to PC1, PC2, PC3
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 1:3, top = 10)

eig.val <- get_eigenvalue(res.pca)
eig.val

# With eigenvalue higher than 1 are included. dim 1:3.

# Most influential variables.
fviz_pca_var(res.pca, select.var= list(cos2 = 5))

# With LM we had: crim + indus + age + dis + rad + ptratio + medv
# With PCA:  medv + indus + dis + lstat and nox

#http://www.win-vector.com/blog/2016/05/pcr_part1_xonly/

######################################################################################
# LM on PCA so PCR

pcr_model <- pcr(nox~., data = boston, scale = TRUE, validation = "CV")
summary(pcr_model)

validationplot(pcr_model)
validationplot(pcr_model, val.type="MSEP")
validationplot(pcr_model, val.type = "R2")

# What you would like to see is a low cross validation error with a lower number of components than the number of variables in your dataset. If this is not the case or if the smalles cross validation error occurs with a number of components close to the number of variables in the original data, then no dimensionality reduction occurs. In the example above, it looks like 3 components are enough to explain more than 90% of the variability in the data although the CV score is a little higher than with 4 or 5 components. Finally, note that 6 components explain all the variability as expected.
predplot(pcr_model)
coefplot(pcr_model)

#Test train
pcr_model <- pcr(nox~., data = trainb ,scale =TRUE, validation = "CV")

pcr_pred <- predict(pcr_model, testb, ncomp = 2)
mean((pcr_pred - testb$nox)^2)

# PCR seems worse than linear regression.

plot(pcr_model, ncomp = 5) # Plot of cross-validated predictions
plot(pcr_model, "scores") # Score plot
plot(pcr_model, "loadings", comps = 1:3) # The three first loadings
plot(pcr_model, "coef", ncomp = 5) # Coefficients
plot(pcr_model, "val") # RMSEP curves
plot(pcr_model, "val", val.type = "MSEP", estimate = "CV") # CV MSEP

library(ggpubr)
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


# 70% of the sample size
sample_size <- floor(0.70 * nrow(bostons))

# set the seed to make your partition reproductible
set.seed(777)
train_ind <- sample(seq_len(nrow(bostons)), size = sample_size)
# splitting data by the vector
trainb <- bostons[train_ind, ]
testb <- bostons[-train_ind, ]

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
# R2 = 78.1 %

# Testing how the predicted values fit with the actual values
pred <- predict(m2, newdata = trainb)

plot(pred,trainb$nox, xlab="predicted",ylab="actual")
abline(a=0,b=1)

#Calculating mean squared error
pred.test <- predict(m2, newdata = testb)
mse.train <- summary(m2)$sigma^2
mse.test  <- sum((pred.test - testb$nox)^2)/(nrow(testb)-length(trainb)-2)

mse.train
mse.test
mse.train/mse.test
# Mean squared error ratio is 92.3 %
# This suggests that the new data can be predicted with good accuracy, scaling improved 0.3 %.
