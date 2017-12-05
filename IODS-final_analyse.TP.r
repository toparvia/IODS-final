# title: "IODS Final Assignment: Boston Dataset"
# author: "Tuure Parviainen"
# date: "December 5, 2017"

# As there is not that much manipulation to be done, for such used test dataset as Boston, I'll look into normality of assumptions in the data as often at least prices are log-normal, so perhaps some transformations are in order. I'll use shapiro test and fitdistrplus package to analyse the distribution of the each variable.

# Loading packages

require("easypackages") # for loading and installing many packages at one time
packages_to_load <- c("dplyr", "tidyverse", "ggplot2", "GGally", "MASS", "reshape2", "FactoMineR", "factoextra")
packages(packages_to_load, prompt = TRUE) # lock'n'load install/load

# rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # to set


# 1) Normal unchanged dataset
boston <- read.table("data/boston.txt")
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

m <- lm(nox ~ ., data = boston)

# Model is able to explain 78 %
stepAIC(m, direction = "both")

m2 <- lm(nox ~ crim + indus + age + dis + rad + ptratio + medv, data = boston)

# Which is better m or m2 by AIC?

AIC(m,m2)
anova(m2, m)

# By AIC more parsimonious model is better.
plot(m2)

# Let's look at the data again with the same variables as in M2
g = ggpairs(boston,columns = c(1,3,7:8,10,14), lower = list(continuous = my_fn))
g
summary(m2)
# We still have quite many dimensions in the data. Maybe we could lower the dimensions using PCA?

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

With LM we had: crim + indus + age + dis + rad + ptratio + medv
With PCA:  medv + indus + dis + lstat and nox

#http://www.win-vector.com/blog/2016/05/pcr_part1_xonly/
