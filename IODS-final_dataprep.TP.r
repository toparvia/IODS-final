# title: "IODS Final Assignment: Boston Dataset"
# author: "Tuure Parviainen"
# date: "December 5, 2017"

# As there is not that much manipulation to be done, for such used test dataset as Boston, I'll look into normality of assumptions in the data as often at least prices are log-normal, so perhaps some transformations are in order. I'll use shapiro test and fitdistrplus package to analyse the distribution of the each variable.

# Loading packages

require("easypackages") # for loading and installing many packages at one time
packages_to_load <- c("dplyr", "tidyverse", "ggplot2", "fitdistrplus", "MASS", "reshape2")
packages(packages_to_load, prompt = TRUE) # lock'n'load install/load

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # to set

data(Boston)

# Histogram of all variables
ggplot(data = melt(Boston), mapping = aes(x = value)) +
  geom_histogram(bins = 10) + facet_wrap(~variable, scales = 'free_x')

# Are these really normal?

shapiro.test(Boston[,1])

sha.test <- do.call(rbind, lapply(Boston, function(x) shapiro.test(x)[c("statistic", "p.value")]))

# Wiki: The null-hypothesis of this test is that the population is normally distributed. Thus, if the p-value is less than the chosen alpha level, then the null hypothesis is rejected and there is evidence that the data tested are not from a normally distributed population; in other words, the data are not normal. On the contrary, if the p-value is greater than the chosen alpha level, then the null hypothesis that the data came from a normally distributed population cannot be rejected (e.g., for an alpha level of 0.05, a data set with a p-value of 0.02 rejects the null hypothesis that the data are from a normally distributed population).

sha.test[,2] >= 0.05

# All false. We have non-normal data. However, Shapiro test does not do well when sample size is large. So we make QQ-plots.

# x11()
# QQ-plots
ggplot(data = melt(Boston), mapping = aes(sample = value)) +
  stat_qq() + facet_wrap(~variable, scales = 'free')

# I sense trouble with normality assumption required for regression.

# We will look at the different distributions fit to look for possible solution for skewdness.
fitW <- fitdist(Boston[,1], "weibull")
fitg <- fitdist(Boston[,1], "gamma")
fitln <- fitdist(Boston[,1], "lnorm")
fitn <- fitdist(Boston[,1], "norm")

cdfcomp(list(fitW, fitg, fitln, fitn), legendtext=c("Weibull", "gamma", "lognormal", "normal"))
denscomp(list(fitW, fitg, fitln, fitn), legendtext=c("Weibull", "gamma", "lognormal", "normal"))
qqcomp(list(fitW, fitg, fitln, fitn), legendtext=c("Weibull", "gamma", "lognormal", "normal"))
ppcomp(list(fitW, fitg, fitln, fitn), legendtext=c("Weibull", "gamma", "lognormal", "normal"))
gofstat(list(fitW, fitg, fitln, fitn), fitnames=c("Weibull", "gamma", "lognormal", "normal"))

descdist(Boston[,1], boot = 100)

# Normal distribution was the least likely for crim by AIC.

# We will compare two different transformations for data. Log and square root transformation.
BostonSq <- sqrt(Boston) %>% scale() %>% as.data.frame() %>%  melt()
BostonSq$group <- "Squared"
BostonL <- log(Boston) %>% scale() %>% as.data.frame() %>%  melt()
BostonL$group <- "Log"
BostonN <- Boston %>% scale() %>% as.data.frame() %>%  melt()
BostonN$group <- "Normal"

BostonGraph <- rbind(BostonSq, BostonL, BostonN)
BostonGraph <-BostonGraph[complete.cases(BostonGraph),]

# Here we plot all the distributions to a QQ plots
ggplot(data = BostonGraph, mapping = aes(sample = value)) +
  stat_qq(aes(color=group)) + facet_wrap(~variable, scales = 'free')

# Here we see the impact of the scaled and transformed crim variable
ggplot(data = subset(BostonGraph,variable == "crim"), mapping = aes(x = value)) +
  geom_histogram(bins = 15, aes(color=group)) + facet_wrap(~group)

# We also look at the median value (medv)
ggplot(data = subset(BostonGraph,variable == "medv"), mapping = aes(x = value)) +
  geom_histogram(bins = 15, aes(color=group)) + facet_wrap(~group)

# No need for transformation (medv)

# Given that the data has consists already of indexes and ratios any further transformation is not probably necessary.
# However, for the assignment.

# We will have three datasets:
# Empty enviroment.
rm(list = ls())
# 1) Normal unchanged dataset
data(Boston)
head(Boston)
BostonN <- Boston # For some reason didn't manage to read directly
# 2) Scaled dataset
BostonS <- scale(BostonN) %>% as.data.frame()
# 3) Scaled dataset with the crim transformed with the square transformation (log can't deal with zeros).
BostonSSq <- Boston
BostonSSq$crim <-  sqrt(Boston$crim)
BostonSSq <- BostonSSq %>% scale() %>% as.data.frame()
# General note: It seems that some of the data has been censored for the high tax category

# Writing the files to Data folder
write.table(Boston, file = "data/boston.txt")
write.table(BostonS, file = "data/bostons.txt")
write.table(BostonSSq, file = "data/bostonssq.txt")
