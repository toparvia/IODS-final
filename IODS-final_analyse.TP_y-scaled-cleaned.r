# title: "IODS Final Assignment: Boston Dataset y-aware scaling"
# author: "Tuure Parviainen"
# date: "December 5, 2017"

require("easypackages") # for loading and installing many packages at one time
packages_to_load <- c("dplyr", "ggplot2", "reshape2","visreg", "tidyverse")
packages(packages_to_load, prompt = TRUE) # lock'n'load install/load


#DataSet Loading and division to train and test sets
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # to set

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

# apply projection for TRAIN
projectedTrain <- as.data.frame(dmTrain %*% proj,
                                stringsAsFactors = FALSE)
# plot data sorted by principal components
projectedTrain$nox <- trainbNTreatedYScaled$nox
ScatterHistN(projectedTrain,'PC1','PC2','nox',
             "Y-Scaled Training Data projected to first two principal components")

model <- lm(nox~PC1+PC2+PC3,data=projectedTrain)
model2 <- lm(nox~PC1+PC2,data=projectedTrain)
AIC(model,model2)
summary(model2)


dev.off()
plot(model, 2)


projectedTrain$estimate <- predict(model,newdata=projectedTrain)
trainrsq = rsq(projectedTrain$estimate,projectedTrain$nox)

ScatterHist(projectedTrain,'estimate','nox','Recovered model versus truth (y aware PCA train)',
            smoothmethod='identity',annot_size=3)
# Upper end of values consistently predicts lower than actual

# apply projection for TEST
projectedTest <- as.data.frame(dmTest %*% proj,
                               stringsAsFactors = FALSE)
# plot data sorted by principal components
projectedTest$nox <- trainbNTreatedYScaled$nox
ScatterHistN(projectedTest,'PC1','PC2','nox',
             "Y-Scaled Test Data projected to first two principal components")

projectedTest$estimate <- predict(model2,newdata=projectedTest)
testrsq = rsq(projectedTest$estimate,projectedTest$nox)
testrsq

# We can explain 73.1 % of the variation on the test dataset.

ScatterHist(projectedTest,'estimate','nox','Recovered model versus truth (y aware PCA test)',
            smoothmethod='identity',annot_size=3)

######################################################################################
# GAM - Generalized additive model on PCA

# Since we saw that there was some data that we coudn't explain with linear models,
# we will try to use more flexible approaches, namely Generilized additive models on the
# y-scaled PCA dimensions. Basically we hope that we are able to use only few dimensions to
# model the data and capture the variability not accounted in the linear model with GAM.

require(mgcv)

gmodel <- gam(nox ~ s(PC1) + PC2 + PC3, family  = gaussian, data=projectedTrain,
              trace=FALSE)

gmodel2 <- gam(nox ~ s(PC1) + PC2, family  = gaussian, data=projectedTrain,
               trace=FALSE)
gmodel3 <- gam(nox ~ s(PC1, bs="cr") + PC2, data=projectedTrain)
gmodel3.5 <- gam(nox ~ s(PC1, bs="cr")+ s(PC2, bs="cr"), data=projectedTrain)
gmodel4 <- gam(nox ~ s(PC1, bs="cr") + PC2 + PC3, data=projectedTrain)
gmodel5 <- gam(nox ~ te(PC1, PC2), data=projectedTrain)
gmodel6 <- gam(nox ~ te(PC1, PC2) + PC3, data=projectedTrain)


# te =  tensor product smooth, and by smoothing the marginal smooths of PC1 and PC2

AIC(gmodel,gmodel2,gmodel3,gmodel3.5, gmodel4,gmodel5,gmodel6)

#models 5 and 6 give most negative (smallest values)

# AIC = -2Ln(L)+ 2k

#par(mfrow=c(1,3)) #to partition the Plotting Window
plot(gmodel5,se = TRUE)
plot(gmodel6,se = TRUE)
#se stands for standard error Bands

summary(gmodel5)$r.sq
summary(gmodel5)$sp.criterion
summary(gmodel6)$r.sq
summary(gmodel6)$sp.criterion

anova(gmodel,gmodel2,gmodel3,gmodel4,gmodel5,gmodel6, test="Chisq")

#Model 5 is selected by AIC and Chi squared test

#vis.gam(gmodel5, type='response', plot.type='contour')
vis.gam(gmodel5, main = "Tensor product smooth of 1:2 PC", plot.type = "contour",
        color = "terrain", contour.col = "black", lwd = 2)
visreg2d(gmodel5, xvar='PC1', yvar='nox', scale='response')
visreg2d(gmodel5, xvar='PC2', yvar='nox', scale='response')
#visreg2d(gmodel5, xvar='PC3', yvar='nox', scale='response')
vis.gam(gmodel5, n.grid = 50, theta = 125, phi = 32, zlab = "",
        ticktype = "detailed", color = "topo", main = "Tensor product smooth of 1:2 PC")


# Lets confirm out model with replications of K (the smoother term)
par(mfrow=c(1,3))
gam.check(gmodel5, k.rep=1000, type=c("pearson"))
qq.gam(gmodel6)
plot(gmodel5)
dev.off()

#visreg2d(gmodel6, xvar='PC3', yvar='nox', scale='response')

#We have perfectly divided are variance to two dimensions, but we can only make use of this with similar data.
# How about predictions on our test data?

projectedTest$g.estimate <- predict(gmodel5, newdata = projectedTest)

trainrsq_gam = rsq(projectedTest$g.estimate ,projectedTest$nox)
trainrsq_gam
testrsq

#GAM model provides better accuracy in predicting test outcome.

ScatterHist(projectedTest,'g.estimate','nox','Recovered model versus truth (y aware PCA test) GAM ',
            smoothmethod='identity',annot_size=3)


