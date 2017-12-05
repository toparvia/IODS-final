# title: "Chapter 4. Clustering and classification"
# author: "Tuure Parviainen"
# date: "November 29, 2017"

# meta http://hdr.undp.org/en/content/human-development-index-hdi
# Info http://hdr.undp.org/sites/default/files/hdr2015_technical_notes.pdf

require("easypackages") # for loading and installing many packages at one time
packages_to_load <- c("dplyr", "tidyverse", "ggplot2", "GGally")
packages(packages_to_load, prompt = TRUE) # lock'n'load install/load

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # to set

hd <- read.csv("http://s3.amazonaws.com/assets.datacamp.com/production/course_2218/datasets/human_development.csv", stringsAsFactors = F)
glimpse(hd)
gii <- read.csv("http://s3.amazonaws.com/assets.datacamp.com/production/course_2218/datasets/gender_inequality.csv", stringsAsFactors = F, na.strings = "..")
glimpse(gii)

hdnames <- c("HDI.rank", "country", "HDI", "LE","edu", "edu.mean", "GNI","GNI.HDI.rank")
colnames(hd) <- hdnames

giinames <- c("GII.rank", "country", "GII", "MMratio", "ABrate", "rep.par", "sec.edu.F", "sec.edu.M", "LFP.F", "LFP.M")
colnames(gii)<- giinames

gii <- gii %>%  mutate(edu.FM.ratio = sec.edu.F / sec.edu.M)
gii <- gii %>%  mutate(LFP.FM.ratio = LFP.F / LFP.M)

human <- inner_join(hd, gii, by = "country")

write.csv(human, file = "data/human.csv")

summary(human)
dim(human)

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) +
    geom_point() +
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}
g = ggpairs(human,columns = c(1,3:6,8:19), lower = list(continuous = my_fn))
g

