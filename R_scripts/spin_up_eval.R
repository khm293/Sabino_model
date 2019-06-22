#script to make WT plot of Sabino domain for one pressure snapshot
#june 21, 2019 

library(ggplot2)

source('PFB-ReadFcn.R')

#read pressure pfb file

setwd('~/Sabino_model/')
pressure=readpfb('sabino.out.press.00009.pfb', verbose = FALSE)

#grab bottom layer pressure

btm_pressure=pressure[,,1]-100
dtw=500-btm_pressure[,]
image(dtw)
contour(dtw)
