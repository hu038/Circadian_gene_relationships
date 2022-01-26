##### Starting from the raw data
rm(list=ls())
library(GeneNet)

G_total=read.csv("FPKMnew.csv")

colnames(G_total)

G_total=as.data.frame(G_total)

ID=G_total[,1]
Symbol=G_total[,483]