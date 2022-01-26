###################### Data: myfinaldata 304genes ###############
############ Methods: GeneNet ########################### 

rm(list=ls())

library(GeneNet)

myfinaldata=readRDS("myfinaldata_dtw.rds")
dim(myfinaldata) #480*304

myfinalid=colnames(myfinaldata)

mydata_long <- as.longitudinal(myfinaldata, repeats=rep(1,480), 
                               time=1:480)

plot(mydata_long,1:9)

G_total[G_total[,1]=="ENSRNOG00000009177",][483]
G_total[G_total[,1]=="ENSRNOG00000004589",][483]
G_total[G_total[,1]=="ENSRNOG00000042978",][483]
G_total[G_total[,1]=="ENSRNOG00000006068",][483]
G_total[G_total[,1]=="ENSRNOG00000018051",][483]
G_total[G_total[,1]=="ENSRNOG00000006470",][483]
G_total[G_total[,1]=="ENSRNOG00000011182",][483]
G_total[G_total[,1]=="ENSRNOG00000017927",][483]
G_total[G_total[,1]=="ENSRNOG00000017333",][483]

par(mfrow=c(3,3))
plot(mydata_long[,1],type="l",cex=0.5,xlab="time",ylab="value",
     main="Fcer1a",col=4)
plot(mydata_long[,2],type="l",cex=0.5,xlab="time",ylab="value",
     main="Galnt16",col=4)
plot(mydata_long[,3],type="l",cex=0.5,xlab="time",ylab="value",
     main="Ncald",col=4)
plot(mydata_long[,4],type="l",cex=0.5,xlab="time",ylab="value",
     main="Tmem117",col=4)
plot(mydata_long[,5],type="l",cex=0.5,xlab="time",ylab="value",
     main="Farp2",col=4)
plot(mydata_long[,6],type="l",cex=0.5,xlab="time",ylab="value",
     main="Camk1g",col=4)
plot(mydata_long[,7],type="l",cex=0.5,xlab="time",ylab="value",
     main="Aanat",col=4)
plot(mydata_long[,8],type="l",cex=0.5,xlab="time",ylab="value",
     main="Drd4",col=4)
plot(mydata_long[,9],type="l",cex=0.5,xlab="time",ylab="value",
     main="Syt4",col=4)
dev.off()
#' # Compute "shrinkage" Partial Dynamic Correlations 
## 304*304 matrix
pcor.dyn = ggm.estimate.pcor(mydata_long, method = "dynamic")
#Specified shrinkage intensity lambda (correlation matrix): 0.0152 
# # G=304, E=304*303/2=46056
# indexes = sm.index(pcor.dyn)
# colnames(indexes) = c("node1", "node2")
# pcor = sm2vec(pcor.dyn) 
# w = cbind(pcor, indexes)
# fdr.cor = fdrtool(w[, 1], statistic = "correlation",plot=TRUE)
# fdr.norm = fdrtool(w[, 1], statistic = "norm",plot=TRUE)
# 
# GG.undir.edges = network.test.edges(pcor.dyn,direct=FALSE,fdr=TRUE)




###########################################################################
arth.edges = network.test.edges(pcor.dyn,direct=TRUE,fdr=TRUE)
dim(arth.edges)
#[1]46056    10
# pval: p-values
# qval: tail area-based FDR (=Fdr)
# prob: probility that edge is nonzero (1-local fdr)
#log.spvar: log ratio of standardized partial variance (determines direction)

# which(myfinalid=="ENSRNOG00000011182")
# [1] 7


mygenenet_aanat1=arth.edges[arth.edges$node1==7,] ##297*10, start from 97th row
mygenenet_aanat2=arth.edges[arth.edges$node2==7,] ## 6*10, start from 9688th row
mygenenet_aanat2

myfinalarth.net[myfinalarth.net$node1==7,]
myfinalarth.net[myfinalarth.net$node2==7,]

###################################################################################
#' We use the strongest 150 edges:
myfinalarth.net = extract.network(arth.edges, method.ggm="number", cutoff.ggm=150)
#Significant edges:  150 
#Corresponding to  0.33 %  of possible edges 

#Significant directions:  1362 
#Corresponding to  2.96 %  of possible directions 
#Significant directions in the network:  2 
#Corresponding to  1.33 %  of possible directions in the network

myfinalarth.net[myfinalarth.net$node1==7,]
myfinalarth.net[myfinalarth.net$node2==7,]

#' # Construct Graph

library("graph")

node.labels = as.character(1:ncol(mydata_long))
gr = network.make.graph(myfinalarth.net, node.labels, drop.singles=TRUE) 
gr
# A graphNEL graph with directed edges
# Number of Nodes = 180 
# Number of Edges = 298
#' Some information about the graph

#' Number of nodes:
num.nodes(gr)
#180

#' Correlations:
str(edge.info(gr)$weight) 

#' Number of directed ("forward") and undirected ("none") edges:
table(edge.info(gr)$dir)
#forward    none 
#2     148 
#'
#' # Well-Connected Nodes

#' Nodes connected with many edges:
sort(node.degree(gr), decreasing=TRUE)[1:10]
#201  55 104 123 142 189 196  41  57  63 
#7   4   4   4   4   4   4   3   3   3 

# myfinalid[267]
# # /* [1] "ENSRNOG00000006950"
# 
# myfinalid[183]
# # /* [1] "ENSRNOG00000018693"
# 
# myfinalid[219]
# # /* [1]"ENSRNOG00000012178"
# 
# myfinalid[242]
# # /* [1] "ENSRNOG00000000902"
# 
# myfinalid[265]
# # /* [1] "ENSRNOG00000050647"
#'
#' # Plot Network

library("Rgraphviz")

#' For a more beautiful plot of the network set node and edge parameters:

#' Set global node and edge attributes:
globalAttrs = list()
globalAttrs$edge = list(color = "black", lty = "solid", lwd = 1, arrowsize=1)
globalAttrs$node = list(fillcolor = gray(.95), shape = "ellipse", fixedsize = FALSE)

#' Set attributes of some particular nodes:
nodeAttrs = list()
nodeAttrs$fillcolor = c('201' = "yellow", "55" = "yellow","7"="red"
) # highlight hub nodes

#' Set edge attributes:
edi = edge.info(gr) # edge directions and correlations
edgeAttrs = list()
edgeAttrs$dir =  edi$dir # set edge directions 
cutoff = quantile(abs(edi$weight), c(0.2, 0.8)) # thresholds for line width / coloring
edgeAttrs$lty = ifelse(edi$weight < 0, "dotted", "solid") # negative correlation
edgeAttrs$color = ifelse( abs(edi$weight <= cutoff[1]), "grey", "black") # lower 20% quantile
edgeAttrs$lwd = ifelse(abs(edi$weight >= cutoff[2]), 2, 1) # upper 20% quantile

#+ fig.width=8, fig.height=7
plot(gr, attrs = globalAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs, "fdp")

myfinalid[7]
#[1] "ENSRNOG00000011182"
myfinalid[201]
#[1] "ENSRNOG00000056153"
G_total[G_total[,1]=="ENSRNOG00000056153",][483]
#Symbol
#16189 Tent5b
myfinalid[43]
#[1] "ENSRNOG00000039856"
G_total[G_total[,1]=="ENSRNOG00000039856",][483]
#Symbol
#7922 Lrrc73
myfinalid[224]
#[1] "ENSRNOG00000006950"
G_total[G_total[,1]=="ENSRNOG00000006950",][483]
#Symbol
#11246  Padi3
myfinalid[91]
#[1] "ENSRNOG00000019737"
G_total[G_total[,1]=="ENSRNOG00000019737",][483]
#Symbol
#14358 Sema4a
myfinalid[55]
#[1] "ENSRNOG00000024717"
G_total[G_total[,1]=="ENSRNOG00000024717",][483]
#Symbol
#4568  Espnl
