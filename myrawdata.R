##### Starting from the raw data
rm(list=ls())
library(GeneNet)

G_total=read.csv("FPKMnew.csv")

colnames(G_total)

G_total=as.data.frame(G_total)

ID=G_total[,1]
Symbol=G_total[,483]

GG=G_total[,2:481] #18670*480
rownames(GG)=ID
GG=G_total[,2:481] #18670*480
rownames(GG)=ID

rowsd=apply(GG,1,sd)
rowmean=rowMeans(GG)  ### include 0

newrowsd=rowsd/rowmean ### include 791 NAs
na_location=which(is.na(newrowsd))

plot(rowsd,type="l",xlab="Genes",ylab="Sd")
plot(rowsd,xlab="Genes",ylab="Sd",cex=0.5,col=4)
abline(h=0,col="red")

orderindex=order(rowsd)
ordersd=rowsd[orderindex]
plot(ordersd,xlab="Ordered genes",ylab="Sd")

rowsd0=rowsd[rowsd==0] #there are 791 genes whose sd is 0

mynewsd=newrowsd[!c(1:18670) %in% na_location] ##17879 genes left

## remove genes whose sd is 0, 29486 genes left
myrowsd=rowsd[!rowsd==0]

plot(myrowsd,type="l",xlab="Genes",ylab="Sd")
abline(h=0,col="blue")

plot(mynewsd,xlab="Genes",ylab="Sd",cex=0.5,col=4)
abline(h=0,col="red")

myorderindex=order(mynewsd)
myordersd=mynewsd[myorderindex]
plot(myordersd,xlab="Ordered genes",ylab="Sd",cex=0.5,col=4)
abline(h=0,col="red")

################### remove genes with sd of 0 ####################
mynames=names(myrowsd)

mydata=GG[mynames,]  ## 17879*480

saveRDS(mydata,"mydata.rds")


#####################################################################
## Fuzzy C means clustering
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Mfuzz")
library(Mfuzz)

mydatacl=GG[mynames,]  #17879*480, mydata

mydf <- new("ExpressionSet", exprs = as.matrix(mydatacl))
mydf <- standardise(mydf)
m <- mestimate(mydf)

#### 9 clusters #################################################
mycl <- mfuzz(mydf, c = 9, m = m)
mfuzz.plot(
  mydf,
  mycl,
  mfrow = c(3, 3),
  new.window = FALSE)

write.csv(mycl$cluster, file = "mynewmfuzz9.csv")

mycl9_6=names(mycl$cluster[mycl$cluster==6]) #1093 genes
mycldata9_6=GG[mycl9_6,]
saveRDS(mycldata9_6,"mycldata9_6.rds")

mycl9_6_mynames=G_total[G_total$Ensembl_ID %in% mycl9_6,483]
write.table(mycl9_6_mynames,"mycl9_6_mynames.txt",row.names = F,col.names = F,quote=FALSE)
