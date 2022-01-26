rm(list=ls())
library(GeneNet)

myfinaldata=readRDS("myfinaldata_dtw.rds")
dim(myfinaldata) #480*277

myfinalid=colnames(myfinaldata)

mydata_long <- as.longitudinal(myfinaldata, repeats=rep(1,480), 
                               time=1:480)
#plot(mydata_long,1:9)

# load functions for estimating shrinkage VAR model
source("mvr.shrink.R") # define mvr.shrink function
source("shrink.VAR.R") 

# functions for plotting the resulting VAR network
source("makedotfile.VAR.R")

################# estimate VAR coefficients ##################### 

# input: n x p data matrix, has to be of type "longitudinal",
# output:  
#  B_est: estimated VAR regression coefficients
#  r.mat: corresponding partial correlations for B_est
#  l.est: shrinkage parameter for correlations
#  l.var.est: shrinkage parameter for variance vector

myVAR.coeff <- shrink.VAR(mydata_long,lambda=0)

# assign local fdr values to each edge (prob=1-local fdr)
myresults.VAR <- test.VAR(myVAR.coeff) # 92416*6



# determine significant edges (local fdr < 0.2)
myresults.sig <- myresults.VAR[myresults.VAR$prob > 0.50,] # 115*6
dim(myresults.sig)[1] # number of significant edges
#[1] 94

myvar_aanat1=myresults.VAR[myresults.VAR$node1_past==7,] #304*6,start from 1588th row
myvar_aanat2=myresults.VAR[myresults.VAR$node2_future==7,] #304*6,start from 92th row

makedotfile_VAR(myvar_aanat2, "finalaanat.dot") 
