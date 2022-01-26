#######################################################
library(dtw)
library(TSdist)

mycldata9_6=readRDS("mycldata9_6.rds")  ###1093 genes


##### calculate the DTW distances for all (17879) genes ###################
mydata=readRDS("mydata.rds")

myGG=scale(t(mydata))
tempy=as.vector(scale(myployaanat))

n=nrow(mydata)
myeuc=mydtw=myacf=mypacf=mycor=myccor=mysts=myfourier=myper=mytam=rep(0,n)

for(i in 1:n){
  myeuc[i]=TSDistances(myGG[,i],tempy,"euclidean")
  mydtw[i]=TSDistances(myGG[,i],tempy,"dtw")
  myacf[i]=TSDistances(myGG[,i],tempy,"acf")
  mypacf[i]=TSDistances(myGG[,i],tempy,"pacf")
  mycor[i]=TSDistances(myGG[,i],tempy,"cor")
  myccor[i]=TSDistances(myGG[,i],tempy,"ccor")
  mysts[i]=TSDistances(myGG[,i],tempy,"sts")
  myfourier[i]=TSDistances(myGG[,i],tempy,"fourier")
  myper[i]=TSDistances(myGG[,i],tempy,"per")
  mytam[i]=TSDistances(myGG[,i],tempy,"tam")
}

mynames=rownames(mydata)
names(myeuc)=names(mydtw)=names(myacf)=names(mypacf)=names(mycor)=names(myccor)=names(mysts)=names(myfourier)=names(myper)=names(mytam)=mynames

## "NOT RUN! Cost too much time ##############
myman=mymin=myinfnorm=mydissim=myfrechet=mycort=myncd=myintper=mycdm=myspec=rep(0,n)
for(i in 1:n){
  myman[i]=TSDistances(myGG[,i],tempy,"manhattan")
  mymin[i]=TSDistances(myGG[,i],tempy,"minkowski",p=3)
  myinfnorm[i]=TSDistances(myGG[,i],tempy,"infnorm")
  mydissim[i]=TSDistances(myGG[,i],tempy,"dissim")
  myfrechet[i]=TSDistances(myGG[,i],tempy,"frechet")
  mycort[i]=TSDistances(myGG[,i],tempy,"cort")
  myncd[i]=TSDistances(myGG[,i],tempy,"ncd")
  myintper[i]=TSDistances(myGG[,i],tempy,"int.per")
  mycdm[i]=TSDistances(myGG[,i],tempy,"cdm")
  myspec[i]=TSDistances(myGG[,i],tempy,"spec.glk")
}
names(myman)=names(mymin)=names(myinfnorm)=names(mydissim)=names(myfrechet)=names(mycort)=names(myncd)=names(myintper)=names(mycdm)=names(myspec)=mynames

## The rank of total 17879 genes
myclnames=rownames(mycldata9_6) ###names of 1093 genes
difclnames=setdiff(mynames,myclnames) ###names of 16786 genes

################# euc ##########################
eucorder=order(myeuc)
myordereuc=myeuc[eucorder]
which(names(myordereuc)=="ENSRNOG00000011182")
# [1] 2

mycleuc=myeuc[myclnames]
cleucorder=order(mycleuc)
myclordereuc=mycleuc[cleucorder]
which(names(myclordereuc)=="ENSRNOG00000011182")
# [1] 2

cleuc=match(names(myclordereuc),names(myordereuc)) ## the location of 1093 genes
names(cleuc)=names(myclordereuc)
cleuc["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 2 
range(cleuc)
#[1]     1 3907
length(which(cleuc<1093))
#[1] 825

m=length(myclnames)
ranky=c(1:m)
mypvalueeuc=rep(0,(m-70+1))
for(j in 70:m){
  mypvalueeuc[j-69]=wilcox.test(cleuc[1:j],ranky[1:j])$p.value
}

which.max(mypvalueeuc<0.05)
#[1] 670
mypvalueeuc[670]
670+69
#[1] 739

eucnames1=names(cleuc[1:739])
myeuc1=myordereuc[eucnames1]
range(myeuc1)
#[1] 10.63821 22.89636

mydifeuc=myeuc[difclnames]  
difeucorder=order(mydifeuc)
mydifordereuc=mydifeuc[difeucorder]
difcleuc=match(names(mydifordereuc),names(myordereuc))  ##the location of 16786 genes  
range(difcleuc)
#[1]   238 17879
length(which(mydifordereuc<max(myeuc1)))
# [1] 168
eucnames2=names(which(mydifordereuc<max(myeuc1)))

myeuc2=myordereuc[eucnames2]

#739+108 genes
myfinaldata_euc=myGG[,c(eucnames1,eucnames2)] ## 907genes, 907*1135


################# DTW ##########################
dtworder=order(mydtw)
myorderdtw=mydtw[dtworder]
which(names(myorderdtw)=="ENSRNOG00000011182")
# [1] 7

myclnames=rownames(mycldata9_6) ###names of  genes
mycldtw=mydtw[myclnames]
cldtworder=order(mycldtw)
myclorderdtw=mycldtw[cldtworder]
which(names(myclorderdtw)=="ENSRNOG00000011182")
# [1] 7

cldtw=match(names(myclorderdtw),names(myorderdtw)) ## the location of 1093 genes
names(cldtw)=names(myclorderdtw)
cldtw["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 7 
range(cldtw)
#[1]     1 13225
length(which(cldtw<1093))
#[1] 626

m=length(myclnames)
ranky=c(1:m)
mypvaluedtw=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluedtw[j-69]=wilcox.test(cldtw[1:j],ranky[1:j])$p.value
}

which.max(mypvaluedtw<0.05)
#[1] 184
mypvaluedtw[184]
184+69
#[1] 253

dtwnames1=names(cldtw[1:253])
mydtw1=myorderdtw[dtwnames1]
range(mydtw1)
#[1] 133.8994 233.8590

mydifdtw=mydtw[difclnames]
difdtworder=order(mydifdtw)
mydiforderdtw=mydifdtw[difdtworder]
difcldtw=match(names(mydiforderdtw),names(myorderdtw))  ##the location of 16786 genes  
range(difcldtw)
#[1]    63 17879
length(which(mydiforderdtw<max(mydtw1)))
# [1] 51
dtwnames2=names(which(mydiforderdtw<max(mydtw1)))

mydtw2=myorderdtw[dtwnames2]

#253+51 genes
myfinaldata_dtw=myGG[,c(dtwnames1,dtwnames2)] ## 304genes, 480*304

max(clman,cldissim,clcort,cleuc,clcor,
    clfourier,clccor,cldtw,clmin)
plot(clmin,cex=0.05,col=3,lty=3,xlab="Genes of the selected cluster",ylab="Rank",main="Different distances")
lines(clman,cex=0.05,col=4,lty=4)
lines(cldissim,cex=0.05,col=5,lty=5)
lines(clcort,cex=0.05,col=6,lty=6)
lines(cleuc,cex=0.05,col=7,lty=7)
lines(cldtw,cex=0.05,col=1,lty=1)
lines(ranky,col=2,cex=5,lty=2)
legend("topleft", legend = c("DTW","y=x","MIN","MAN","DISSIM",
                             "CORT","EUC"), 
       col = 1:7,lty=1:7)

saveRDS(myfinaldata_euc,"myfinaldata_euc.rds")
saveRDS(myfinaldata_dtw,"myfinaldata_dtw.rds")


################# acf ##########################
acforder=order(myacf)
myorderacf=myacf[acforder]
which(names(myorderacf)=="ENSRNOG00000011182")
# [1] 30

myclnames=rownames(mycldata9_6) ###names of  genes
myclacf=myacf[myclnames]
clacforder=order(myclacf)
myclorderacf=myclacf[clacforder]
which(names(myclorderacf)=="ENSRNOG00000011182")
# [1] 21

clacf=match(names(myclorderacf),names(myorderacf)) ## the location of 1093 genes
names(clacf)=names(myclorderacf)
clacf["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 30 
range(clacf)
#[1]     1 10955
length(which(clacf<1093))
#[1] 515

m=length(myclnames)
ranky=c(1:m)
mypvalueacf=rep(0,(m-70+1))
for(j in 70:m){
  mypvalueacf[j-69]=wilcox.test(clacf[1:j],ranky[1:j])$p.value
}

which.max(mypvalueacf<0.05)
#[1] 1
mypvalueacf[1]
wilcox.test(clacf[1:35],ranky[1:35])


acfnames1=names(clacf[1:35])
myacf1=myorderacf[acfnames1]
range(myacf1)
#[1] 0.6866002 1.5463831

mydifacf=myacf[difclnames]
difacforder=order(mydifacf)
mydiforderacf=mydifacf[difacforder]
difclacf=match(names(mydiforderacf),names(myorderacf))  ##the location of 16786 genes  
range(difclacf)
#[1]    2 17879
length(which(mydiforderacf<max(myacf1)))
# [1] 15
acfnames2=names(which(mydiforderacf<max(myacf1)))

myacf2=myorderacf[acfnames2]

#35+15 genes
myfinaldata_acf=myGG[,c(acfnames1,acfnames2)] ## 50genes, 480*304

################# pacf ##########################
pacforder=order(mypacf)
myorderpacf=mypacf[pacforder]
which(names(myorderpacf)=="ENSRNOG00000011182")
# [1] 17049

myclnames=rownames(mycldata9_6) ###names of  genes
myclpacf=mypacf[myclnames]
clpacforder=order(myclpacf)
myclorderpacf=myclpacf[clpacforder]
which(names(myclorderpacf)=="ENSRNOG00000011182")
# [1] 955

clpacf=match(names(myclorderpacf),names(myorderpacf)) ## the location of 1093 genes
names(clpacf)=names(myclorderpacf)
clpacf["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 17049 
range(clpacf)
#[1]     83 17877
length(which(clpacf<1093))
#[1] 21

m=length(myclnames)
ranky=c(1:m)
mypvaluepacf=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluepacf[j-69]=wilcox.test(clpacf[1:j],ranky[1:j])$p.value
}

which.max(mypvaluepacf<0.05)
#[1] 1
wilcox.test(clpacf[1:3],ranky[1:3])

pacfnames1=names(clpacf[1:3])
mypacf1=myorderpacf[pacfnames1]
range(mypacf1)
#[1] 1.013852 1.020497

mydifpacf=mypacf[difclnames]
difpacforder=order(mydifpacf)
mydiforderpacf=mydifpacf[difpacforder]
difclpacf=match(names(mydiforderpacf),names(myorderpacf))  ##the location of 16786 genes  
range(difclpacf)
#[1]    1 17879
length(which(mydiforderpacf<max(mypacf1)))
# [1] 123
pacfnames2=names(which(mydiforderpacf<max(mypacf1)))

mypacf2=myorderpacf[pacfnames2]

#3+123 genes
myfinaldata_pacf=myGG[,c(pacfnames1,pacfnames2)] ## 126genes

################# cor ##########################
cororder=order(mycor)
myordercor=mycor[cororder]
which(names(myordercor)=="ENSRNOG00000011182")
# [1] 2

myclnames=rownames(mycldata9_6) ###names of  genes
myclcor=mycor[myclnames]
clcororder=order(myclcor)
myclordercor=myclcor[clcororder]
which(names(myclordercor)=="ENSRNOG00000011182")
# [1] 2

clcor=match(names(myclordercor),names(myordercor)) ## the location of 1093 genes
names(clcor)=names(myclordercor)
clcor["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 2 
range(clcor)
#[1]     1 3907
length(which(clcor<1093))
#[1] 825

m=length(myclnames)
ranky=c(1:m)
mypvaluecor=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluecor[j-69]=wilcox.test(clcor[1:j],ranky[1:j])$p.value
}

which.max(mypvaluecor<0.05)
#[1] 670
mypvaluecor[670]
670+69
#[1] 739

cornames1=names(clcor[1:739])
mycor1=myordercor[cornames1]
range(mycor1)
#[1] 0.4860722 1.0461612

mydifcor=mycor[difclnames]
difcororder=order(mydifcor)
mydifordercor=mydifcor[difcororder]
difclcor=match(names(mydifordercor),names(myordercor))  ##the location of 16786 genes  
range(difclcor)
#[1]    238 17879
length(which(mydifordercor<max(mycor1)))
# [1] 168
cornames2=names(which(mydifordercor<max(mycor1)))

mycor2=myordercor[cornames2]

#739+168 genes
myfinaldata_cor=myGG[,c(cornames1,cornames2)] ## 907genes

################# ccor ##########################
ccororder=order(myccor)
myorderccor=myccor[ccororder]
which(names(myorderccor)=="ENSRNOG00000011182")
# [1] 5

myclnames=rownames(mycldata9_6) ###names of  genes
myclccor=myccor[myclnames]
clccororder=order(myclccor)
myclorderccor=myclccor[clccororder]
which(names(myclorderccor)=="ENSRNOG00000011182")
# [1] 5

clccor=match(names(myclorderccor),names(myorderccor)) ## the location of 1093 genes
names(clccor)=names(myclorderccor)
clccor["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 5 
range(clccor)
#[1]     1 9903
length(which(clccor<1093))
#[1] 474

m=length(myclnames)
ranky=c(1:m)
mypvalueccor=rep(0,(m-70+1))
for(j in 70:m){
  mypvalueccor[j-69]=wilcox.test(clccor[1:j],ranky[1:j])$p.value
}

which.max(mypvalueccor<0.05)
#[1] 1
wilcox.test(clccor[1:63],ranky[1:63])


ccornames1=names(clccor[1:63])
myccor1=myorderccor[ccornames1]
range(myccor1)
#[1] 0.05897944 0.09506795

mydifccor=myccor[difclnames]
difccororder=order(mydifccor)
mydiforderccor=mydifccor[difccororder]
difclccor=match(names(mydiforderccor),names(myorderccor))  ##the location of 16786 genes  
range(difclccor)
#[1]    11 17879
length(which(mydiforderccor<max(myccor1)))
# [1] 18
ccornames2=names(which(mydiforderccor<max(myccor1)))

myccor2=myorderccor[ccornames2]

#63+18 genes
myfinaldata_ccor=myGG[,c(ccornames1,ccornames2)] ## 81genes


################# sts ##########################
stsorder=order(mysts)
myordersts=mysts[stsorder]
which(names(myordersts)=="ENSRNOG00000011182")
# [1] 34

myclnames=rownames(mycldata9_6) ###names of  genes
myclsts=mysts[myclnames]
clstsorder=order(myclsts)
myclordersts=myclsts[clstsorder]
which(names(myclordersts)=="ENSRNOG00000011182")
# [1] 24

clsts=match(names(myclordersts),names(myordersts)) ## the location of 1093 genes
names(clsts)=names(myclordersts)
clsts["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 34 
range(clsts)
#[1]     1 14506
length(which(clsts<1093))
#[1] 509

m=length(myclnames)
ranky=c(1:m)
mypvaluests=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluests[j-69]=wilcox.test(clsts[1:j],ranky[1:j])$p.value
}

which.max(mypvaluests<0.05)
#[1] 1
mypvaluests[1]
wilcox.test(clsts[1:34],ranky[1:34])

stsnames1=names(clsts[1:34])
mysts1=myordersts[stsnames1]
range(mysts1)
#[1] 13.30450 16.43646

mydifsts=mysts[difclnames]
difstsorder=order(mydifsts)
mydifordersts=mydifsts[difstsorder]
difclsts=match(names(mydifordersts),names(myordersts))  ##the location of 16786 genes  
range(difclsts)
#[1]    3 17879
length(which(mydifordersts<max(mysts1)))
# [1] 16
stsnames2=names(which(mydifordersts<max(mysts1)))

mysts2=myordersts[stsnames2]

#34+16 genes
myfinaldata_sts=myGG[,c(stsnames1,stsnames2)] ## 50genes


################# fourier ##########################
fourierorder=order(myfourier)
myorderfourier=myfourier[fourierorder]
which(names(myorderfourier)=="ENSRNOG00000011182")
# [1] 2

myclnames=rownames(mycldata9_6) ###names of  genes
myclfourier=myfourier[myclnames]
clfourierorder=order(myclfourier)
myclorderfourier=myclfourier[clfourierorder]
which(names(myclorderfourier)=="ENSRNOG00000011182")
# [1] 2

clfourier=match(names(myclorderfourier),names(myorderfourier)) ## the location of 1093 genes
names(clfourier)=names(myclorderfourier)
clfourier["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 2 
range(clfourier)
#[1]     1 3900
length(which(clfourier<1093))
#[1] 823

m=length(myclnames)
ranky=c(1:m)
mypvaluefourier=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluefourier[j-69]=wilcox.test(clfourier[1:j],ranky[1:j])$p.value
}

which.max(mypvaluefourier<0.05)
#[1] 672
mypvaluefourier[672]
672+69
#[1] 741

fouriernames1=names(clfourier[1:741])
myfourier1=myorderfourier[fouriernames1]
range(myfourier1)
#[1] 164.8395 355.3202

mydiffourier=myfourier[difclnames]
diffourierorder=order(mydiffourier)
mydiforderfourier=mydiffourier[diffourierorder]
difclfourier=match(names(mydiforderfourier),names(myorderfourier))  ##the location of 16786 genes  
range(difclfourier)
#[1]    238 17879
length(which(mydiforderfourier<max(myfourier1)))
# [1] 171
fouriernames2=names(which(mydiforderfourier<max(myfourier1)))

myfourier2=myorderfourier[fouriernames2]

#741+171 genes
myfinaldata_fourier=myGG[,c(fouriernames1,fouriernames2)] ## 912genes

################# per ##########################
perorder=order(myper)
myorderper=myper[perorder]
which(names(myorderper)=="ENSRNOG00000011182")
# [1] 64

myclnames=rownames(mycldata9_6) ###names of  genes
myclper=myper[myclnames]
clperorder=order(myclper)
myclorderper=myclper[clperorder]
which(names(myclorderper)=="ENSRNOG00000011182")
# [1] 31

clper=match(names(myclorderper),names(myorderper)) ## the location of 1093 genes
names(clper)=names(myclorderper)
clper["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 64 
range(clper)
#[1]     1 14041
length(which(clper<1093))
#[1] 247

m=length(myclnames)
ranky=c(1:m)
mypvalueper=rep(0,(m-70+1))
for(j in 70:m){
  mypvalueper[j-69]=wilcox.test(clper[1:j],ranky[1:j])$p.value
}

which.max(mypvalueper<0.05)
#[1] 1
mypvalueper[1]
#[1] 3.05434e-07
wilcox.test(clper[1:9],ranky[1:9])

pernames1=names(clper[1:9])
myper1=myorderper[pernames1]
range(myper1)
#[1] 0.03449179 0.04820186

mydifper=myper[difclnames]
difperorder=order(mydifper)
mydiforderper=mydifper[difperorder]
difclper=match(names(mydiforderper),names(myorderper))  ##the location of 16786 genes  
range(difclper)
#[1]    3 17879
length(which(mydiforderper<max(myper1)))
# [1] 3
pernames2=names(which(mydiforderper<max(myper1)))

myper2=myorderper[pernames2]

#9+3 genes
myfinaldata_per=myGG[,c(pernames1,pernames2)] ## 12genes

################# tam ##########################
tamorder=order(mytam)
myordertam=mytam[tamorder]
which(names(myordertam)=="ENSRNOG00000011182")
# [1] 32

myclnames=rownames(mycldata9_6) ###names of  genes
mycltam=mytam[myclnames]
cltamorder=order(mycltam)
myclordertam=mycltam[cltamorder]
which(names(myclordertam)=="ENSRNOG00000011182")
# [1] 31

cltam=match(names(myclordertam),names(myordertam)) ## the location of 1093 genes
names(cltam)=names(myclordertam)
cltam["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 32 
range(cltam)
#[1]     1 16299
length(which(cltam<1093))
#[1] 644

m=length(myclnames)
ranky=c(1:m)
mypvaluetam=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluetam[j-69]=wilcox.test(cltam[1:j],ranky[1:j])$p.value
}

which.max(mypvaluetam<0.05)
#[1] 195
mypvaluetam[195]
195+69
#[1] 264

tamnames1=names(cltam[1:264])
mytam1=myordertam[tamnames1]
range(mytam1)
#[1] 2.411273 2.668058

mydiftam=mytam[difclnames]
diftamorder=order(mydiftam)
mydifordertam=mydiftam[diftamorder]
difcltam=match(names(mydifordertam),names(myordertam))  ##the location of 16786 genes  
range(difcltam)
#[1]    20 17879
length(which(mydifordertam<max(mytam1)))
# [1] 40
tamnames2=names(which(mydifordertam<max(mytam1)))

mytam2=myordertam[tamnames2]

#264+40 genes
myfinaldata_tam=myGG[,c(tamnames1,tamnames2)] ## 304genes


################# man ##########################
manorder=order(myman)
myorderman=myman[manorder]
which(names(myorderman)=="ENSRNOG00000011182")
# [1] 1

myclnames=rownames(mycldata9_6) ###names of  genes
myclman=myman[myclnames]
clmanorder=order(myclman)
myclorderman=myclman[clmanorder]
which(names(myclorderman)=="ENSRNOG00000011182")
# [1] 1

clman=match(names(myclorderman),names(myorderman)) ## the location of 1093 genes
names(clman)=names(myclorderman)
clman["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 1 
range(clman)
#[1]     1 8464
length(which(clman<1093))
#[1] 855

m=length(myclnames)
ranky=c(1:m)
mypvalueman=rep(0,(m-70+1))
for(j in 70:m){
  mypvalueman[j-69]=wilcox.test(clman[1:j],ranky[1:j])$p.value
}

which.max(mypvalueman<0.05)
#[1] 754
mypvalueman[754]
754+69
#[1] 823

mannames1=names(clman[1:823])
myman1=myorderman[mannames1]
range(myman1)
#[1] 144.0834 389.9579

mydifman=myman[difclnames]
difmanorder=order(mydifman)
mydiforderman=mydifman[difmanorder]
difclman=match(names(mydiforderman),names(myorderman))  ##the location of 16786 genes  
range(difclman)
#[1]    174 17879
length(which(mydiforderman<max(myman1)))
# [1] 187
mannames2=names(which(mydiforderman<max(myman1)))

myman2=myorderman[mannames2]

#823+187 genes
myfinaldata_man=myGG[,c(mannames1,mannames2)] ## 1010genes


################# min ##########################
minorder=order(mymin)
myordermin=mymin[minorder]
which(names(myordermin)=="ENSRNOG00000011182")
# [1] 7

myclnames=rownames(mycldata9_6) ###names of  genes
myclmin=mymin[myclnames]
clminorder=order(myclmin)
myclordermin=myclmin[clminorder]
which(names(myclordermin)=="ENSRNOG00000011182")
# [1] 7

clmin=match(names(myclordermin),names(myordermin)) ## the location of 1093 genes
names(clmin)=names(myclordermin)
clmin["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 7 
range(clmin)
#[1]     1 14484
length(which(clmin<1093))
#[1] 737

m=length(myclnames)
ranky=c(1:m)
mypvaluemin=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluemin[j-69]=wilcox.test(clmin[1:j],ranky[1:j])$p.value
}

which.max(mypvaluemin<0.05)
#[1] 515
mypvaluemin[515]
515+69
#[1] 584

minnames1=names(clmin[1:584])
mymin1=myordermin[minnames1]
range(mymin1)
#[1] 4.913826 9.588129

mydifmin=mymin[difclnames]
difminorder=order(mydifmin)
mydifordermin=mydifmin[difminorder]
difclmin=match(names(mydifordermin),names(myordermin))  ##the location of 16786 genes  
range(difclmin)
#[1]    149 17879
length(which(mydifordermin<max(mymin1)))
# [1] 133
minnames2=names(which(mydifordermin<max(mymin1)))

mymin2=myordermin[minnames2]

#584+133 genes
myfinaldata_min=myGG[,c(minnames1,minnames2)] ## 717genes


################# infnorm ##########################
infnormorder=order(myinfnorm)
myorderinfnorm=myinfnorm[infnormorder]
which(names(myorderinfnorm)=="ENSRNOG00000011182")
# [1] 97

myclnames=rownames(mycldata9_6) ###names of  genes
myclinfnorm=myinfnorm[myclnames]
clinfnormorder=order(myclinfnorm)
myclorderinfnorm=myclinfnorm[clinfnormorder]
which(names(myclorderinfnorm)=="ENSRNOG00000011182")
# [1] 96

clinfnorm=match(names(myclorderinfnorm),names(myorderinfnorm)) ## the location of 1093 genes
names(clinfnorm)=names(myclorderinfnorm)
clinfnorm["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 97 
range(clinfnorm)
#[1]     1 15776
length(which(clinfnorm<1093))
#[1] 427

m=length(myclnames)
ranky=c(1:m)
mypvalueinfnorm=rep(0,(m-70+1))
for(j in 70:m){
  mypvalueinfnorm[j-69]=wilcox.test(clinfnorm[1:j],ranky[1:j])$p.value
}

which.max(mypvalueinfnorm<0.05)
#[1] 244
mypvalueinfnorm[244]
244+69
#[1] 313

infnormnames1=names(clinfnorm[1:313])
myinfnorm1=myorderinfnorm[infnormnames1]
range(myinfnorm1)
#[1] 1.770082 3.499937

mydifinfnorm=myinfnorm[difclnames]
difinfnormorder=order(mydifinfnorm)
mydiforderinfnorm=mydifinfnorm[difinfnormorder]
difclinfnorm=match(names(mydiforderinfnorm),names(myorderinfnorm))  ##the location of 16786 genes  
range(difclinfnorm)
#[1]    92 17879
length(which(mydiforderinfnorm<max(myinfnorm1)))
# [1] 174
infnormnames2=names(which(mydiforderinfnorm<max(myinfnorm1)))

myinfnorm2=myorderinfnorm[infnormnames2]

#313+174 genes
myfinaldata_infnorm=myGG[,c(infnormnames1,infnormnames2)] ## 487genes


################# dissim ##########################
dissimorder=order(mydissim)
myorderdissim=mydissim[dissimorder]
which(names(myorderdissim)=="ENSRNOG00000011182")
# [1] 1

myclnames=rownames(mycldata9_6) ###names of  genes
mycldissim=mydissim[myclnames]
cldissimorder=order(mycldissim)
myclorderdissim=mycldissim[cldissimorder]
which(names(myclorderdissim)=="ENSRNOG00000011182")
# [1] 1

cldissim=match(names(myclorderdissim),names(myorderdissim)) ## the location of 1093 genes
names(cldissim)=names(myclorderdissim)
cldissim["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 1 
range(cldissim)
#[1]     1 8481
length(which(cldissim<1093))
#[1] 840

m=length(myclnames)
ranky=c(1:m)
mypvaluedissim=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluedissim[j-69]=wilcox.test(cldissim[1:j],ranky[1:j])$p.value
}

which.max(mypvaluedissim<0.05)
#[1] 711
mypvaluedissim[711]
711+69
#[1] 780

dissimnames1=names(cldissim[1:780])
mydissim1=myorderdissim[dissimnames1]
range(mydissim1)
#[1] 118.5691 349.6542

mydifdissim=mydissim[difclnames]
difdissimorder=order(mydifdissim)
mydiforderdissim=mydifdissim[difdissimorder]
difcldissim=match(names(mydiforderdissim),names(myorderdissim))  ##the location of 16786 genes  
range(difcldissim)
#[1]    281 17879
length(which(mydiforderdissim<max(mydissim1)))
# [1] 173
dissimnames2=names(which(mydiforderdissim<max(mydissim1)))

mydissim2=myorderdissim[dissimnames2]

#780+173 genes
myfinaldata_dissim=myGG[,c(dissimnames1,dissimnames2)] ## 953genes


################# frechet ##########################
frechetorder=order(myfrechet)
myorderfrechet=myfrechet[frechetorder]
which(names(myorderfrechet)=="ENSRNOG00000011182")
# [1] 166

myclnames=rownames(mycldata9_6) ###names of  genes
myclfrechet=myfrechet[myclnames]
clfrechetorder=order(myclfrechet)
myclorderfrechet=myclfrechet[clfrechetorder]
which(names(myclorderfrechet)=="ENSRNOG00000011182")
# [1] 156

clfrechet=match(names(myclorderfrechet),names(myorderfrechet)) ## the location of 1093 genes
names(clfrechet)=names(myclorderfrechet)
clfrechet["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 166 
range(clfrechet)
#[1]     1 15821
length(which(clfrechet<1093))
#[1] 390

m=length(myclnames)
ranky=c(1:m)
mypvaluefrechet=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluefrechet[j-69]=wilcox.test(clfrechet[1:j],ranky[1:j])$p.value
}

which.max(mypvaluefrechet<0.05)
#[1] 209
mypvaluefrechet[209]
209+69
#[1] 278

frechetnames1=names(clfrechet[1:278])
myfrechet1=myorderfrechet[frechetnames1]
range(myfrechet1)
#[1] 1.703942 3.068691

mydiffrechet=myfrechet[difclnames]
diffrechetorder=order(mydiffrechet)
mydiforderfrechet=mydiffrechet[diffrechetorder]
difclfrechet=match(names(mydiforderfrechet),names(myorderfrechet))  ##the location of 16786 genes  
range(difclfrechet)
#[1]    94 17879
length(which(mydiforderfrechet<max(myfrechet1)))
# [1] 135
frechetnames2=names(which(mydiforderfrechet<max(myfrechet1)))

myfrechet2=myorderfrechet[frechetnames2]

#278+135 genes
myfinaldata_frechet=myGG[,c(frechetnames1,frechetnames2)] ## 413genes


################# cort ##########################
cortorder=order(mycort)
myordercort=mycort[cortorder]
which(names(myordercort)=="ENSRNOG00000011182")
# [1] 1

myclnames=rownames(mycldata9_6) ###names of  genes
myclcort=mycort[myclnames]
clcortorder=order(myclcort)
myclordercort=myclcort[clcortorder]
which(names(myclordercort)=="ENSRNOG00000011182")
# [1] 1

clcort=match(names(myclordercort),names(myordercort)) ## the location of 1093 genes
names(clcort)=names(myclordercort)
clcort["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 1 
range(clcort)
#[1]     1 3710
length(which(clcort<1093))
#[1] 835

m=length(myclnames)
ranky=c(1:m)
mypvaluecort=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluecort[j-69]=wilcox.test(clcort[1:j],ranky[1:j])$p.value
}

which.max(mypvaluecort<0.05)
#[1] 722
mypvaluecort[722]
722+69
#[1] 791

cortnames1=names(clcort[1:791])
mycort1=myordercort[cortnames1]
range(mycort1)
#[1] 10.05900 23.05384

mydifcort=mycort[difclnames]
difcortorder=order(mydifcort)
mydifordercort=mydifcort[difcortorder]
difclcort=match(names(mydifordercort),names(myordercort))  ##the location of 16786 genes  
range(difclcort)
#[1]    253 17879
length(which(mydifordercort<max(mycort1)))
# [1] 177
cortnames2=names(which(mydifordercort<max(mycort1)))

mycort2=myordercort[cortnames2]

#791+177 genes
myfinaldata_cort=myGG[,c(cortnames1,cortnames2)] ## 968genes


################# ncd ##########################
ncdorder=order(myncd)
myorderncd=myncd[ncdorder]
which(names(myorderncd)=="ENSRNOG00000011182")
# [1] 4219

myclnames=rownames(mycldata9_6) ###names of  genes
myclncd=myncd[myclnames]
clncdorder=order(myclncd)
myclorderncd=myclncd[clncdorder]
which(names(myclorderncd)=="ENSRNOG00000011182")
# [1] 391

clncd=match(names(myclorderncd),names(myorderncd)) ## the location of 1093 genes
names(clncd)=names(myclorderncd)
clncd["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 4219 
range(clncd)
#[1]     10 17183
length(which(clncd<1093))
#[1] 120

m=length(myclnames)
ranky=c(1:m)
mypvaluencd=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluencd[j-69]=wilcox.test(clncd[1:j],ranky[1:j])$p.value
}

which.max(mypvaluencd<0.05)
#[1] 1
mypvaluencd[1]
wilcox.test(clncd[1:3],ranky[1:3])

ncdnames1=names(clncd[1:3])
myncd1=myorderncd[ncdnames1]
range(myncd1)
#[1] 0.9781860 0.9793956

mydifncd=myncd[difclnames]
difncdorder=order(mydifncd)
mydiforderncd=mydifncd[difncdorder]
difclncd=match(names(mydiforderncd),names(myorderncd))  ##the location of 16786 genes  
range(difclncd)
#[1]    1 17879
length(which(mydiforderncd<max(myncd1)))
# [1] 16
ncdnames2=names(which(mydiforderncd<max(myncd1)))

myncd2=myorderncd[ncdnames2]

#3+16 genes
myfinaldata_ncd=myGG[,c(ncdnames1,ncdnames2)] ## 19genes


################# intper ##########################
intperorder=order(myintper)
myorderintper=myintper[intperorder]
which(names(myorderintper)=="ENSRNOG00000011182")
# [1] 121

myclnames=rownames(mycldata9_6) ###names of  genes
myclintper=myintper[myclnames]
clintperorder=order(myclintper)
myclorderintper=myclintper[clintperorder]
which(names(myclorderintper)=="ENSRNOG00000011182")
# [1] 72

clintper=match(names(myclorderintper),names(myorderintper)) ## the location of 1093 genes
names(clintper)=names(myclorderintper)
clintper["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 121 
range(clintper)
#[1]     6 17660
length(which(clintper<1093))
#[1] 448

m=length(myclnames)
ranky=c(1:m)
mypvalueintper=rep(0,(m-70+1))
for(j in 70:m){
  mypvalueintper[j-69]=wilcox.test(clintper[1:j],ranky[1:j])$p.value
}

which.max(mypvalueintper<0.05)
#[1] 1
mypvalueintper[1]
wilcox.test(clintper[1:3],ranky[1:3])

intpernames1=names(clintper[1:3])
myintper1=myorderintper[intpernames1]
range(myintper1)
#[1] 25.46848 33.38029

mydifintper=myintper[difclnames]
difintperorder=order(mydifintper)
mydiforderintper=mydifintper[difintperorder]
difclintper=match(names(mydiforderintper),names(myorderintper))  ##the location of 16786 genes  
range(difclintper)
#[1]    1 17879
length(which(mydiforderintper<max(myintper1)))
# [1] 6
intpernames2=names(which(mydiforderintper<max(myintper1)))

myintper2=myorderintper[intpernames2]

#3+6 genes
myfinaldata_intper=myGG[,c(intpernames1,intpernames2)] ## 9genes


################# cdm ##########################
cdmorder=order(mycdm)
myordercdm=mycdm[cdmorder]
which(names(myordercdm)=="ENSRNOG00000011182")
# [1] 4952

myclnames=rownames(mycldata9_6) ###names of  genes
myclcdm=mycdm[myclnames]
clcdmorder=order(myclcdm)
myclordercdm=myclcdm[clcdmorder]
which(names(myclordercdm)=="ENSRNOG00000011182")
# [1] 396

clcdm=match(names(myclordercdm),names(myordercdm)) ## the location of 1093 genes
names(clcdm)=names(myclordercdm)
clcdm["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 4952 
range(clcdm)
#[1]  26 17266
length(which(clcdm<1093))
#[1] 95

m=length(myclnames)
ranky=c(1:m)
mypvaluecdm=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluecdm[j-69]=wilcox.test(clcdm[1:j],ranky[1:j])$p.value
}

which.max(mypvaluecdm<0.05)
#[1] 1
mypvaluecdm[1]
wilcox.test(clcdm[1:3],ranky[1:3])

cdmnames1=names(clcdm[1:3])
mycdm1=myordercdm[cdmnames1]
range(mycdm1)
#[1] 0.9883916 0.9894722

mydifcdm=mycdm[difclnames]
difcdmorder=order(mydifcdm)
mydifordercdm=mydifcdm[difcdmorder]
difclcdm=match(names(mydifordercdm),names(myordercdm))  ##the location of 16786 genes  
range(difclcdm)
#[1]    1 17879
length(which(mydifordercdm<max(mycdm1)))
# [1] 65
cdmnames2=names(which(mydifordercdm<max(mycdm1)))

mycdm2=myordercdm[cdmnames2]

#3+65 genes
myfinaldata_cdm=myGG[,c(cdmnames1,cdmnames2)] ## 68genes


################# spec ##########################
specorder=order(myspec)
myorderspec=myspec[specorder]
which(names(myorderspec)=="ENSRNOG00000011182")
# [1] 63

myclnames=rownames(mycldata9_6) ###names of  genes
myclspec=myspec[myclnames]
clspecorder=order(myclspec)
myclorderspec=myclspec[clspecorder]
which(names(myclorderspec)=="ENSRNOG00000011182")
# [1] 39

clspec=match(names(myclorderspec),names(myorderspec)) ## the location of 1093 genes
names(clspec)=names(myclorderspec)
clspec["ENSRNOG00000011182"]
# ENSRNOG00000011182 
# 63 
range(clspec)
#[1]   18 15084
length(which(clspec<1093))
#[1] 506

m=length(myclnames)
ranky=c(1:m)
mypvaluespec=rep(0,(m-70+1))
for(j in 70:m){
  mypvaluespec[j-69]=wilcox.test(clspec[1:j],ranky[1:j])$p.value
}

which.max(mypvaluespec<0.05)
#[1] 1
mypvaluespec[1]
wilcox.test(clspec[1:3],ranky[1:3])

specnames1=names(clspec[1:3])
myspec1=myorderspec[specnames1]
range(myspec1)
#[1] 1855.494 1875.220

mydifspec=myspec[difclnames]
difspecorder=order(mydifspec)
mydiforderspec=mydifspec[difspecorder]
difclspec=match(names(mydiforderspec),names(myorderspec))  ##the location of 16786 genes  
range(difclspec)
#[1]    1 17879
length(which(mydiforderspec<max(myspec1)))
# [1] 17
specnames2=names(which(mydiforderspec<max(myspec1)))

myspec2=myorderspec[specnames2]

#3+17 genes
myfinaldata_spec=myGG[,c(specnames1,specnames2)] ## 20genes

