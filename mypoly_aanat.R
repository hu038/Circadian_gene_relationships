### continue with "myrawdata.R"

G11182=GG["ENSRNOG00000011182",]
G_total[G_total[,1]=="ENSRNOG00000011182",][483]
plot(t(G11182),type="l",xlab="Time points",ylab="Expression value",
     main="Aanat")
sd(G11182)

library(fANCOVA)

index <- 1:480
GG11182=data.frame(t(G11182),index)
loessMod <- loess(ENSRNOG00000011182 ~ index, degree=2, data=GG11182)
smoothed <- predict(loessMod)


loessMod8 <- loess(ENSRNOG00000011182 ~ index, degree=2, data=GG11182, span=0.07970886)
loessMod10 <- loess(ENSRNOG00000011182 ~ index, degree=2,data=GG11182, span=0.10) # 10% smoothing span
loessMod25 <- loess(ENSRNOG00000011182 ~ index,degree=2, data=GG11182, span=0.25) # 25% smoothing span
loessMod50 <- loess(ENSRNOG00000011182 ~ index,degree=2, data=GG11182, span=0.50) # 50% smoothing span


# get smoothed output
smoothed8 <- predict(loessMod8)
smoothed10 <- predict(loessMod10) 
smoothed25 <- predict(loessMod25) 
smoothed50 <- predict(loessMod50)

# Plot it
plot(t(G11182),type="l",xlab="Time points",ylab="Expression value",
     main="Aanat")
lines(smoothed8,  col=2,lwd=3)
lines(smoothed10,  col=3,lwd=3)
lines(smoothed25, col=4,lwd=3)
lines(smoothed50,  col=5,lwd=3)
lines(smoothed, col=6,lwd=3)
legend("topleft",legend=c("8%","10%","25%","50%","75%"),col=c(2,3,4,5,6),
       lty=rep(1,5))

FTSE.lo3 <- loess.as(index, t(G11182), degree = 2, 
                     criterion = c("aicc", "gcv")[2], user.span = NULL, plot = T)
### summary(FTSE.lo3)
### span     :  0.05693009 degree   :  1
FTSE.lo.predict3 <- predict(FTSE.lo3, data.frame(Index=index))
plot(t(G11182),type="l",xlab="Time points",ylab="Expression value",
     main="Aanat")
lines(FTSE.lo.predict3, col=2,lwd=3)
myployaanat=FTSE.lo.predict3


