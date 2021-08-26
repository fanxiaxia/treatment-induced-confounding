library("VineCopula")
library("beepr")
library("MASS")
set.seed(336)


#helper functions
expit <- function(x) exp(x)/(1+exp(x))
n <- 5000
#parameters

#parameters for A
a0 <- 0.4
aU <- c(0.3,-0.2,0.4) 
aint12 <- 0.1
aint23 <- 0.1
aint13 <- 0.2
#marginal for mediators
#C marginal parameters
c0 <- 1
cA <- -3
cU <- c(-1,2,1)
cint12 <- 3
cint23 <- 2
cint13 <- 3

sdc.mar <- 1
Varc_au <- sdc.mar^2

#M marginal parameters
m0 <- -1
mA <- 2
mU <- c(2,-3,2)
mint12 <- 2
mint23 <- 3
mint13 <- 2

sdm.mar <- 5
Varm_au <- sdm.mar^2



#generate Y
y0 <- -2
yA <- 2
yM <- 1
yC <- -2
yU <- c(1,2,2)
yAC <- 3
yAM <- -2
yint12 <- 3
yint23 <- -1
yint13 <- 2

sdy <- 1

par.a <- data.frame(a0,aU,aint12,aint23,aint13)
par.c <- data.frame(c0,cA,cU,cint12,cint13,cint23,sdc.mar,Varc_au)
par.m <- data.frame(m0,mA,mU,mint12,mint13,mint23,sdm.mar,Varm_au)
par.y <- data.frame(y0,yA,yM,yC,yU,yAC,yAM,yint12,yint13,yint23,sdy)

corr <- 0.2
par.set <- data.frame(cbind(pu1,par.a,par.c,par.m,par.y,corr))


#calculating the true value exactly

truth.exact <- yA + yC*(cA)+yAC*(c0+cA)+yAM*(m0)


truth.exact 
