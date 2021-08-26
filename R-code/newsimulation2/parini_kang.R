library("VineCopula")
library("beepr")
library("MASS")
set.seed(336)
corr <- 0.2
#helper functions
expit <- function(x) exp(x)/(1+exp(x))
n <- 5000
#parameters
tau_0 <- 1.2
tau_A <- 0
tau_X <- c(0,0,0,0)
#baseline covariates
#parameters for A
aU <- c(-1,0.5,-0.25,-0.1) 
#marginal for mediators
#C marginal parameters
c0 <- -1.6
cA <- 2
cU <- c(1,-0.8,0.6,-1)


#M marginal parameters
m0 <- -1.5
mA <- 2
mU <- c(1,-0.5,0.9,-1)


#generate Y
y0 <- 210
yA <- 1
yM <- 1
yC <--50
yU <- c(27.4,13.7,13.7,13.7)
sdy <- 30

n <- 5000
#calculating the true value exactly
U1 <- U2 <- U3 <- U4 <- rep(NA,n)
for(i in 1:n)
{
  U.temp <- mvrnorm(1,mu=cbind(0,0,0,0), Sigma=diag(c(1,1,1,1)))
  U1[i] <- U.temp[1]
  U2[i] <- U.temp[2]
  U3[i] <- U.temp[3]
  U4[i] <- U.temp[4]
}
U.row <- cbind(U1,U2,U3,U4)
U <- t(U.row)


Ey <- function(a,m,c)
{
  return(y0+yA*a+yM*m+yC*c+yU %*% U)
}

EY1cmu <- mean(Ey(1,1,1)*expit(m0+mU%*%U)*expit(c0+cA+cU%*%U)+
                 Ey(1,1,0)*expit(m0+mU%*%U)*(1-expit(c0+cA+cU%*%U))+
                 Ey(1,0,1)*(1-expit(m0+mU%*%U))*expit(c0+cA+cU%*%U)+
                 Ey(1,0,0)*(1-expit(m0+mU%*%U))*(1-expit(c0+cA+cU%*%U))
                 )

EY0cmu <- mean(Ey(0,1,1)*expit(m0+mU%*%U)*expit(c0+cU%*%U)+
                 Ey(0,1,0)*expit(m0+mU%*%U)*(1-expit(c0+cU%*%U))+
                 Ey(0,0,1)*(1-expit(m0+mU%*%U))*expit(c0+cU%*%U)+
                 Ey(0,0,0)*(1-expit(m0+mU%*%U))*(1-expit(c0+cU%*%U))
)

truth<- EY1cmu - EY0cmu

truth


dataset.Z <- data_genZ(n)
Xz <- cbind(dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
Xu <- cbind(dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
#
fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
fita <- glm(A~Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
fity.tru <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)

cor(fitted(fity),fitted(fity.tru))

############
data.a <- function(a)
{
  datac0m0 <<- data.frame(cbind(a,0,0,Xz))
  colnames(datac0m0) <<- c("A","C","M","Z1","Z2","Z3","Z4")
  datac1m1 <<- data.frame(cbind(a,1,1,Xz))
  colnames(datac1m1) <<- c("A","C","M","Z1","Z2","Z3","Z4")
  datac1m0 <<- data.frame(cbind(a,1,0,Xz))
  colnames(datac1m0) <<- c("A","C","M","Z1","Z2","Z3","Z4")
  datac0m1 <<- data.frame(cbind(a,0,1,Xz))
  colnames(datac0m1) <<- c("A","C","M","Z1","Z2","Z3","Z4")
  
  datac0m0.u <<- data.frame(cbind(a,0,0,Xu))
  colnames(datac0m0.u) <<- c("A","C","M","U1","U2","U3","U4")
  datac1m1.u <<- data.frame(cbind(a,1,1,Xu))
  colnames(datac1m1.u) <<- c("A","C","M","U1","U2","U3","U4")
  datac1m0.u <<- data.frame(cbind(a,1,0,Xu))
  colnames(datac1m0.u) <<- c("A","C","M","U1","U2","U3","U4")
  datac0m1.u <<- data.frame(cbind(a,0,1,Xu))
  colnames(datac0m1.u) <<- c("A","C","M","U1","U2","U3","U4")
  
  return("I hate my life")
}
data.a(1)
cor(predict(fity,newdata = datac0m0),predict(fity.tru,newdata = datac0m0.u))
cor(predict(fity,newdata = datac0m1),predict(fity.tru,newdata = datac0m1.u))
cor(predict(fity,newdata = datac1m0),predict(fity.tru,newdata = datac1m0.u))
cor(predict(fity,newdata = datac1m1),predict(fity.tru,newdata = datac1m1.u))

data.a(0)
cor(predict(fity,newdata = datac0m0),predict(fity.tru,newdata = datac0m0.u))
cor(predict(fity,newdata = datac0m1),predict(fity.tru,newdata = datac0m1.u))
cor(predict(fity,newdata = datac1m0),predict(fity.tru,newdata = datac1m0.u))
cor(predict(fity,newdata = datac1m1),predict(fity.tru,newdata = datac1m1.u))



data.a(A)
cor(predict(fity,newdata = datac0m0),predict(fity.tru,newdata = datac0m0.u))
cor(predict(fity,newdata = datac0m1),predict(fity.tru,newdata = datac0m1.u))
cor(predict(fity,newdata = datac1m0),predict(fity.tru,newdata = datac1m0.u))
cor(predict(fity,newdata = datac1m1),predict(fity.tru,newdata = datac1m1.u))


set.seed(310)
library(beepr)
library(parallel)
library(nleqslv)
logit <- function(x){log(x/(1-x))}
result.Z.wrong <- mclapply(trials,do.oneZ.wrong,mc.cores=numCores)
result2.Z.wrong <- matrix(unlist(result.Z.wrong),ncol=B)
resZ.wrong <- apply(result2.Z.wrong,1,mean)*100
resZ.wrong
resZ.wrong.sd <- apply(result2.Z.wrong,1,sd)*100
resZ.wrong.sd
beep()



