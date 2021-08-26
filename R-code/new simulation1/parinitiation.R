library("VineCopula")
library("beepr")
library("MASS")
set.seed(336)

#helper functions
expit <- function(x) exp(x)/(1+exp(x))
n <- 5
#parameters
tau_0 <- 1.2

#parameters for A
a0 <- 0.3
aU <- c(0.3,-0.2,0.4) 
aint12 <- -0.1
aint23 <- 0.3
aint13 <- -0.2

#marginal for mediators
#C marginal parameters
c0 <- 0.3
cA <- -0.2
cU <- c(-0.1,0.1,0.2)
cint12 <- -0.2
cint23 <- -0.2
cint13 <- 0.1


#M marginal parameters
m0 <- -0.3
mA <- -0.3
mU <- c(0.1,0.3,-0.2)
mint12 <- 0.2
mint23 <- -0.2
mint13 <- -0.1





#generate Y
y0 <- -10
yA <- 6
yM <- 5
yC <- -3
yU <- c(3,-2,2)
yAC <- 3
yAM <- -3
yint12 <- 2
yint23 <- 2
yint13 <- 1

sdy <- 5


#calculate the true value empirically
n <- 5000
#calculating the true value exactly
U1 <- U2 <- U3 <- rep(NA,n)
for(i in 1:n)
{
  U.temp <- mvrnorm(1,mu=cbind(0,0,0), Sigma=diag(c(1,1,1)))
  U1[i] <- U.temp[1]
  U2[i] <- U.temp[2]
  U3[i] <- U.temp[3]
}
U.row <- cbind(U1,U2,U3)
U <- t(U.row)


Ey <- function(a,m,c)
{
  return(y0 + yA*a + yM*m + yC*c + yU %*% U + yAC * a*c + yAM*a*m + yint12 * U1 * U2 + yint23 * U2 * U3 +
           yint13 * U1 * U3)
}

EY1cmu <- mean(Ey(1,1,1)*expit(m0+mU%*%U+ mint12 * U1 * U2 + mint13 * U1 * U3 +
                                 mint23 * U2 * U3)*expit(c0+cA+cU%*%U +cint12 * U1 * U2 + cint13 * U1 * U3 +
                                                           cint23 * U2 * U3)+
                 Ey(1,1,0)*expit(m0+mU%*%U+ mint12 * U1 * U2 + mint13 * U1 * U3 +
                                   mint23 * U2 * U3)*(1-expit(c0+cA+cU%*%U+cint12 * U1 * U2 + cint13 * U1 * U3 +
                                                                cint23 * U2 * U3))+
                 Ey(1,0,1)*(1-expit(m0+mU%*%U+ mint12 * U1 * U2 + mint13 * U1 * U3 +
                                      mint23 * U2 * U3))*expit(c0+cA+cU%*%U+cint12 * U1 * U2 + cint13 * U1 * U3 +
                                                                 cint23 * U2 * U3)+
                 Ey(1,0,0)*(1-expit(m0+mU%*%U+ mint12 * U1 * U2 + mint13 * U1 * U3 +
                                      mint23 * U2 * U3))*(1-expit(c0+cA+cU%*%U+cint12 * U1 * U2 + cint13 * U1 * U3 +
                                                                    cint23 * U2 * U3))
)

EY0cmu <- mean(Ey(0,1,1)*expit(m0+mU%*%U+ mint12 * U1 * U2 + mint13 * U1 * U3 +
                                 mint23 * U2 * U3)*expit(c0+cU%*%U+ cint12 * U1 * U2 + cint13 * U1 * U3 +
                                                           cint23 * U2 * U3)+
                 Ey(0,1,0)*expit(m0+mU%*%U+ mint12 * U1 * U2 + mint13 * U1 * U3 +
                                   mint23 * U2 * U3)*(1-expit(c0+cU%*%U+ cint12 * U1 * U2 + cint13 * U1 * U3 +
                                                                cint23 * U2 * U3))+
                 Ey(0,0,1)*(1-expit(m0+mU%*%U+ mint12 * U1 * U2 + mint13 * U1 * U3 +
                                      mint23 * U2 * U3))*expit(c0+cU%*%U+ cint12 * U1 * U2 + cint13 * U1 * U3 +
                                                                 cint23 * U2 * U3)+
                 Ey(0,0,0)*(1-expit(m0+mU%*%U+ mint12 * U1 * U2 + mint13 * U1 * U3 +
                                      mint23 * U2 * U3))*(1-expit(c0+cU%*%U+ cint12 * U1 * U2 + cint13 * U1 * U3 +
                                                                    cint23 * U2 * U3))
)

truth<- EY1cmu - EY0cmu

truth
