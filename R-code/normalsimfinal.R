library("VineCopula")
library("beepr")
set.seed(336)
#parameters
#parameters for A
delta_0 <- -0.4
delta_X <- 0.6
#parameters for C
alpha_0 <- 1
alpha_A <- 2
alpha_X <- 2
sdC.mar <- 4
#parameters for M
gamma_0 <- 3
gamma_A <- 2
gamma_X <- 4
sdM.mar <- 5
#parameters for Y
beta_0 <- 1
beta_A <- 2
beta_M <- 2
beta_C <- 3
beta_X <- 5
beta_AC <- 4
beta_AM <- 2
sdy <- 4

#copula
copula <- BiCop(1,0.2)


expit <- function(x) exp(x)/(1+exp(x))
library(MASS)

quad <- function(fitm.mar,fitc.mar,copulafamily=1,fity,fita,A,M,C,X,Y)
{
  ##M and C
  fEM.mar <- predict(fitm.mar)
  fsdm.mar <- sigma(fitm.mar)
  dm.mar <- dnorm(M,mean = fEM.mar,sd = fsdm.mar)
  u1 <- pnorm(M,mean = fEM.mar,sd=fsdm.mar)
  
  fEC.mar <- predict(fitc.mar)
  fsdc.mar <- sigma(fitc.mar)
  dc.mar <- dnorm(C,mean = fEC.mar, sd=fsdc.mar)
  u2 <- pnorm(C,mean = fEC.mar,sd=fsdc.mar)
  rho_est <- BiCopEst(u1, u2, family = copulafamily, method = "mle")$par
  copula <- BiCop(copulafamily,rho_est)
  copulaPDF <- BiCopPDF(u1, u2,copula)
  
  #f(M|A=0,X)
  datam0x <- data.frame(cbind(0,X))
  colnames(datam0x) <- c("A","X")
  EM0x <- predict(fitm.mar,newdata = datam0x)
  sdm0x <- sigma(fitm.mar)
  dm0x <- dnorm(M,mean = EM0x,sdm0x)
  
  ##A and Y
  proba <- A*expit(predict(fita))+(1-A)*(1-expit(predict(fita)))
  Ey <- predict(fity)
  
  #estimate
  #delta1
  ratio1 <- ((2*A-1)/proba)*(dm0x/(dm.mar*copulaPDF))
  est_delta1 <- ratio1*Y
  est_comp1 <- ratio1*(Y-Ey)
  
  #delta2
  #eta
  etafun <- function(a)
  {
    datam <- data.frame(cbind(0,X))
    colnames(datam) <- c("A","X")
    EM <- predict(fitm.mar,newdata = datam)
    data <- data.frame(cbind(a,C,EM,X))
    colnames(data) <- c("A","C","M","X")
    eta <- predict(fity,newdata = data)
    eta
  }
  ratio2 <- (2*A-1)/proba
  est_delta2 <- est_comp2 <- ratio2*etafun(A)
  
  #delta3
  gammafun <-function(a)
  {
    datac <- data.frame(cbind(a,X))
    colnames(datac) <- c("A","X")
    EC <- predict(fitc.mar,newdata = datac)
    data <- data.frame(cbind(a,EC,M,X))
    colnames(data) <- c("A","C","M","X")
    gamma <- predict(fity,newdata = data)
    gamma
  }
  ratio3 <- (1-A)/proba
  est_delta3 <- est_comp3 <- ratio3*(gammafun(1)-gammafun(0))
  
  #delta4
  taufun <- function(a)
  {
    datac <- data.frame(cbind(a,X))
    colnames(datac) <- c("A","X")
    EC <- predict(fitc.mar,newdata=datac)
    datam <- data.frame(cbind(0,X))
    colnames(datam) <- c("A","X")
    EM <- predict(fitm.mar,newdata=datam)
    data <- data.frame(cbind(a,EC,EM,X))
    colnames(data) <- c("A","C","M","X")
    tau <- predict(fity,newdata = data)
    tau
  }
  est_delta4 <- taufun(1)-taufun(0)
  est_comp4 <- (1-ratio3)*est_delta4
  
  est_comp5 <- -ratio2*taufun(A)
  
  quad <- mean(est_comp1 + est_comp2 + est_comp3 + est_comp4 + est_comp5)
  
  result <- data.frame(mean(est_delta1),mean(est_delta2),mean(est_delta3),mean(est_delta4),mean(quad))
  
  result
}

gen <- function(n,caltruth=T)
{
  X <- rnorm(n)
  prob_A <- expit(delta_0+delta_X*X)
  A <- rbinom(n,1,prob = prob_A)
  
  #marginals
  Ec_ax <- alpha_0+alpha_A*A+alpha_X*X
  Varc_ax <- sdC.mar^2
  sdc_ax <- sqrt(Varc_ax)
  Em_ax <- gamma_0+gamma_A*A+gamma_X*X
  Varm_ax <- sdM.mar^2
  sdm_ax <- sqrt(Varm_ax)
  
  #####
  a <- BiCopSim(n,copula)
  C <- qnorm(a[,1],mean=Ec_ax,sd=sdc_ax)
  M <- qnorm(a[,2],mean=Em_ax,sd=sdm_ax)
  
  mu_y <- beta_0 + beta_A*A + beta_M*M + beta_C*C + beta_X*X + beta_AC * A*C + beta_AM*A*M
  Y <- rnorm(n,mu_y,sdy)
  
  ###truth
  M0X <- gamma_0 + gamma_X*X
  C1X <- alpha_0 + alpha_A + alpha_X*X 
  C0X <- alpha_0 + alpha_X*X
  EY1cmx <- beta_0 + beta_A + beta_M*M0X + beta_C*C1X + beta_X*X + beta_AC * C1X + beta_AM *M0X
  EY0cmx <- beta_0 + beta_M*M0X + beta_C*C0X + beta_X*X
  truth <- mean(EY1cmx) - mean(EY0cmx)
  if(caltruth)
  {
    return(truth)
  }
  if(!caltruth)
  {
    return(data.frame(Y,A,X,M,C))
  }
}

#all correct
#quad(fitm.mar,fitc.mar,1,fity,fita,A,M,C,X,Y)-truth
#wrong EY

do.one<-function(n)
{
  truth <- gen(10000,T)
  dataset <- dataset.wrong <- gen(n,F)
  X2 <- rnorm(n)
  dataset.wrong$X <- X2
  
  fitm.mar <- lm(M~A+X,data = dataset)
  fitc.mar <- lm(C~A+X,data = dataset)
  fita <- glm(A~X,family = binomial(link = "logit"),data = dataset)
  fity <- lm(Y~A+X+C+M+A*C+A*M,data = dataset)
  
  quad(fitm.mar,fitc.mar,1,fity,fita,dataset$A,dataset$M,dataset$C,dataset$X,dataset$Y)-truth
}
B <- 1000
n <- 1500
result.allcorrect <- replicate(B,do.one(n))
result2.allcorrect <- matrix(unlist(result.allcorrect),ncol=B)
correct <- apply(result2.allcorrect,1,mean)*100
beep()


do.one.M1 <-function(n)
{
  truth <- gen(10000,T)
  dataset <- dataset.wrong <- gen(n,F)
  X2 <- rnorm(n)
  dataset.wrong$X <- X2
  
  fitm.mar <- lm(M~A+X,data = dataset)
  fitc.mar <- lm(C~A+X,data = dataset)
  fita <- glm(A~X,family = binomial(link = "logit"),data = dataset)
  fity <- lm(Y~A+X+C+M+A*C+A*M,data = dataset.wrong)

  quad(fitm.mar,fitc.mar,1,fity,fita,dataset$A,dataset$M,dataset$C,dataset$X,dataset$Y)-truth
}

result.M1 <- replicate(B,do.one.M1(n))
result2.M1 <- matrix(unlist(result.M1),ncol=B)
M1 <- apply(result2.M1,1,mean)*100
beep()

do.one.M2 <-function(n)
{
  truth <- gen(10000,T)
  dataset <- dataset.wrong <- gen(n,F)
  X2 <- rnorm(n)
  dataset.wrong$X <- X2
  
  fitm.mar <- lm(M~A+X,data = dataset)
  fitc.mar <- lm(C~A+X,data = dataset.wrong)
  fita <- glm(A~X,family = binomial(link = "logit"),data = dataset)
  fity <- lm(Y~A+X+C+M+A*C+A*M,data = dataset)

  quad(fitm.mar,fitc.mar,1,fity,fita,dataset$A,dataset$M,dataset$C,dataset$X,dataset$Y)-truth
}

result.M2 <- replicate(B,do.one.M2(n))
result2.M2 <- matrix(unlist(result.M2),ncol=B)
M2 <- apply(result2.M2,1,mean)*100
beep()

do.one.M3<-function(n)
{
  truth <- gen(10000,T)
  dataset <- dataset.wrong <- gen(n,F)
  X2 <- rnorm(n)
  dataset.wrong$X <- X2
  
  fitm.mar <- lm(M~A+X,data = dataset.wrong)
  fitc.mar <- lm(C~A+X,data = dataset)
  fita <- glm(A~X,family = binomial(link = "logit"),data = dataset)
  fity <- lm(Y~A+X+C+M+A*C+A*M,data = dataset)
  
  quad(fitm.mar,fitc.mar,1,fity,fita,dataset$A,dataset$M,dataset$C,dataset$X,dataset$Y)-truth
}

result.M3 <- replicate(B,do.one.M3(n))
result2.M3 <- matrix(unlist(result.M3),ncol=B)
M3 <- apply(result2.M3,1,mean)*100
beep()

do.one.M4<-function(n)
{
  truth <- gen(10000,T)
  dataset <- dataset.wrong <- gen(n,F)
  X2 <- rnorm(n)
  dataset.wrong$X <- X2
  
  fitm.mar <- lm(M~A+X,data = dataset)
  fitc.mar <- lm(C~A+X,data = dataset)
  fita <- glm(A~X,family = binomial(link = "logit"),data = dataset.wrong)
  fity <- lm(Y~A+X+C+M+A*C+A*M,data = dataset)
  
  quad(fitm.mar,fitc.mar,1,fity,fita,dataset$A,dataset$M,dataset$C,dataset$X,dataset$Y)-truth
}

result.M4 <- replicate(B,do.one.M4(n))
result2.M4 <- matrix(unlist(result.M4),ncol=B)
M4 <- apply(result2.M4,1,mean)*100
beep()


do.one.Mnon<-function(n)
{
  truth <- gen(10000,T)
  dataset <- dataset.wrong <- gen(n,F)
  X2 <- rnorm(n)
  dataset.wrong$X <- X2
  
  fitm.mar <- lm(M~A+X,data = dataset.wrong)
  fitc.mar <- lm(C~A+X,data = dataset.wrong)
  fita <- glm(A~X,family = binomial(link = "logit"),data = dataset.wrong)
  fity <- lm(Y~A+X+C+M+A*C+A*M,data = dataset.wrong)
  
  quad(fitm.mar,fitc.mar,1,fity,fita,dataset$A,dataset$M,dataset$C,dataset$X,dataset$Y)-truth
}

result.Mnon <- replicate(B,do.one.Mnon(n))
result2.Mnon <- matrix(unlist(result.Mnon),ncol=B)
Mnon <- apply(result2.Mnon,1,mean)*100
beep()

correct.se <- apply(result2.allcorrect,1,sd)
M1.se <- apply(result2.M1,1,sd)
M2.se <- apply(result2.M2,1,sd)
M3.se <- apply(result2.M3,1,sd)
M4.se <- apply(result2.M4,1,sd)

table1 <- rbind(correct,M1,M2,M3,M4)
table1.se <- rbind(correct.se,M1.se,M2.se,M3.se,M4.se)*100
write.csv(table1,"normalsim1500final.csv")
write.csv(table1.se,"normalsimse1500final.csv")

library(xtable)
xtable(table1,digits = 0)
xtable(table1.se,digits = 0)


