set.seed(233)
library(beepr)
library(parallel)
library(nleqslv)
logit <- function(x){log(x/(1-x))}
B <- 1000
n <- 1000
numCores <- detectCores()
trials <- seq(1:B)
######do not observe U
do.oneZ.corr <- function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)

  ###########  
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc

  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                mu=T,cu=T,au=T,yu=T,p11.est.stable)[5]-truth
  #print(summary(fita))
  
  ###
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=T,au=T,yu=T,p11.est.stable)[5]-truth
  
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=T,au=T,yu=T,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.corr <- replicate(B,do.oneZ.corr(n))
result.Z.corr <- mclapply(trials,do.oneZ.corr,mc.cores=numCores)
result2.Z.corr <- matrix(unlist(result.Z.corr),nrow=7)
resZ.corr <- c(apply(result2.Z.corr,1,mean)*100)
resZ.corr
resZ.corr.sd <- apply(result2.Z.corr,1,sd)*100
beep()


do.oneZ.m<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  ######
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                        mu=F,cu=T,au=T,yu=T,p11.est.stable)[5]-truth
  #######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=T,au=T,yu=T,p11.est.stable)[5]-truth
  
  ######
  corr <- cor(fitted(fitm.mar),fitted(glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=T,au=T,yu=T,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.m <- replicate(B,do.oneZ.m(n))
result.Z.m <- mclapply(trials,do.oneZ.m,mc.cores=numCores)
result2.Z.m <- matrix(unlist(result.Z.m),ncol=B)
resZ.m <- c(apply(result2.Z.m,1,mean)*100)
resZ.m
resZ.m.sd <- apply(result2.Z.m,1,sd)*100
beep()


do.oneZ.mc <-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  ######
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=F,au=T,yu=T,p11.est.stable)[5]-truth
  #######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=F,au=T,yu=T,p11.est.stable)[5]-truth
  
  ######
#  corr <- cor(fitted(fitm.mar),fitted(glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=F,au=T,yu=T,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.m <- replicate(B,do.oneZ.m(n))
result.Z.mc <- mclapply(trials,do.oneZ.mc,mc.cores=numCores)
result2.Z.mc <- matrix(unlist(result.Z.mc),ncol=B)
resZ.mc <- c(apply(result2.Z.mc,1,mean)*100)
resZ.mc
resZ.mc.sd <- apply(result2.Z.mc,1,sd)*100
beep()


do.oneZ.mca <-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  ######
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=F,au=F,yu=T,p11.est.stable)[5]-truth
  #######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=F,au=F,yu=T,p11.est.stable)[5]-truth
  
  ######
  #  corr <- cor(fitted(fitm.mar),fitted(glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=F,au=F,yu=T,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.m <- replicate(B,do.oneZ.m(n))
result.Z.mca <- mclapply(trials,do.oneZ.mca,mc.cores=numCores)
result2.Z.mca <- matrix(unlist(result.Z.mca),ncol=B)
resZ.mca <- c(apply(result2.Z.mca,1,mean)*100)
resZ.mca
resZ.mca.sd <- apply(result2.Z.mca,1,sd)*100
beep()



do.oneZ.mcy <-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  
  ######
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=F,au=T,yu=F,p11.est.stable)[5]-truth
  #######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=F,au=T,yu=F,p11.est.stable)[5]-truth
  
  ######
  #  corr <- cor(fitted(fitm.mar),fitted(glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=F,au=T,yu=F,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.m <- replicate(B,do.oneZ.m(n))
result.Z.mcy <- mclapply(trials,do.oneZ.mcy,mc.cores=numCores)
result2.Z.mcy <- matrix(unlist(result.Z.mcy),nrow=7)
resZ.mcy <- apply(result2.Z.mcy,1,mean)*100
resZ.mcy
resZ.mcy.sd <- apply(result2.Z.mcy,1,sd)*100
beep()



do.oneZ.ma<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  ######
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=T,au=F,yu=T,p11.est.stable)[5]-truth
  #######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=T,au=F,yu=T,p11.est.stable)[5]-truth
  
  ######
  corr <- cor(fitted(fitm.mar),fitted(glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=T,au=F,yu=T,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.m <- replicate(B,do.oneZ.m(n))
result.Z.ma <- mclapply(trials,do.oneZ.ma,mc.cores=numCores)
result2.Z.ma <- matrix(unlist(result.Z.ma),ncol=B)
resZ.ma <- c(apply(result2.Z.ma,1,mean)*100)
resZ.ma
resZ.ma.sd <- apply(result2.Z.ma,1,sd)*100
beep()


do.oneZ.may<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  
  ######
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=T,au=F,yu=F,p11.est.stable)[5]-truth
  #######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=T,au=F,yu=F,p11.est.stable)[5]-truth
  
  ######
  corr <- cor(fitted(fitm.mar),fitted(glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=T,au=F,yu=F,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.m <- replicate(B,do.oneZ.m(n))
result.Z.may <- mclapply(trials,do.oneZ.may,mc.cores=numCores)
result2.Z.may <- matrix(unlist(result.Z.may),ncol=B)
resZ.may <- c(apply(result2.Z.may,1,mean)*100)
resZ.may
resZ.may.sd <- apply(result2.Z.may,1,sd)*100
beep()



do.oneZ.my<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  
  ######
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=T,au=T,yu=F,p11.est.stable)[5]-truth
  #######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=T,au=T,yu=F,p11.est.stable)[5]-truth
  
  ######
  corr <- cor(fitted(fitm.mar),fitted(glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=T,au=T,yu=F,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.m <- replicate(B,do.oneZ.m(n))
result.Z.my <- mclapply(trials,do.oneZ.my,mc.cores=numCores)
result2.Z.my <- matrix(unlist(result.Z.my),ncol=B)
resZ.my <- c(apply(result2.Z.my,1,mean)*100)
resZ.my
resZ.my.sd <- apply(result2.Z.my,1,sd)*100
beep()


do.oneZ.c<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  ##########
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                        mu=T,cu=F,au=T,yu=T,p11.est.stable)[5]-truth
  
  ######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=F,au=T,yu=T,p11.est.stable)[5]-truth
  
  corr <- cor(fitted(fitc.mar),fitted(glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=F,au=T,yu=T,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.c <- replicate(B,do.oneZ.c(n))
result.Z.c <- mclapply(trials,do.oneZ.c,mc.cores=numCores)
result2.Z.c <- matrix(unlist(result.Z.c),ncol=B)
resZ.c <- c(apply(result2.Z.c,1,mean)*100)
resZ.c
resZ.c.sd <- apply(result2.Z.c,1,sd)*100
beep()

do.oneZ.ca<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  ##########
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=F,au=F,yu=T,p11.est.stable)[5]-truth
  
  ######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=F,au=F,yu=T,p11.est.stable)[5]-truth
  
  corr <- cor(fitted(fitc.mar),fitted(glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=F,au=F,yu=T,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.c <- replicate(B,do.oneZ.c(n))
result.Z.ca <- mclapply(trials,do.oneZ.ca,mc.cores=numCores)
result2.Z.ca <- matrix(unlist(result.Z.ca),ncol=B)
resZ.ca <- c(apply(result2.Z.ca,1,mean)*100)
resZ.ca
resZ.ca.sd <- apply(result2.Z.ca,1,sd)*100
beep()

do.oneZ.cay<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  
  ##########
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=F,au=F,yu=F,p11.est.stable)[5]-truth
  
  ######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=F,au=F,yu=F,p11.est.stable)[5]-truth
  
  corr <- cor(fitted(fitc.mar),fitted(glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=F,au=F,yu=F,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.c <- replicate(B,do.oneZ.c(n))
result.Z.cay <- mclapply(trials,do.oneZ.cay,mc.cores=numCores)
result2.Z.cay <- matrix(unlist(result.Z.cay),ncol=B)
resZ.cay <- c(apply(result2.Z.cay,1,mean)*100)
resZ.cay
resZ.cay.sd <- apply(result2.Z.cay,1,sd)*100
beep()

do.oneZ.cy<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  
  ##########
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=F,au=T,yu=F,p11.est.stable)[5]-truth
  
  ######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=F,au=T,yu=F,p11.est.stable)[5]-truth
  
  corr <- cor(fitted(fitc.mar),fitted(glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=F,au=T,yu=F,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.c <- replicate(B,do.oneZ.c(n))
result.Z.cy <- mclapply(trials,do.oneZ.cy,mc.cores=numCores)
result2.Z.cy <- matrix(unlist(result.Z.cy),ncol=B)
resZ.cy <- c(apply(result2.Z.cy,1,mean)*100)
resZ.cy
resZ.cy.sd <- apply(result2.Z.cy,1,sd)*100
beep()

do.oneZ.a<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  #########
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                        mu=T,cu=T,au=F,yu=T,p11.est.stable)[5]-truth
  ########
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=T,au=F,yu=T,p11.est.stable)[5]-truth
  #print(summary(fita))
  corr <- cor(fitted(fita),fitted(glm(A~U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=T,au=F,yu=T,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.a <- replicate(B,do.oneZ.a(n))
result.Z.a <- mclapply(trials,do.oneZ.a,mc.cores=numCores)
result2.Z.a <- matrix(unlist(result.Z.a),ncol=B)
resZ.a <- c(apply(result2.Z.a,1,mean)*100)
resZ.a
resZ.a.sd <- apply(result2.Z.a,1,sd)*100
beep()

do.oneZ.ay<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  
  #########
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=T,au=F,yu=F,p11.est.stable)[5]-truth
  ########
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- dataset.Z$A*newpa^(-1) + (1-dataset.Z$A)*(1-newpa)^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=T,au=F,yu=F,p11.est.stable)[5]-truth
  #print(summary(fita))
  corr <- cor(fitted(fita),fitted(glm(A~U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=T,au=F,yu=F,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.a <- replicate(B,do.oneZ.a(n))
result.Z.ay <- mclapply(trials,do.oneZ.ay,mc.cores=numCores)
result2.Z.ay <- matrix(unlist(result.Z.ay),ncol=B)
resZ.ay <- c(apply(result2.Z.ay,1,mean)*100)
resZ.ay
resZ.ay.sd <- apply(result2.Z.ay,1,sd)*100
beep()

do.oneZ.y<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  
  ######
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                        mu=T,cu=T,au=T,yu=F,p11.est.stable)[5]-truth
  
  #######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- (dataset.Z$A*newpa + (1-dataset.Z$A)*(1-newpa))^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$U1,dataset.Z$U2,dataset.Z$U3,dataset.Z$U4)
  colnames(datam0) <- c("A","U1","U2","U3","U4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=T,cu=T,au=T,yu=F,p11.est.stable)[5]-truth
  
  #print(summary(fita))
  corr <- cor(fitted(fity),fitted(glm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=T,au=T,yu=F,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)
}
#result.Z.y <- replicate(B,do.oneZ.y(n))
result.Z.y <- mclapply(trials,do.oneZ.y,mc.cores=numCores)
result2.Z.y <- matrix(unlist(result.Z.y),ncol=B)
resZ.y <- c(apply(result2.Z.y,1,mean)*100)
resZ.y
resZ.y.sd <- apply(result2.Z.y,1,sd)*100
beep()

do.oneZ.wrong<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  
  ##weight
  
  ########
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable1 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                        mu=F,cu=F,au=F,yu=F,p11.est.stable)[5]-truth
  #######
  hatc1 <- -log(1-mean(dataset.Z$A))+log(mean(dataset.Z$A*(1-fitted(fita))/fitted(fita)))
  newpa <- expit(logit(fitted(fita))+hatc1)
  weighta <- (dataset.Z$A*newpa + (1-dataset.Z$A)*(1-newpa))^(-1)
  
  fitc.mar.stable <- glm(cbind(C,1-C)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  fitm.mar.stable <- glm(cbind(M,1-M)~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z,
                         weights = weighta)
  
  fc.stable <- fitted(fitc.mar.stable)
  fm.stable <- fitted(fitm.mar.stable)
  datam0 <- data.frame(0,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4)
  colnames(datam0) <- c("A","Z1","Z2","Z3","Z4")
  pm0 <- expit(predict(fitm.mar.stable,newdata=datam0))
  weightm <- dataset.Z$M*pm0+(1-dataset.Z$M)*(1-pm0)
  weightc <- dataset.Z$C*fc.stable+(1-dataset.Z$C)*(1-fc.stable)
  p11.est.stable <- estp11(fitm.mar.stable,fitc.mar.stable,dataset.Z$M,dataset.Z$C,n)
  weightmc <- dataset.Z$M*dataset.Z$C*p11.est.stable+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-fc.stable-fm.stable+p11.est.stable)+
    (1-dataset.Z$M)*dataset.Z$C*(fc.stable-p11.est.stable)+(1-dataset.Z$C)*dataset.Z$M*(fm.stable-p11.est.stable)
  weight1 <- ((fitted(fita))^(-1))*weightm*weightc/weightmc
  
  fity.stable <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z,weights=weight1)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  
  quad.stable2 <- quadZu(fitm.mar.stable,fitc.mar.stable,fity.stable,fita,dataset = dataset.Z,
                         mu=F,cu=F,au=F,yu=F,p11.est.stable)[5]-truth
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=F,au=F,yu=F,p11.est)-truth
  # print(res)
  c(res,quad.stable1,quad.stable2)

}
#result.Z.wrong <- replicate(B,do.oneZ.wrong(n))
result.Z.wrong <- mclapply(trials,do.oneZ.wrong,mc.cores=numCores)
result2.Z.wrong <- matrix(unlist(result.Z.wrong),ncol=B)
resZ.wrong <- apply(result2.Z.wrong,1,mean)*100
resZ.wrong
resZ.wrong.sd <- apply(result2.Z.wrong,1,sd)*100
resZ.wrong.sd
beep()
##############################


res.final <- rbind(resZ.corr,resZ.a,resZ.c,resZ.m,resZ.y,
                   resZ.ma,resZ.mc,resZ.my,resZ.ca,resZ.cy,resZ.ay,
                   resZ.mcy,resZ.mca,resZ.cay,resZ.may,
                   resZ.wrong)
res.final.sd <- rbind(resZ.corr.sd,resZ.a.sd,resZ.c.sd,resZ.m.sd,resZ.y.sd,
                      resZ.ma.sd,resZ.mc.sd,resZ.my.sd,resZ.ca.sd,resZ.cy.sd,resZ.ay.sd,
                      resZ.mcy.sd,resZ.mca.sd,resZ.cay.sd,resZ.may.sd,
                      resZ.wrong.sd)
colnames(res.final) <- c("delta1","delta2","delta3","delta4","delta_quad","stable.quad1","stable.quad2")
# rownames(res.final) <- c("All correct","incorrect $A$", "incorrect $C$", "incorrect $M$",
#                          "incorrect $EY$", "incorrect $M$, $A$","incorrect $M$ , $C$",
#                          "incorrect $M$ , $EY$","incorrect $C$ , $A$","incorrect $C$ , $EY$",
#                          "incorrect $A$ , $EY$", "incorrect $M$, $C$, EY$","incorrect $M$, $C$, $A$",
#                          "incorrect $C$, $A$, $EY$","incorrect $M$, $A$, $EY$",
#                          "incorrect all wrong")
res.final <- data.frame(res.final)
res.final$type <- 1
colnames(res.final.sd) <- c("delta1","delta2","delta3","delta4","delta_quad","stable.quad1","stable.quad2")
# rownames(res.final.sd) <- c("All correct","incorrect $A$", "incorrect $C$", "incorrect $M$",
#                          "incorrect $EY$", "incorrect $M$, $A$","incorrect $M$ , $C$",
#                          "incorrect $M$ , $EY$","incorrect $C$ , $A$","incorrect $C$ , $EY$",
#                          "incorrect $A$ , $EY$", "incorrect $M$, $C$, EY$","incorrect $M$, $C$, $A$",
#                          "incorrect $C$, $A$, $EY$","incorrect $M$, $A$, $EY$",
#                          "incorrect all wrong")
res.final.sd <- data.frame(res.final.sd)
res.final.sd$type <- 2
write.csv(res.final,"stable1tchetgen1000.csv")
write.csv(res.final.sd,"stable1sdtchetgen1000.csv")

res.final <- read.csv("stable1tchetgen500.csv",header = T)
res.final.sd <- read.csv("stable1sdtchetgen500.csv",header = T)
library("xtable")

colnames(res.final) <- c("scenario","delta1","delta2","delta3","delta4","delta_quad","stable.quad1","stable.quad2")
# rownames(res.final) <- c("All correct","incorrect $A$", "incorrect $C$", "incorrect $M$",
#                          "incorrect $EY$", "incorrect $M$, $A$","incorrect $M$ , $C$",
#                          "incorrect $M$ , $EY$","incorrect $C$ , $A$","incorrect $C$ , $EY$",
#                          "incorrect $A$ , $EY$", "incorrect $M$, $C$, EY$","incorrect $M$, $C$, $A$",
#                          "incorrect $C$, $A$, $EY$","incorrect $M$, $A$, $EY$",
#                          "incorrect all wrong")
res.final <- data.frame(res.final)
res.final$type <- 1
colnames(res.final.sd) <- c("scenario","delta1","delta2","delta3","delta4","delta_quad","stable.quad1","stable.quad2")
# rownames(res.final.sd) <- c("All correct","incorrect $A$", "incorrect $C$", "incorrect $M$",
#                          "incorrect $EY$", "incorrect $M$, $A$","incorrect $M$ , $C$",
#                          "incorrect $M$ , $EY$","incorrect $C$ , $A$","incorrect $C$ , $EY$",
#                          "incorrect $A$ , $EY$", "incorrect $M$, $C$, EY$","incorrect $M$, $C$, $A$",
#                          "incorrect $C$, $A$, $EY$","incorrect $M$, $A$, $EY$",
#                          "incorrect all wrong")
res.final.sd <- data.frame(res.final.sd)
res.final.sd$type <- 2

table.out <- res.final[,-8]
table.out.sd <- res.final.sd[,-8]

index <- rep(NA,32)
for(i in 1:32)
{
  if(i %% 2 ==0)
  {
    index[i] <- i/2+16
  }
  if(i %% 2 !=0)
  {
    index[i] <- (i+1)/2
  }
}
  
try <- rbind(table.out,table.out.sd)
new_df <- try[index,]
finalout <- new_df[,c(1,8,2:7)]

xtable(finalout,digits=1)
