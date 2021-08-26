set.seed(9560)
library(beepr)
######do not observe U
do.oneZ.corr <- function(n)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  #print(summary(fita))
  
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=T,au=T,yu=T,p11.est)-truth
  # print(res)
  c(res,min(fitted(fita)),NA)
}
B <- 1000
n <- 1000
result.Z.corr <- replicate(B,do.oneZ.corr(n))
result2.Z.corr <- matrix(unlist(result.Z.corr),ncol=B)
resZ.corr <- c(apply(result2.Z.corr,1,mean)[1:5]*100, apply(result2.Z.corr,1,min)[6:7])
resZ.corr
resZ.corr.sd <- apply(result2.Z.corr,1,var)[1:5]*100
beep()

do.oneZ.m<-function(n)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  corr <- cor(fitted(fitm.mar),fitted(glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=T,au=T,yu=T,p11.est)-truth
  # print(res)
  c(res,min(fitted(fita)),corr)
}
B <- 1500
n <- 1000
result.Z.m <- replicate(B,do.oneZ.m(n))
result2.Z.m <- matrix(unlist(result.Z.m),ncol=B)
resZ.m <- c(apply(result2.Z.m,1,mean)[1:5]*100, apply(result2.Z.m,1,min)[6:7])
resZ.m

resZ.m.sd <- apply(result2.Z.m,1,sd)[1:5]*100

do.oneZ.c<-function(n)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  corr <- cor(fitted(fitc.mar),fitted(glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=F,au=T,yu=T,p11.est)-truth
  # print(res)
  c(res,min(fitted(fita)),corr)
}
B <- 1500
n <- 1000
result.Z.c <- replicate(B,do.oneZ.c(n))
result2.Z.c <- matrix(unlist(result.Z.c),ncol=B)
resZ.c <- c(apply(result2.Z.c,1,mean)[1:5]*100, apply(result2.Z.c,1,min)[6:7])
resZ.c

resZ.c.sd <- apply(result2.Z.c,1,sd)[1:5]*100

do.oneZ.a<-function(n)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  #print(summary(fita))
  corr <- cor(fitted(fita),fitted(glm(A~U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=T,au=F,yu=T,p11.est)-truth
  # print(res)
  c(res,min(fitted(fita)),corr)
}
B <- 1500
n <- 1000
result.Z.a <- replicate(B,do.oneZ.a(n))
result2.Z.a <- matrix(unlist(result.Z.a),ncol=B)
resZ.a <- c(apply(result2.Z.a,1,mean)[1:5]*100, apply(result2.Z.a,1,min)[6:7])
resZ.a

resZ.a.sd <- apply(result2.Z.a,1,sd)[1:5]*100

do.oneZ.y<-function(n)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+U1+U2+U3+U4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  
  #print(summary(fita))
  corr <- cor(fitted(fity),fitted(glm(Y~A+U1+U2+U3+U4+C+M,family = binomial(link = "logit"),data = dataset.Z)))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=T,cu=T,au=T,yu=F,p11.est)-truth
  # print(res)
  c(res,min(fitted(fita)),corr)
}
B <- 1500
n <- 1000
result.Z.y <- replicate(B,do.oneZ.y(n))
result2.Z.y <- matrix(unlist(result.Z.y),ncol=B)
resZ.y <- c(apply(result2.Z.y,1,mean)[1:5]*100, apply(result2.Z.y,1,min)[6:7])
resZ.y

resZ.y.sd <- apply(result2.Z.y,1,sd)[1:5]*100

do.oneZ.wrong<-function(n)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  
  #print(summary(fita))
  p11.est <- estp11(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=F,au=F,yu=F,p11.est)-truth
  # print(res)
  c(res,min(fitted(fita)),NA)
}
B <- 1500
n <- 1000
result.Z.wrong <- replicate(B,do.oneZ.wrong(n))
result2.Z.wrong <- matrix(unlist(result.Z.wrong),ncol=B)
resZ.wrong <- c(apply(result2.Z.wrong,1,mean)[1:5]*100, apply(result2.Z.wrong,1,min)[6:7])
resZ.wrong

resZ.wrong.sd <- apply(result2.Z.wrong,1,sd)[1:5]*100

res.final <- rbind(resZ.corr,resZ.a,resZ.c,resZ.m,resZ.y,resZ.wrong)
res.final.sd <- rbind(resZ.corr.sd,resZ.a.sd,resZ.c.sd,resZ.m.sd,resZ.y.sd,resZ.wrong.sd)
colnames(res.final) <- c("delta1","delta2","delta3","delta4","delta_quad","min(fitted(fita))","correlation")
rownames(res.final) <- c("correct","incorrect fita", "incorrect fitc", "incorrect fitm","incorrect fity","all wrong")
colnames(res.final.sd) <- c("delta1","delta2","delta3","delta4","delta_quad")
rownames(res.final.sd) <- c("correct","incorrect fita", "incorrect fitc", "incorrect fitm","incorrect fity","all wrong")

