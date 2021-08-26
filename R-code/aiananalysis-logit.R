library(BB)
usedata <- read.csv("dataset1.csv",header = T)
dataset1 <- usedata
N <- nrow(dataset1)
dataset <- dataset1
#datainad <- dataset[which(dataset$A=="inadequate" || dataset$A=="adequate"),]
fitm.mar <- glm(M~A+mager+meduc+dmar,data=dataset,family = binomial(link="logit"))
fitc.mar <- glm(C~A+mager+meduc+dmar,family = binomial(link="logit"),data = dataset)
fita <- glm(A~mager+meduc+dmar,family = binomial(link = "logit"),data = dataset)
fity <- glm(Y~A+mager+meduc+dmar+C+M,family = binomial(link = "logit"),data = dataset)
V <- cbind(dataset$mager,dataset$meduc,dataset$dmar)
X <- V
A <- dataset$A
Y <- dataset$Y
C <- dataset$C
M <- dataset$M
fc <- fitted(fitc.mar)
fm <- fitted(fitm.mar)
n <- nrow(dataset)
esttau <- read.csv("parfinal.csv",header = T)[,2]
phiest <- exp(esttau[1]+esttau[2]*A+X%*%esttau[-c(1,2)])
p11.est <- rep(NA,n)
p11.fun <- function(phi,pc,pm)
{
  if (phi==1)
  {
    return(pc*pm)
  }
  if(phi!=1)
  {
    s <- sqrt((1+(pc+pm)*(phi-1))^2+4*phi*(1-phi)*pc*pm)
    numer <- 1+(pm+pc)*(phi-1)-s
    dnom <- 2*(phi-1)
    return(numer/dnom)
  }
}
for(i in 1:n)
{
  p11.est[i] <- p11.fun(phiest[i],fc[i],fm[i])
}
#help functions:
logit <- function(x){log(x/(1 - x))}
expit <- function(x){exp(x)/(1 + exp(x))}



Fmc <- M*C*p11.est+(1-M)*(1-C)*(1-fc-fm+p11.est)+(1-M)*C*(fc-p11.est)+(1-C)*M*(fm-p11.est)
dc.mar <- dbinom(C,1,prob=fitted(fitc.mar))
dm.mar <- dbinom(M,1,prob=fitted(fitm.mar))


#f(M|A=0,X)
datam0x <- data.frame(cbind(0,dataset$mager,dataset$meduc,dataset$dmar))
colnames(datam0x) <- c("A","mager","meduc","dmar")
pM0x <- expit(predict(fitm.mar,newdata = datam0x))
dm0x <- dbinom(M,1,pM0x)

##A and Y
proba <- dbinom(A,1,fitted(fita))
Ey <- fitted(fity)


#estimate
#delta1*******see if it still has this form?I think so
ratio1 <- ((2*A-1)/proba)*(dm0x*dc.mar/Fmc)
est_delta1 <- ratio1*Y
est_comp1 <- ratio1*(Y-Ey)

#delta2
#eta
a <- 1
etafun <- function(a)
{
  datam <- data.frame(cbind(0,dataset$mager,dataset$meduc,dataset$dmar))
  colnames(datam) <- c("A","mager","meduc","dmar")
  pM <- expit(predict(fitm.mar,newdata = datam))
  data1 <- data.frame(cbind(a,C,1,dataset$mager,dataset$meduc,dataset$dmar))
  colnames(data1) <- c("A","C","M","mager","meduc","dmar")
  data0 <- data.frame(cbind(a,C,0,dataset$mager,dataset$meduc,dataset$dmar))
  colnames(data0) <- c("A","C","M","mager","meduc","dmar")
  #
  eta <- expit(predict(fity,newdata = data1))*pM+expit(predict(fity,newdata = data0))*(1-pM)
  #
  
  eta
}
ratio2 <- (2*A-1)/proba
est_delta2 <- est_comp2 <- ratio2*etafun(A)

#delta3
gammafun <-function(a)
{
  datac <- data.frame(cbind(a,dataset$mager,dataset$meduc,dataset$dmar))
  colnames(datac) <- c("A","mager","meduc","dmar")
  pC <- expit(predict(fitc.mar,newdata = datac))
  data <- data.frame(cbind(a,pC,M,dataset$mager,dataset$meduc,dataset$dmar))
  colnames(data) <- c("A","C","M","mager","meduc","dmar")
  gamma <- expit(predict(fity,newdata = data))
  gamma
}
ratio3 <- (1-A)/proba
est_delta3 <- est_comp3 <- ratio3*(gammafun(1)-gammafun(0))

#delta4
taufun <- function(a)
{
  datac <- data.frame(cbind(a,dataset$mager,dataset$meduc,dataset$dmar))
  colnames(datac) <- c("A","mager","meduc","dmar")
  pC <- expit(predict(fitc.mar,newdata=datac))
  datam <- data.frame(cbind(0,dataset$mager,dataset$meduc,dataset$dmar))
  colnames(datam) <- c("A","mager","meduc","dmar")
  pM <- expit(predict(fitm.mar,newdata=datam))
  data <- data.frame(cbind(a,pC,pM,dataset$mager,dataset$meduc,dataset$dmar))
  colnames(data) <- c("A","C","M","mager","meduc","dmar")
  tau <- expit(predict(fity,newdata = data))
  tau
}
est_delta4 <- taufun(1)-taufun(0)
est_comp4 <- (1-ratio3)*est_delta4

est_comp5 <- -ratio2*taufun(A)

quad <- mean(est_comp1 + est_comp2 + est_comp3 + est_comp4 + est_comp5)

result <- data.frame(mean(est_delta1),mean(est_delta2),mean(est_delta3),mean(est_delta4),mean(quad))

result
write.csv(result,"pointestimate-logit.csv")
pi <- fitted(fita)
ymar.fun <- function(a)
{
  phiest.temp <- exp(esttau[1]+esttau[2]*a+X%*%esttau[-c(1,2)])
  p11.est.temp <- rep(NA,n)
  for(i in 1:n)
  {
    p11.est.temp[i] <- p11.fun(phiest.temp[i],fc[i],fm[i])
  }
  
  data11 <- cbind(a,dataset$mager,dataset$meduc,dataset$dmar,1,1)
  colnames(data11) <- c("A","mager","meduc","dmar","M","C")
  data11 <- data.frame(data11)
  ymar.fitted11 <- expit(predict(fity,newdata = data11))
  data10 <- cbind(a,dataset$mager,dataset$meduc,dataset$dmar,1,0)
  colnames(data10) <- c("A","mager","meduc","dmar","M","C")
  data10 <- data.frame(data10)
  ymar.fitted10 <- expit(predict(fity,newdata = data10))
  data01 <- cbind(a,dataset$mager,dataset$meduc,dataset$dmar,0,1)
  colnames(data01) <- c("A","mager","meduc","dmar","M","C")
  data01 <- data.frame(data01)
  ymar.fitted01 <- expit(predict(fity,newdata = data01))
  data00 <- cbind(a,dataset$mager,dataset$meduc,dataset$dmar,0,0)
  colnames(data00) <- c("A","mager","meduc","dmar","M","C")
  data00 <- data.frame(data00)
  ymar.fitted00 <- expit(predict(fity,newdata = data00))
  
  ymar.fitted <- ymar.fitted11*p11.est.temp+ymar.fitted00*(1-fc-fm+p11.est.temp)+
    ymar.fitted01*(fc-p11.est.temp)+ymar.fitted10*(fm-p11.est.temp)
  ymar.fitted
}
ateest <- mean(A*Y/pi-(A-pi)*ymar.fun(1)/pi-(1-A)*Y/(1-pi)-(A-pi)*ymar.fun(0)/(1-pi))

write.csv(ateest,"ateest-logit.csv")
ateest <- read.csv("ateest-logit.csv")[1,2]
quadest <- read.csv("pointestimate-logit.csv")[1,6]
ateest-quadest
