stable<-function(dataset.Z)
{
  fita <- glm(A~Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  weighta <- (dataset.Z$A*fitted(fita)+(1-dataset.Z$A)*(1-fitted(fita)))^(-1)
  
  ##weight
  w.fitc0 <- glm(cbind(C,1-C)~Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z[which(dataset.Z$A==0),],
               weights = weighta[which(dataset.Z$A==0)])
  
  pc0.dagger <- expit(predict(w.fitc0,newdata = dataset.Z)) 
  fc0.dagger <- dataset.Z$C*expit(predict(w.fitc0,newdata = dataset.Z))+
                (1-dataset.Z$C)*(1-expit(predict(w.fitc0,newdata = dataset.Z)))
  
  w.fitc1 <- glm(cbind(C,1-C)~Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z[which(dataset.Z$A==1),],
                 weights = weighta[which(dataset.Z$A==1)])
  pc1.dagger <- expit(predict(w.fitc1,newdata = dataset.Z))
  fc1.dagger <- dataset.Z$C*expit(predict(w.fitc1,newdata = dataset.Z))+
    (1-dataset.Z$C)*(1-expit(predict(w.fitc1,newdata = dataset.Z)))
  
  w.fitm0 <- glm(cbind(M,1-M)~Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z[which(dataset.Z$A==0),],
                 weights = weighta[which(dataset.Z$A==0)])
  pm0.dagger <- expit(predict(w.fitm0,newdata = dataset.Z))
  fm0.dagger <- dataset.Z$M*expit(predict(w.fitm0,newdata = dataset.Z))+
    (1-dataset.Z$M)*(1-expit(predict(w.fitm0,newdata = dataset.Z)))
  
  w.fitm1 <- glm(cbind(M,1-M)~Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z[which(dataset.Z$A==1),],
                 weights = weighta[which(dataset.Z$A==1)])
  pm1.dagger <- expit(predict(w.fitm1,newdata = dataset.Z))
  fm1.dagger <- dataset.Z$M*expit(predict(w.fitm1,newdata = dataset.Z))+
    (1-dataset.Z$M)*(1-expit(predict(w.fitm1,newdata = dataset.Z)))
  
  p11.est.stable1 <- estp11.w(w.fitm1,w.fitc1,dataset.Z$M,dataset.Z$C,n,dataset.Z)
  p11.est.stable0 <- estp11.w(w.fitm0,w.fitc0,dataset.Z$M,dataset.Z$C,n,dataset.Z)
  
  fmc1 <- dataset.Z$M*dataset.Z$C*p11.est.stable1+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-pc1.dagger-pm1.dagger+p11.est.stable1)+
    (1-dataset.Z$M)*dataset.Z$C*(pc1.dagger-p11.est.stable1)+(1-dataset.Z$C)*dataset.Z$M*(pm1.dagger-p11.est.stable1) 
  fmc0 <- dataset.Z$M*dataset.Z$C*p11.est.stable0+(1-dataset.Z$M)*(1-dataset.Z$C)*(1-pc0.dagger-pm0.dagger+p11.est.stable0)+
    (1-dataset.Z$M)*dataset.Z$C*(pc0.dagger-p11.est.stable0)+(1-dataset.Z$C)*dataset.Z$M*(pm0.dagger-p11.est.stable0) 
  
  weight1 <- (fc1.dagger*fm0.dagger/fmc1)*(fitted(fita))^(-1)
  weight0 <- (fc0.dagger*fm0.dagger/fmc0)*(1-fitted(fita))^(-1)
  
  w.fity1 <- lm(Y~Z1+Z2+Z3+Z4+C+M,data = dataset.Z[which(dataset.Z$A==1),],weights=weight1[which(dataset.Z$A==1)])
  w.fity0 <- lm(Y~Z1+Z2+Z3+Z4+C+M,data = dataset.Z[which(dataset.Z$A==0),],weights=weight1[which(dataset.Z$A==0)])

  Ey1 <- predict(w.fity1,newdata = dataset.Z)
  Ey0 <- predict(w.fity0,newdata = dataset.Z)
  
  eta <- function(a)
  {
    if(a==1)
    {
      return((cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,dataset.Z$C,1) %*% coef(w.fity1)) * pm0.dagger +
        (cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,dataset.Z$C,0) %*% coef(w.fity1)) * (1-pm0.dagger))
    }
    if(a==0)
    {
      return((cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,dataset.Z$C,1) %*% coef(w.fity0)) * pm0.dagger +
        (cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,dataset.Z$C,0) %*% coef(w.fity0)) * (1-pm0.dagger))
    }
  }
  
  gamma <- function(a)
  {
    if(a==1)
    {
      return((cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,1,dataset.Z$M) %*% coef(w.fity1)) * pc1.dagger +
               (cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,0,dataset.Z$M) %*% coef(w.fity1)) * (1-pc1.dagger))
    }
    if(a==0)
    {
      return((cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,1,dataset.Z$M) %*% coef(w.fity0)) * pc0.dagger +
               (cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,0,dataset.Z$M) %*% coef(w.fity0)) * (1-pc0.dagger))
    }
  }
  
  tau <- function(a)
  {
    if(a==1)
    {
      return((cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,1,1) %*% coef(w.fity1)) * pc1.dagger * pm0.dagger +
               (cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,1,0) %*% coef(w.fity1)) * pc1.dagger * (1-pm0.dagger)+
               (cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,0,1) %*% coef(w.fity1)) * (1-pc1.dagger) * pm0.dagger +
        (cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,0,0) %*% coef(w.fity1)) * (1-pc1.dagger) * (1-pm0.dagger) )
        
    }
    if(a==0)
    {
      return((cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,1,1) %*% coef(w.fity0)) * pc0.dagger * pm0.dagger +
               (cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,1,0) %*% coef(w.fity0)) * pc0.dagger * (1-pm0.dagger)+
             (cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,0,1) %*% coef(w.fity0)) * (1-pc0.dagger) * pm0.dagger +
        (cbind(1,dataset.Z$Z1,dataset.Z$Z2,dataset.Z$Z3,dataset.Z$Z4,0,0) %*% coef(w.fity0)) * (1-pc0.dagger) * (1-pm0.dagger) )
    }
  }

  mu1 <- ((dataset.Z$A/fitted(fita))*fm0.dagger*fc1.dagger/fmc1)*(dataset.Z$Y - Ey1) +
         (dataset.Z$A/fitted(fita))*(eta(1)-tau(1)) + ((1-dataset.Z$A)/(1-fitted(fita)))*(gamma(1)-tau(1)) + tau(1)
  mu0 <- (((1-dataset.Z$A)/(1-fitted(fita)))*fm0.dagger*fc0.dagger/fmc0)*(dataset.Z$Y - Ey0) +
  ((1-dataset.Z$A)/(1-fitted(fita)))*(eta(0)-tau(0)) + ((1-dataset.Z$A)/(1-fitted(fita)))*(gamma(0)-tau(0)) + tau(0)
  theta <- mean(mu1 - mu0)
  
  return(theta)
}

do.oneZ.wrong.w<-function(trials)
{
  dataset.Z <- data_genZ(n)
  #
  fitm.mar <- glm(M~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data=dataset.Z)
  fitc.mar <- glm(C~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+C+M,data = dataset.Z)
  fity.tru <- lm(Y~A+U1+U2+U3+U4+C+M,data = dataset.Z)
  
  
  quad.stable1 <- stable(dataset.Z) - truth
  
  #print(summary(fita))
  p11.est <- estp11.w(fitm.mar,fitc.mar,dataset.Z$M,dataset.Z$C,n,dataset.Z)
  res <- quadZu(fitm.mar,fitc.mar,fity,fita,dataset = dataset.Z,
                mu=F,cu=F,au=F,yu=F,p11.est)-truth
  # print(res)
  c(res,quad.stable1)
  
}

do.oneZ.wrong(2)

##############################
#result.Z.wrong <- replicate(B,do.oneZ.wrong(n))
n <- 500
B <- 200
trials <- seq(1:B)
numCores <- detectCores()
set.seed(201)
result.Z.wrong <- mclapply(trials,do.oneZ.wrong,mc.cores=numCores)
result2.Z.wrong <- matrix(unlist(result.Z.wrong),ncol=B)
resZ.wrong <- apply(result2.Z.wrong,1,mean)*100
resZ.wrong
resZ.wrong.sd <- apply(result2.Z.wrong,1,sd)*100
resZ.wrong.sd
beep()
