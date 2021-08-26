quadZu <- function(fitm.mar,fitc.mar,fity,fita,dataset,p11.est)
{
  Xu <- cbind(dataset$U1, dataset$U2, dataset$U3)
  A <- dataset$A
  C <- dataset$C
  M <- dataset$M
  Y <- dataset$Y
  
  fc <- fitted(fitc.mar)
  fm <- fitted(fitm.mar)
  
  Fmc <- M*C*p11.est+(1-M)*(1-C)*(1-fc-fm+p11.est)+(1-M)*C*(fc-p11.est)+(1-C)*M*(fm-p11.est)
  dc.mar <- dbinom(C,1,prob=fitted(fitc.mar))
  dm.mar <- dbinom(M,1,prob=fitted(fitm.mar))
  
  
  #f(M|A=0,X)

    datam0x <- data.frame(cbind(0,Xu))
    colnames(datam0x) <- c("A","U1","U2","U3")
  
  pM0x <- expit(predict(fitm.mar,newdata = datam0x))
  dm0x <- dbinom(M,1,pM0x)
  
  ##A and Y
  proba <- dbinom(A,1,fitted(fita))
  Ey <- predict(fity)
  
  
  #estimate
  #delta1*******check
  ratio1 <- ((2*A-1)/proba)*(dm0x*dc.mar/Fmc)
  est_delta1 <- ratio1*Y
  est_comp1 <- ratio1*(Y-Ey)
  
  #delta2
  #eta
  #a <- 1
  etafun <- function(a)
  {
      datam <- data.frame(cbind(0,Xu))
      colnames(datam) <- c("A","U1","U2","U3")

    pM <- expit(predict(fitm.mar,newdata = datam))

      data1 <- data.frame(cbind(a,C,1,Xu))
      colnames(data1) <- c("A","C","M","U1","U2","U3")
      data0 <- data.frame(cbind(a,C,0,Xu))
      colnames(data0) <- c("A","C","M","U1","U2","U3")


    eta <- predict(fity,newdata = data1)*pM+predict(fity,newdata = data0)*(1-pM)
    eta
  }
  ratio2 <- (2*A-1)/proba
  est_delta2 <- est_comp2 <- ratio2*etafun(A)
  
  #delta3
  gammafun <-function(a)
  {

      datac <- data.frame(cbind(a,Xu))
      colnames(datac) <- c("A","U1","U2","U3")
    
    pC <- expit(predict(fitc.mar,newdata = datac))

      datac1 <- data.frame(cbind(a,1,M,Xu))
      colnames(datac1) <- c("A","C","M","U1","U2","U3")
      datac0 <- data.frame(cbind(a,0,M,Xu))
      colnames(datac0) <- c("A","C","M","U1","U2","U3")
  
    gamma <- predict(fity,newdata = datac0)*(1-pC)+predict(fity,newdata = datac1)*pC
    gamma
  }
  ratio3 <- (1-A)/proba
  est_delta3 <- est_comp3 <- ratio3*(gammafun(1)-gammafun(0))
  
  #delta4
  taufun <- function(a)
  {
      datac <- data.frame(cbind(a,Xu))
      colnames(datac) <- c("A","U1","U2","U3")
    
    pC <- expit(predict(fitc.mar,newdata=datac))

      datam <- data.frame(cbind(0,Xu))
      colnames(datam) <- c("A","U1","U2","U3")
    
    pM <- expit(predict(fitm.mar,newdata=datam))

      datac0m0 <- data.frame(cbind(a,0,0,Xu))
      colnames(datac0m0) <- c("A","C","M","U1","U2","U3")
      datac1m1 <- data.frame(cbind(a,1,1,Xu))
      colnames(datac1m1) <- c("A","C","M","U1","U2","U3")
      datac1m0 <- data.frame(cbind(a,1,0,Xu))
      colnames(datac1m0) <- c("A","C","M","U1","U2","U3")
      datac0m1 <- data.frame(cbind(a,0,1,Xu))
      colnames(datac0m1) <- c("A","C","M","U1","U2","U3")
   
    tau <- predict(fity,newdata = datac0m0)*(1-pC)*(1-pM)+
      predict(fity,newdata = datac1m1)*pC*pM+
      predict(fity,newdata = datac1m0)*pC*(1-pM)+
      predict(fity,newdata = datac0m1)*(1-pC)*pM
    tau
  }
  est_delta4 <- taufun(1)-taufun(0)
  est_comp4 <- (1-ratio3)*est_delta4
  
  est_comp5 <- -ratio2*taufun(A)
  
  quad <- mean(est_comp1 + est_comp2 + est_comp3 + est_comp4 + est_comp5)
  
  result <- data.frame(mean(est_delta1),mean(est_delta2),mean(est_delta3),mean(est_delta4),mean(quad))
  
  result
}