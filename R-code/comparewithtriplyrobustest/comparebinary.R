set.seed(336)
library(nleqslv)
library(MASS)
library(beepr)
library(pbapply)
library(parallel)
#logit and expit
logit <- function(x) log(x/(1 - x))
expit <- function(x) exp(x)/(1 + exp(x))

n <- 1500
X <- rnorm(n)
#phi
#
tau_0 <- 1
tau_A <- -2
tau_X <- 3
#phi <- exp(tau_0+tau_A*A+tau_X*X)
#A
delta_0 <- -0.2
delta_X <- 0.3
pa <- expit(delta_0 + delta_X*X)
A <- rbinom(n,1,pa)
#C
alpha_0 <- -0.2
alpha_A <- -0.1
alpha_X <- 0.3
#pc <- expit(alpha_0 + alpha_A*A + alpha_X*X)
#mean(expit(pc-alpha_A*A))
#M.mar
gamma_0 <- -0.3
gamma_A <- -0.2
gamma_X <- 0.5
#pm <- expit(gamma_0+gamma_A*A+gamma_X*X)
#mean(expit(pm-gamma_A*A+0.3))
expit(gamma_0+gamma_A+gamma_X)
####



####
#parameters for Y
beta_0 <- 1
beta_A <- 3
beta_M <- 6
beta_C <- 3
beta_X <- 6
beta_AC <- 4
beta_AM <- 2
sdy <- 4



library(MASS)

##treat the treatment-induced confounding as pretreatment confounder
triple.pre <- function(A,M,C,X,Y)
{
  dataset <- data.frame(A,M,C,X,Y)
  fita <- glm(A~X+C,data = dataset,family = binomial(link = "logit"))
  fitm <- glm(M~A+X+C, data = dataset,family = binomial(link = "logit"))
  fity <- lm(Y~A+X+C+M+A*C+A*M, data = dataset)

  
  etafun <- function(a)
  {
    datam <- data.frame(cbind(0,X,C))
    colnames(datam) <- c("A","X","C")
    pM <- expit(predict(fitm,newdata = datam))
    data1 <- data.frame(cbind(a,C,1,X))
    colnames(data1) <- c("A","C","M","X")
    data0 <- data.frame(cbind(a,C,0,X))
    colnames(data0) <- c("A","C","M","X")
    eta <- predict(fity,newdata = data1)*pM+predict(fity,newdata = data0)*(1-pM)
    eta
  }
  
  
  dm.mar <- dbinom(M,1,prob=fitted(fitm))
  
  datam0x <- data.frame(cbind(0,X,C))
  colnames(datam0x) <- c("A","X","C")
  pM0x <- expit(predict(fitm,newdata = datam0x))
  dm0x <- dbinom(M,1,pM0x)
  Ey <- predict(fity)
  ##ignoring the treatment-induced confounding
  triple <- (A*dm0x/(dm.mar*fitted(fita)))*(Y-Ey) + ((1-A)/(1-fitted(fita)))*(Ey-etafun(1)) +
    etafun(1) - 
    (A/fitted(fita))*(Y-etafun(0))+etafun(0)
  
  triple
  
}

quad <- function(fitm.mar,fitc.mar,p11.est,fity,fita,A,M,C,X,Y,n)
{
  fc <- fitted(fitc.mar)
  fm <- fitted(fitm.mar)
  
  Fmc <- M*C*p11.est+(1-M)*(1-C)*(1-fc-fm+p11.est)+(1-M)*C*(fc-p11.est)+(1-C)*M*(fm-p11.est)
  dc.mar <- dbinom(C,1,prob=fitted(fitc.mar))
  dm.mar <- dbinom(M,1,prob=fitted(fitm.mar))
  
  
  #f(M|A=0,X)
  datam0x <- data.frame(cbind(0,X))
  colnames(datam0x) <- c("A","X")
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
    datam <- data.frame(cbind(0,X))
    colnames(datam) <- c("A","X")
    pM <- expit(predict(fitm.mar,newdata = datam))
    data1 <- data.frame(cbind(a,C,1,X))
    colnames(data1) <- c("A","C","M","X")
    data0 <- data.frame(cbind(a,C,0,X))
    colnames(data0) <- c("A","C","M","X")
    eta <- predict(fity,newdata = data1)*pM+predict(fity,newdata = data0)*(1-pM)
    eta
  }
  ratio2 <- (2*A-1)/proba
  est_delta2 <- est_comp2 <- ratio2*etafun(A)
  
  ##ignoring the treatment-induced confounding
  triple <- (A*dm0x/(dm.mar*fitted(fita)))*(Y-Ey) + ((1-A)/(1-fitted(fita)))*(Ey-etafun(1)) +
    etafun(1) - (A/fitted(fita))*(Y-etafun(0))+etafun(0)
  
  ##treat C as a pretreatment confounder
  triple.aspre <- triple.pre(A,M,C,X,Y)
  

  #delta3
  gammafun <-function(a)
  {
    datac <- data.frame(cbind(a,X))
    colnames(datac) <- c("A","X")
    pC <- expit(predict(fitc.mar,newdata = datac))
    datac1 <- data.frame(cbind(a,1,M,X))
    colnames(datac1) <- c("A","C","M","X")
    datac0 <- data.frame(cbind(a,0,M,X))
    colnames(datac0) <- c("A","C","M","X")
    gamma <- predict(fity,newdata = datac0)*(1-pC)+predict(fity,newdata = datac1)*pC
    gamma
  }
  ratio3 <- (1-A)/proba
  est_delta3 <- est_comp3 <- ratio3*(gammafun(1)-gammafun(0))
  
  #delta4
  taufun <- function(a)
  {
    datac <- data.frame(cbind(a,X))
    colnames(datac) <- c("A","X")
    pC <- expit(predict(fitc.mar,newdata=datac))
    datam <- data.frame(cbind(0,X))
    colnames(datam) <- c("A","X")
    pM <- expit(predict(fitm.mar,newdata=datam))
    datac0m0 <- data.frame(cbind(a,0,0,X))
    colnames(datac0m0) <- c("A","C","M","X")
    datac1m1 <- data.frame(cbind(a,1,1,X))
    colnames(datac1m1) <- c("A","C","M","X")
    datac1m0 <- data.frame(cbind(a,1,0,X))
    colnames(datac1m0) <- c("A","C","M","X")
    datac0m1 <- data.frame(cbind(a,0,1,X))
    colnames(datac0m1) <- c("A","C","M","X")
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
  
  result <- data.frame(mean(est_delta1),mean(est_delta2),mean(est_delta3),mean(est_delta4),mean(quad),
                       mean(triple),mean(triple.aspre))
  
  result
}

gen <- function(n,caltruth=T)
{
  X <- rnorm(n)
  prob_A <- expit(delta_0+delta_X*X)
  A <- rbinom(n,1,prob = prob_A)
  pc <- expit(alpha_0 + alpha_A*A + alpha_X*X)
  pm <- expit(gamma_0 + gamma_A*A + gamma_X*X)
  phi <- exp(tau_0+tau_A*A+tau_X*X)
  
  ##generate M and C
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
  
  p11.vector <- rep(NA,n)
  for(i in 1:n)
  {
    p11.vector[i] <- p11.fun(phi[i],pc[i],pm[i])
  }
  
  p10.vector <- rep(NA,n)
  for(i in 1:n)
  {
    p10.vector[i] <- pc[i]-p11.vector[i]
  }
  p01.vector <- rep(NA,n)
  for(i in 1:n)
  {
    p01.vector[i] <- pm[i]-p11.vector[i]
  }
  p00.vector <- rep(NA,n)
  for(i in 1:n)
  {
    p00.vector[i] <- 1-pm[i]-pc[i]+p11.vector[i]
  }
  
  joint <- matrix(NA, nrow = n, ncol=2)
  
  for(i in 1:n)
  {
    prob <- c(p11.vector[i],p10.vector[i],p01.vector[i],p00.vector[i])
    whichbox <- rmultinom(1,1,prob)
    index <- which(whichbox==1)
    if(index==1)
    {
      joint[i,]= c(1,1)
    }
    if(index==2)
    {
      joint[i,]= c(1,0)
    }
    if(index==3)
    {
      joint[i,]= c(0,1)
    }
    if(index==4)
    {
      joint[i,]= c(0,0)
    }
  }
  
  C <- joint[,1]
  M <- joint[,2]
  
  
  mu_y <- beta_0 + beta_A*A + beta_M*M + beta_C*C + beta_X*X + beta_AC * A*C + beta_AM*A*M
  Y <- rnorm(n,mu_y,sdy)
  
  ###truth
  M0X <- expit(gamma_0 + gamma_X*X)
  C1X <- expit(alpha_0 + alpha_A + alpha_X*X)
  C0X <- expit(alpha_0 + alpha_X*X)
  EY1cmx <- beta_0 + beta_A + beta_M*M0X + beta_C*C1X + beta_X*X + beta_AC * C1X +beta_AM * M0X
  EY0cmx <- beta_0 + beta_M*M0X + beta_C*C0X + beta_X*X
  truth <- mean(EY1cmx) - mean(EY0cmx)
  ##
  if(caltruth)
  {
    return(truth)
  }
  if(!caltruth)
  {
    return(data.frame(Y,A,X,M,C))
  }
}

estp11 <- function(fitm.mar,fitc.mar,A,M,C,X,Y,n)
{
  #psuedo MLE
  
  fc <- fitted(fitc.mar)
  fm <- fitted(fitm.mar)
  a <- b <- c <- d <- rep(NA, n)
  
  for(i in 1:n)
  {
    if (M[i]==1 & C[i]==1)
    {
      a[i] <- 1
      b[i] <- c[i] <- d[i] <- 0
    }
    if (M[i]==0 & C[i]==1)
    {
      b[i] <- 1
      a[i] <- c[i] <- d[i] <- 0
    }
    if (M[i]==1 & C[i]==0)
    {
      c[i] <- 1
      a[i] <- b[i] <- d[i] <- 0
    }
    if (M[i]==0 & C[i]==0)
    {
      d[i] <- 1
      a[i] <- b[i] <- c[i] <- 0
    }
  }
  
  fn <- function(par)
  {
    
    f1 <- function(par)
    {
      tau0 <- par[1]
      tauA <- par[2]
      tauX <- par[3]
      
      thing <- rep(NA,n)
      
      for(i in 1:n)
      {
        #phi
        phi <- exp(tau0+tauA*A+tauX*X)
        #s  
        s <- sqrt((1+(fm[i]+fc[i])*(phi[i]-1))^2+4*(1-phi[i])*phi[i]*fm[i]*fc[i])
        
        if(phi[i]==1)
        {
          p11 <- fm[i]*fc[i]
          partialphi <- 0
        }
        if(phi[i]!=1)
        {
          p11 <- (1+(fm[i]+fc[i])*(phi[i]-1)-s)/(2*(phi[i]-1))
          term11 <- fm[i]+fc[i]
          term12 <- 0.5*((1+(fm[i]+fc[i])*(phi[i]-1))^2+4*phi[i]*(1-phi[i])*fm[i]*fc[i])^(-0.5)
          term13 <- 2*(1+(fc[i]+fm[i])*(phi[i]-1))*(fm[i]+fc[i])+4*fc[i]*fm[i]-8*fc[i]*fm[i]*phi[i]
          term1 <- (term11-term12*term13)/(2*(phi[i]-1))
          term2 <- p11/(phi[i]-1)
          partialphi <- term1-term2
        }
        
        partialp11tau0 <- partialphi*phi[i]
        
        t <- a[i]/p11-b[i]/(fc[i]-p11)-c[i]/(fm[i]-p11)+d[i]/(1-fc[i]-fm[i]+p11)
        thing[i] <- t * partialp11tau0
      }
      
      sum(thing)
    }
    
    f2 <- function(par)
    {
      tau0 <- par[1]
      tauA <- par[2]
      tauX <- par[3]
      thing <- rep(NA,n)
      
      for(i in 1:n)
      {
        #phi
        phi <- exp(tau0+tauA*A+tauX*X)
        #s  
        s <- sqrt((1+(fm[i]+fc[i])*(phi[i]-1))^2+4*(1-phi[i])*phi[i]*fm[i]*fc[i])
        
        if(phi[i]==1)
        {
          p11 <- fm[i]*fc[i]
          partialphi <- 0
        }
        if(phi[i]!=1)
        {
          p11 <- (1+(fm[i]+fc[i])*(phi[i]-1)-s)/(2*(phi[i]-1))
          term11 <- fm[i]+fc[i]
          term12 <- 0.5*((1+(fm[i]+fc[i])*(phi[i]-1))^2+4*phi[i]*(1-phi[i])*fm[i]*fc[i])^(-0.5)
          term13 <- 2*(1+(fc[i]+fm[i])*(phi[i]-1))*(fm[i]+fc[i])+4*fc[i]*fm[i]-8*fc[i]*fm[i]*phi[i]
          term1 <- (term11-term12*term13)/(2*(phi[i]-1))
          term2 <- p11/(phi[i]-1)
          partialphi <- term1-term2
        }
        
        partialp11tauA <- partialphi*phi[i]*A[i]
        
        t <- a[i]/p11-b[i]/(fc[i]-p11)-c[i]/(fm[i]-p11)+d[i]/(1-fc[i]-fm[i]+p11)
        thing[i] <- t * partialp11tauA
      }
      sum(thing)
      
    }
    
    f3 <- function(par)
    {
      tau0 <- par[1]
      tauA <- par[2]
      tauX <- par[3]
      #phi
      thing <- rep(NA,n)
      for(i in 1:n)
      {
        phi <- exp(tau0+tauA*A+tauX*X)
        #s  
        s <- sqrt((1+(fm[i]+fc[i])*(phi[i]-1))^2+4*(1-phi[i])*phi[i]*fm[i]*fc[i])
        
        if(phi[i]==1)
        {
          p11 <- fm[i]*fc[i]
          partialphi <- 0
        }
        if(phi[i]!=1)
        {
          p11 <- (1+(fm[i]+fc[i])*(phi[i]-1)-s)/(2*(phi[i]-1))
          term11 <- fm[i]+fc[i]
          term12 <- 0.5*((1+(fm[i]+fc[i])*(phi[i]-1))^2+4*phi[i]*(1-phi[i])*fm[i]*fc[i])^(-0.5)
          term13 <- 2*(1+(fc[i]+fm[i])*(phi[i]-1))*(fm[i]+fc[i])+4*fc[i]*fm[i]-8*fc[i]*fm[i]*phi[i]
          term1 <- (term11-term12*term13)/(2*(phi[i]-1))
          term2 <- p11/(phi[i]-1)
          partialphi <- term1-term2
        }
        partialp11tauX <- partialphi*phi[i]*X[i]
        t <- a[i]/p11-b[i]/(fc[i]-p11)-c[i]/(fm[i]-p11)+d[i]/(1-fc[i]-fm[i]+p11)
        thing[i] <- t * partialp11tauX
      }
      sum(thing)
    }
    
    f <- rep(NA,3)
    
    f[1] <- f1(par)
    f[2] <- f2(par)
    f[3] <- f3(par)
    
    f
  }
  
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
  
  p0 <- c(0.5,-1,2)
  solution <- nleqslv(p0, fn,control = list(trace=T)) 
  
  esttau <- solution$x
  phiest <- exp(esttau[1]+esttau[2]*A+esttau[3]*X)
  p11.est <- rep(NA,n)
  for(i in 1:n)
  {
    p11.est[i] <- p11.fun(phiest[i],fc[i],fm[i])
  }
  print(c(esttau[1]-tau_0,esttau[2]-tau_A,esttau[3]-tau_X))
  return(p11.est)
}

#all correct
#quad(fitm.mar,fitc.mar,1,fity,fita,A,M,C,X,Y)-truth
#wrong EY

#number of replicates
B <- 1000

do.one<-function(trials)
{
  truth <- gen(10000,T)
  dataset <- gen(n,F)
  
  fitm.mar <- glm(M~A+X,data = dataset,family = binomial(link="logit"))
  fitc.mar <- glm(C~A+X,data = dataset,family = binomial(link="logit"))
  fita <- glm(A~X,family = binomial(link = "logit"),data = dataset)
  fity <- lm(Y~A+X+C+M+A*C+A*M,data = dataset)
  p11.est <- estp11(fitm.mar,fitc.mar,dataset$A,dataset$M,dataset$C,dataset$X,dataset$Y,n)
  quad(fitm.mar,fitc.mar,p11.est,fity,fita,dataset$A,dataset$M,dataset$C,dataset$X,dataset$Y,n)-truth
}


library(parallel)
numCores <- detectCores()
trials <- seq(1:B)
result.allcorrect <- mclapply(trials,do.one,mc.cores=numCores)




#result.allcorrect <- pbreplicate(1000,do.one(n))
result2.allcorrect <- matrix(unlist(result.allcorrect),ncol=B)
correct <- apply(result2.allcorrect,1,mean)*100
correct.se <- apply(result2.allcorrect,1,sd)*100
beep()

correct <- matrix(correct,nrow = 1)
colnames(correct) <- c("delta1","delta2","delta3","delta4","delta_quad","triple","tripleaspre")
xtable(correct,digits=1)
correct.se <- matrix(correct.se,nrow = 1)
colnames(correct.se) <- c("delta1","delta2","delta3","delta4","delta_quad","triple","tripleaspre")
xtable(correct.se,digits=1)

write.csv(correct,"comparebinary.csv")
write.csv(correct.se,"comparebinary_se.csv")



