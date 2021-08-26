estp11.w <- function(fitm.mar,fitc.mar,M,C,n,dataset.Z)
{
  #psuedo MLE
  
  fc <- expit(predict(fitc.mar,newdata = dataset.Z))
  fm <- expit(predict(fitm.mar,newdata = dataset.Z))
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
      tau0 <- par

      thing <- rep(NA,n)
      
      for(i in 1:n)
      {
        #phi
        phi <- exp(tau0)
        #s  
        s <- sqrt((1+(fm[i]+fc[i])*(phi-1))^2+4*(1-phi)*phi*fm[i]*fc[i])
        
        if(phi==1)
        {
          p11 <- fm[i]*fc[i]
          partialphi <- 0
        }
        if(phi!=1)
        {
          p11 <- (1+(fm[i]+fc[i])*(phi-1)-s)/(2*(phi-1))
          term11 <- fm[i]+fc[i]
          term12 <- 0.5*((1+(fm[i]+fc[i])*(phi-1))^2+4*phi*(1-phi)*fm[i]*fc[i])^(-0.5)
          term13 <- 2*(1+(fc[i]+fm[i])*(phi-1))*(fm[i]+fc[i])+4*fc[i]*fm[i]-8*fc[i]*fm[i]*phi
          term1 <- (term11-term12*term13)/(2*(phi-1))
          term2 <- p11/(phi-1)
          partialphi <- term1-term2
        }
        
        partialp11tau0 <- partialphi*phi
        
        t <- a[i]/p11-b[i]/(fc[i]-p11)-c[i]/(fm[i]-p11)+d[i]/(1-fc[i]-fm[i]+p11)
        thing[i] <- t * partialp11tau0
      }
      
      sum(thing)
    }
    
    f <- f1(par)
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
  
  p0 <- 0.5
  solution <- nleqslv(p0, fn,control = list(trace=T)) 
  
  esttau <- solution$x
  phiest <- exp(esttau)
  p11.est <- rep(NA,n)
  for(i in 1:n)
  {
    p11.est[i] <- p11.fun(phiest,fc[i],fm[i])
  }
 # print(c(esttau[1]-tau_0))
  return(p11.est)
}
estp11 <- function(fitm.mar,fitc.mar,M,C,n)
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
      tau0 <- par
      
      thing <- rep(NA,n)
      
      for(i in 1:n)
      {
        #phi
        phi <- exp(tau0)
        #s  
        s <- sqrt((1+(fm[i]+fc[i])*(phi-1))^2+4*(1-phi)*phi*fm[i]*fc[i])
        
        if(phi==1)
        {
          p11 <- fm[i]*fc[i]
          partialphi <- 0
        }
        if(phi!=1)
        {
          p11 <- (1+(fm[i]+fc[i])*(phi-1)-s)/(2*(phi-1))
          term11 <- fm[i]+fc[i]
          term12 <- 0.5*((1+(fm[i]+fc[i])*(phi-1))^2+4*phi*(1-phi)*fm[i]*fc[i])^(-0.5)
          term13 <- 2*(1+(fc[i]+fm[i])*(phi-1))*(fm[i]+fc[i])+4*fc[i]*fm[i]-8*fc[i]*fm[i]*phi
          term1 <- (term11-term12*term13)/(2*(phi-1))
          term2 <- p11/(phi-1)
          partialphi <- term1-term2
        }
        
        partialp11tau0 <- partialphi*phi
        
        t <- a[i]/p11-b[i]/(fc[i]-p11)-c[i]/(fm[i]-p11)+d[i]/(1-fc[i]-fm[i]+p11)
        thing[i] <- t * partialp11tau0
      }
      
      sum(thing)
    }
    
    f <- f1(par)
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
  
  p0 <- 0.5
  solution <- nleqslv(p0, fn,control = list(trace=T)) 
  
  esttau <- solution$x
  phiest <- exp(esttau)
  p11.est <- rep(NA,n)
  for(i in 1:n)
  {
    p11.est[i] <- p11.fun(phiest,fc[i],fm[i])
  }
  # print(c(esttau[1]-tau_0))
  return(p11.est)
}