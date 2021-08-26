
data_genZ <- function(n)
{
  U1 <- U2 <- U3 <- U4 <- rep(NA,n)
  for(i in 1:n)
  {
    U.temp <- mvrnorm(1,mu=cbind(0,0,0,0), Sigma=diag(c(1,1,1,1)))
    U1[i] <- U.temp[1]
    U2[i] <- U.temp[2]
    U3[i] <- U.temp[3]
    U4[i] <- U.temp[4]
  }
  U <- cbind(U1,U2,U3,U4)
  pa <- expit(U %*% aU)
  A <- rbinom(n,1,prob = pa)
  
  pc <- expit(c0 + cA*A + U %*% cU)
  pm <- expit(m0 + mA*A + U %*% mU)
  phi <- exp(tau_0+tau_A*A+U %*% tau_X)
  
  
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
  
  EY <- y0 + yA*A + yM*M + yC*C + U %*% yU 
  
  Y <- rnorm(n,EY,sdy)
  
  Z1 <- exp(U1/2)
  Z2 <- U2/(1+exp(U1))+10
  Z3 <- (U1*U3/25+0.6)^3
  Z4 <- (U2+U4+20)^2
  
  
  dataset <- data.frame(Z1,Z2,Z3,Z4,A,C,M,Y,U1,U2,U3,U4)
  
  dataset
}
