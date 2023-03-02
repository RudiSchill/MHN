

#Convert a state from a bit-vector to a natural number.
State.to.Int <- function(x){
  x <- as.logical(rev(x))
  packBits(rev(c(rep(FALSE, 32 - length(x)%%32), x)), type="integer") + 1
}

#Convert a data matrix, where each row is the bit-vector of a state,
#to a probability distribution as a vector.
Data.to.pD <- function(Data){
  Data <- as.matrix(Data)
  n <- ncol(Data)
  N <- 2^n
  
  Data <- apply(Data, 1, State.to.Int)

  pD <- tabulate(Data, nbins=N)
  pD <- pD/sum(pD)
  
  return(pD)
}

#Simulate an empirical sample from a probability distribution.
Finite.Sample <- function(pTh, k){
  N <- length(pTh)
  tabulate(sample(1:N, k, prob=pTh, replace=T), nbins=N) / k
}

#Kullback-Leibler divergence from model distribution q to true distribution p
KL.Div <- function(p,q){
  as.numeric(p%*%log(p) - p%*%log(q))
}

#Convenient wrapper function to learn an MHN directly from data
MHN <- function(Dat, lambda=NULL, logarithmic=TRUE){
  if(is.null(lambda)) lambda <- 1/nrow(Dat)
  
  pD <- Data.to.pD(Dat)
  Theta <- Learn.MHN(pD, lambda=lambda)
  
  colnames(Theta) <- colnames(Dat)
  rownames(Theta) <- colnames(Dat)
  
  return(Theta)
}

#Generate a single synthetic data point from a given MHN
SingleDatapoint <- function(theta){
  n <- ncol(theta)
  x <- rep(0, n)
  b <- diag(theta)
  diag(theta) <- 0  
  
  observation_time <- rexp(1,1)
  
  t <- 0
  for(i in 1:n){
    rates <- exp(as.vector(theta%*%x) + b)
    rates[which(x==1)] <- 0
    Z <- sum(rates)
    
    t <- t + rexp(1, rate=Z)
    if(t >= observation_time) break;
    
    next_event <- sample(1:n, size=1, prob=rates/Z)
    x[next_event] <- 1
  }
  
  return(x)
}

#Generate a synthetic dataset from an MHN
SampleMHN <- function(theta, size=1000){
  Dat <- matrix(0, nrow=size, ncol=ncol(theta))
  colnames(Dat) <- colnames(theta)
  
  for(i in 1:size){
    Dat[i,] <- SingleDatapoint(theta)
  }
  return(Dat)
}



