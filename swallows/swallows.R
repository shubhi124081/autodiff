library(RTMB)

data <- readRDS('data.RDS')
inits <- function()
  list(sigmayearphi=runif(1, 0, 2),
       sigmaphi=runif(1,0,2),
       sigmap=runif(1,0,2),
       a=rnorm(data$K-1, 3, 1), a1=rnorm(1, .5, 1),
       b0=rnorm(4, 3, sd=1), b1=rnorm(4, 0, .15),
       fameffphi_raw=rnorm(data$nfam,0,1),
       fameffp_raw=rnorm(data$nfam,0,1),
       yeareffphi_raw=rnorm(4, 0,1))
set.seed(23254)
inits <- inits()
random <- c("fameffphi_raw", "fameffp_raw", "yeareffphi_raw")


func <- function(pars){
  getAll(pars,data)
  nlp <- 0
  p <- matrix(1, I, K)
  phi <- matrix(0, I, K-1)
  chi <- matrix(1, I, K+1)
  sigmayearphi2 <- exp(sigmayearphi)
  sigmaphi2 <- exp(sigmaphi)
  sigmap2 <- exp(sigmap)
  
  # priors
  nlp <- -sum(dnorm(b0, 0, 5, TRUE))-
    sum(dnorm(b1, 0, 5, TRUE))-
    sum(dnorm(a, 0, 1.5, TRUE))-
    dnorm(a1, 0, 5, TRUE)-
    dcauchy(sigmaphi2, 0, 1, TRUE)-
    dnorm(sigmayearphi2, 0, 2, TRUE)-
    dcauchy(sigmap2, 0, 1, TRUE)
  # Jacobians and hyper distributions
  nll <- -sigmaphi-sigmayearphi-sigmap-
    sum(dnorm(fameffphi_raw, 0, 1, TRUE))-
    sum(dnorm(fameffp_raw,0, 1, TRUE))-
    sum(dnorm(yeareffphi_raw, 0, 1, TRUE))
  for(i in 1:I){
    phi[i,1:(K-1)] <- plogis(a[1:(K-1)] + a1*carez[i] +
                              sigmayearphi2*yeareffphi_raw[year[i]] +
                              sigmaphi2*fameffphi_raw[family[i]])
    p[i,2:K] <- plogis(b0[year[i]]+ b1[year[i]]*agec[2:K]+
                         sigmap2*fameffp_raw[family[i]])
    k <- K
    while (k > 1) {
      chi[i,k] <- (1 - phi[i,k-1]) + phi[i,k-1] * (1 - p[i,k])*chi[i,k+1]
      k <- k - 1
    }
    chi[i,1] <- (1 - p[i,1]) * chi[i,2]
    # likelihood of the data  
    nll <- nll -
      # probability of survival, known alive since k<last
      ifelse(last[i]>1, sum(dbinom(1,1,phi[i,2:last[i]-1],TRUE)),0)-
      # probability of observation given known alive
      sum(dbinom(CH[i,1:last[i]], 1, prob=p[i,1:last[i]], TRUE))-
      #  probability of no observations after time period last
      log(chi[i,last[i]+1])
  }
  nld <- nll+nlp # negative log density
  REPORT(nld)
  REPORT(chi)
  REPORT(p)
  REPORT(phi)
  return(nld)
}

obj <- MakeADFun(func, inits, random=random, silent=TRUE)
obj$report()$nld - func(inits)
opt <- with(obj, nlminb(par,fn,gr))
sdreport(obj)
