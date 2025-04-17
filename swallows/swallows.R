library(RTMB)

data <- readRDS('data.RDS')
# Pick some arbitrary birds to examine survival probabilities.
# Makes sense to pick one from each family from a different year
# and survives more than a few periods
library(dplyr)
set.seed(1021)
meta <- with(data, data.frame(id=1:I, family, last, year))
data$id <- group_by(meta, family) |> filter(last > 5) |> slice_sample(n=1) |>
  group_by(year) |> slice_sample(n=1) |> pull(id)

# initial parameter list
pars <- list(logsigmayearphi=1,
       logsigmaphi=-1,
       logsigmap=0,
       a=rep(1,data$K-1), a1=.5,
       b0=rep(3,4), b1=rep(0,4),
       fameffphi_raw=rep(.1,data$nfam),
       fameffp_raw=rep(.1, data$nfam),
       yeareffphi_raw=rep(.1,4))

# The RTMB model
func <- function(pars){
  getAll(pars,data)
  p <- matrix(1, I, K)
  phi <- matrix(0, I, K-1)
  chi <- matrix(1, I, K+1)
  sigmayearphi <- exp(logsigmayearphi)
  sigmaphi <- exp(logsigmaphi)
  sigmap <- exp(logsigmap)
  
  # priors
  nlp <- -sum(dnorm(b0, 0, 5, TRUE))-
    sum(dnorm(b1, 0, 5, TRUE))-
    sum(dnorm(a, 0, 1.5, TRUE))-
    dnorm(a1, 0, 5, TRUE)-
    dcauchy(sigmaphi, 0, 1, TRUE)-
    dnorm(sigmayearphi, 0, 2, TRUE)-
    dcauchy(sigmap, 0, 1, TRUE)
  # Jacobians and hyper distributions
  nll <- -logsigmaphi-logsigmayearphi-logsigmap-
    sum(dnorm(fameffphi_raw, 0, 1, TRUE))-
    sum(dnorm(fameffp_raw,0, 1, TRUE))-
    sum(dnorm(yeareffphi_raw, 0, 1, TRUE))
  for(i in 1:I){
    phi[i,1:(K-1)] <- plogis(a[1:(K-1)] + a1*carez[i] +
                              sigmayearphi*yeareffphi_raw[year[i]] +
                              sigmaphi*fameffphi_raw[family[i]])
    p[i,2:K] <- plogis(b0[year[i]]+ b1[year[i]]*agec[2:K]+
                         sigmap*fameffp_raw[family[i]])
    for(k in K:2){
      chi[i,k] <- (1-phi[i,k-1]) + phi[i,k-1] * (1-p[i,k])*chi[i,k+1]
    }
    chi[i,1] <- (1-p[i,1]) * chi[i,2]
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
  prep <- p[id,] # report a subset of individuals
  REPORT(prep)
  logit_prep <- log(prep)-log(1-prep)
  ADREPORT(logit_prep)
  REPORT(nld)
  return(nld)
}
func(pars)

obj <- MakeADFun(func, pars, silent=TRUE,
                 random=c("fameffphi_raw", "fameffp_raw", "yeareffphi_raw"))
opt <- with(obj, nlminb(par,fn,gr))
sum(opt$par) #  14.41553
# confidence intervals of survival probailities of selected birds
prep.ml <- summary(sdreport(obj), 'report') |> as.data.frame() |>
  setNames(c('Estimate', 'SE')) |> 
  # CI in logit space
  mutate(lwr=Estimate-1.96*SE, upr=Estimate+1.96*SE) |>
  # CI in natural space
  mutate(Estimate=plogis(Estimate), lwr=plogis(lwr), upr=plogis(upr)) |>
  # add ID data for plotting
  mutate(id=rep(data$id, times=data$K), k=rep(1:data$K, each=length(data$id)))


# Run MCMC on the model
library(adnuts)
mcmc <- sample_sparse_tmb(obj, iter=1500, warmup=500,  seed=1231,
                          control=list(adapt_delta=.99),
                          globals=list(data=data, pars=pars))

# loop through each posterior sample and call report to get
# posterior distributions
post <- extract_samples(mcmc)
prep.all <- list()
for(i in 1:nrow(post)){
  tmp <- obj$report(as.numeric(post[i,]))
  prep.all[[i]] <- reshape2::melt(tmp$prep) |> 
    setNames(c('id', 'k', 'surv_prob')) |>
    cbind(iter=i)
}
prep <- do.call(rbind, prep.all) 
prep$id <- data$id[prep$id]

prep.bayes <- group_by(prep, id, k) |>
  summarize(Estimate=median(surv_prob),
         lwr=quantile(surv_prob, .025),
         upr=quantile(surv_prob, 0.975), .groups='drop')

library(ggplot2)
theme_set(theme_bw())
# Bayesian credible intervals vs asymptotic confidence intervals
g <- ggplot(prep.bayes, aes(k, y=Estimate, ymin=lwr, ymax=upr)) + 
  geom_line(color=gray(.5)) +
  geom_ribbon(alpha=.5) + facet_wrap('id') +
  geom_pointrange(data=prep.ml, fatten = 1, color=1) + 
  labs(x='Capture period', y='Survival probability')
ggsave('swallows.png', g, width=6, height=4, units='in')

