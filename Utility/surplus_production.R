

setwd( R'(C:\Users\James.Thorson\Desktop\Git\autodiff\Utility)')

library(RTMB)
set.seed(101)

#####################
# Functions
#####################

policy <-
function( B,
          list_fixed ){

  Brate = exp(list_fixed$ln_Brate)
  Bmid = plogis(list_fixed$invlogis_Kratio) * list_bio$K
  Fmax = exp(list_fixed$ln_Fmax)
  F_out = Fmax * plogis( Brate * (B - Bmid) / list_bio$K )
  return(F_out)
}

#
project <-
function( list_bio,
          list_both ){

  # Storage
  Bstart_t = Bend_t = F_t = C_t = P_t = rep(0, n_years)

  # Given Bstart at beginning of year
  for( t_index in seq_len(n_years) ){
    if(t_index==1){
      Bstart_t[t_index] = list_bio$K
    }else{
      Bstart_t[t_index] = Bend_t[t_index-1]
    }
    # Production
    P_t[t_index] = list_bio$r * Bstart_t[t_index] * (1 - Bstart_t[t_index] / list_bio$K)
    P_t[t_index] = P_t[t_index] * exp(list_bio$sigmaB * list_both$eps[t_index] - list_bio$sigmaB^2 / 2)
    # Fishing mortality
    F_t[t_index] = policy( Bstart_t[t_index] + P_t[t_index], list_both )
    F_t[t_index] = F_t[t_index] * exp(list_bio$sigmaF * list_both$delta[t_index] - list_bio$sigmaF^2 / 2)
    # Catch
    C_t[t_index] = (Bstart_t[t_index] + P_t[t_index]) * (1 - exp(-F_t[t_index]))
    # Bend
    Bend_t[t_index] = Bstart_t[t_index] + P_t[t_index] - C_t[t_index]
  }

  return(C_t)
}

#
sample_catch <-
function( n_eval,
          list_both,
          seed = NULL ){

  rerun = function(list_both, run_num, seed ){
    set.seed(seed + run_num)
    list_both$eps = rnorm(n_years)
    list_both$delta = rnorm(n_years)
    C_t = project( list_bio,
                   list_both )
    mean(C_t)
  }
  C_z = sapply( X = seq_len(n_eval), FUN = rerun, list_both = list_both, seed = seed )
  neg_log_mean_C = -1 * log(mean(C_z))
  return( neg_log_mean_C )
}

# RTMB tape and gradient
negJ_wrt_list_both = function( list_both ){
  sample_catch( n_eval = 100, list_both, seed = 101 )
}

#####################
# Run model
#####################


n_years = 100

# Simulation parameters (fixed)
list_bio = list(
  r = 0.2,       # MSY = rK / 4
  K = 1000,
  sigmaB = 0.1,
  sigmaF = 0.1
)

# Policy parameters (to estimate)
list_fixed = list(
  ln_Fmax = log(0.05),
  invlogis_Kratio = qlogis(0.5),
  ln_Brate = log(5)
)

# Add depending on version
list_both = list_fixed
list_both$eps = rnorm(n_years)     # Process error innovation
list_both$delta = rnorm(n_years)   # Implementation error innovation

# Map off eps & delta (overwritten in sample_catch)
map = list(
  eps = factor(rep(NA,n_years)),
  delta = factor(rep(NA,n_years))
)

#
list_bio$sigmaB = 0.1
list_bio$sigmaF = 0.1
obj <- MakeADFun( negJ_wrt_list_both,
                                     map = map,
                                     list_both )
opt = nlminb( start = obj$par,
              objective = obj$fn,
              gradient = obj$gr,
              control = list(trace=1) )

# SE doesn't mean anything (objective is not likelihood!), but Hessian should be posdef
H = optimHess( opt$par, obj$fn, obj$gr )
eigen(H)$values

#
list_bio$sigmaB = 0.1
list_bio$sigmaF = 0.5
obj_2 <- MakeADFun( negJ_wrt_list_both,
                                     map = map,
                                     list_both)
opt_2 = nlminb( start = obj_2$par,
              objective = obj_2$fn,
              gradient = obj_2$gr,
              control = list(trace=1) )

# Visualize production and policy functions
png( file = "utility.png", width=5, height=5, res=200, units="in" )
  par( mfrow=c(2,1), mgp=c(2,0.5,0), mar=c(1,1,0,1), oma=c(2,2,1,2) )
  Pmax = 70
  Fmax = 1
  do_plot = function( obj, sigmaB ){
    parhat = obj$env$parList()
    B = seq(1, 1 * list_bio$K, length=100)
    F = sapply( B, parhat, FUN = policy )
    P = list_bio$r * B * (1 - B / list_bio$K)
    Plow = P * exp(1.95 * sigmaB)
    Phigh = P * exp(-1.95 * sigmaB)
    plot( x=B, y=P, type = "l", ylim = c(0,Pmax), col = "blue", lwd=2,
          xlab="", ylab = "", xaxt="n", yaxs="i", xaxs="i" )
    polygon( x = c(B,rev(B)), y = c(Plow,rev(Phigh)), col = rgb(0,0,1,0.2), border=NA )
    lines( x=B, y=F / Fmax * Pmax, type = "l", lwd=2, col = "red" )
    axis( 4, labels = pretty(seq(0,Fmax,length=5)),
             at = pretty(seq(0,Fmax,length=5)) / Fmax * Pmax )
  }
  do_plot( obj, sigmaB = 0.1 )
  do_plot( obj_2, sigmaB = 0.1 )
  axis(1)
  mtext( side=c(2,4), outer=TRUE, line=c(1,1), text=c("Production","Fishing mortality") )
dev.off()
