

setwd( R'(C:\Users\James.Thorson\Desktop\Git\autodiff\Utility)')

library(RTMB)
set.seed(101)

#####################
# Functions
#####################

# Function to calculate fishing rate policy given biomass
policy <-
function( B,
          list_policy ){

  Brate = exp(list_policy$ln_Brate)
  Bmid = plogis(list_policy$invlogis_Kratio) * list_bio$K
  Fmax = exp(list_policy$ln_Fmax)
  F_out = Fmax * plogis( Brate * (B - Bmid) / list_bio$K )
  return(F_out)
}

# Project dynamics given biology, policy, process and implementation errors
project <-
function( list_bio,
          list_both ){

  # Store sequence of updates
  P_t = F_t = C_t = B_t = rep(0, n_years)
  Bstart = list_bio$K

  # Given Bstart at beginning of year
  for( t_index in seq_len(n_years) ){
    # Production
    P_t[t_index] = list_bio$r * Bstart * (1 - Bstart / list_bio$K)
    P_t[t_index] = P_t[t_index] * exp(list_bio$sigmaB * list_both$eps[t_index] - list_bio$sigmaB^2 / 2)
    # Fishing mortality
    F_t[t_index] = policy( Bstart + P_t[t_index], list_both )
    F_t[t_index] = F_t[t_index] * exp(list_bio$sigmaF * list_both$delta[t_index] - list_bio$sigmaF^2 / 2)
    # Catch
    C_t[t_index] = (Bstart + P_t[t_index]) * (1 - exp(-F_t[t_index]))
    # Bend
    B_t[t_index] = Bstart + P_t[t_index] - C_t[t_index]
    Bstart = B_t[t_index]
  }

  return(C_t)
}

# Sample realized catches
sample_catch <-
function( n_eval,
          list_policy,
          seed = NULL ){

  rerun = function(list_policy, run_num, seed ){
    set.seed(seed + run_num)
    list_policy$eps = rnorm(n_years)
    list_policy$delta = rnorm(n_years)
    C_t = project( list_bio,
                   list_policy )
    mean(C_t)
  }
  C_z = sapply( X = seq_len(n_eval), FUN = rerun, list_policy = list_policy, seed = seed )
  neg_log_mean_C = -1 * log(mean(C_z))
  return( neg_log_mean_C )
}

# Calculate neg-log-utility given policy parameters
negJ_wrt_list_both = function( list_policy ){
  sample_catch( n_eval = 100, list_policy, seed = 101 )
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
list_policy = list(
  ln_Fmax = log(0.05),
  invlogis_Kratio = qlogis(0.5),
  ln_Brate = log(5)
)

####################
# Comparison
####################

sigmaB_z = c( 0.1, 0.1 )
sigmaF_z = c( 0.1, 0.5 )

#
list_bio$sigmaB = sigmaB_z[1]
list_bio$sigmaF = sigmaF_z[1]
obj <- MakeADFun( negJ_wrt_list_both,
                  list_policy )
opt = nlminb( start = obj$par,
              objective = obj$fn,
              gradient = obj$gr,
              control = list(trace=1) )

# SE doesn't mean anything (objective is not likelihood!), but Hessian should be posdef
H = optimHess( opt$par, obj$fn, obj$gr )
eigen(H)$values

#
list_bio$sigmaB = sigmaB_z[2]
list_bio$sigmaF = sigmaF_z[2]
obj_2 <- MakeADFun( negJ_wrt_list_both,
                    list_policy)
opt_2 = nlminb( start = obj_2$par,
              objective = obj_2$fn,
              gradient = obj_2$gr,
              control = list(trace=1) )

## Visualize curve ... not very interesting
#J_z2 = rep(NA, 10)
#sigmaB_z2 = exp( seq(log(0.01),log(1), length=length(J_z2)) )
#for(z2 in seq_along(J_z2)){
#  list_bio$sigmaF = 0.1
#  list_bio$sigmaB = sigmaB_z2[z2]
#  obj_z2 <- MakeADFun( negJ_wrt_list_both,
#                       map = map,
#                       list_both,
#                       silent = TRUE )
#  opt_z2 = nlminb( start = obj_z2$par,
#                objective = obj_z2$fn,
#                gradient = obj_z2$gr )
#  J_z2[z2] = exp( opt_z2$obj )
#}

# Visualize production and policy functions
png( file = "utility.png", width=5, height=5, res=200, units="in" )
  par( mfrow=c(2,1), mgp=c(2,0.5,0), mar=c(1,1,0,1), oma=c(2,2,1,2) )
  Pmax = 70
  Fmax = 1
  do_plot = function( obj, sigmaB, sigmaF ){
    parhat = obj$env$parList()
    B = seq(1, 1 * list_bio$K, length=100)
    F = sapply( B, parhat, FUN = policy )
    Flow = F * exp(-1 * sigmaF - sigmaF^2 / 2)
    Fhigh = F * exp(1 * sigmaF - sigmaF^2 / 2)
    P = list_bio$r * B * (1 - B / list_bio$K)
    Plow = P * exp(-1 * sigmaB - sigmaB^2 / 2)
    Phigh = P * exp(1 * sigmaB - sigmaB^2 / 2)
    plot( x=B, y=P, type = "l", ylim = c(0,Pmax), col = "blue", lwd=2,
          xlab="", ylab = "", xaxt="n", yaxs="i", xaxs="i" )
    polygon( x = c(B,rev(B)), y = c(Plow,rev(Phigh)), col = rgb(0,0,1,0.2), border=NA )
    lines( x=B, y=F / Fmax * Pmax, type = "l", lwd=2, col = "red" )
    polygon( x = c(B,rev(B)), y = c(Flow,rev(Fhigh)) / Fmax * Pmax, col = rgb(1,0,0,0.2), border=NA )
    axis( 4, labels = pretty(seq(0,Fmax,length=5)),
             at = pretty(seq(0,Fmax,length=5)) / Fmax * Pmax )
  }
  do_plot( obj, sigmaB = sigmaB_z[1], sigmaF = sigmaF_z[1] )
  do_plot( obj_2, sigmaB = sigmaB_z[2], sigmaF = sigmaF_z[2] )
  axis(1)
  margin = c("Population biomass", "Production (+/- process error)", "Fishing rate policy (+/- implementation error)")
  mtext( side=c(1,2,4), outer=TRUE, line=c(1,1,1), text=margin )
  legend( "top", bty="n", fill = c("blue","red"), legend = c("Production", "Fishing rate"), ncol = 2 )
dev.off()
