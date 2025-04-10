

setwd( R'(C:\Users\James.Thorson\Desktop\Git\autodiff\Utility)')

library(RTMB)

#####################
# Functions
#####################

# Function to calculate fishing rate policy given biomass
policy <-
function( B,
          list_policy,
          list_bio ){

  Brate = exp(list_policy$ln_Brate)
  Bmid = plogis(list_policy$invlogis_Kratio) * list_bio$K
  Fmax = exp(list_policy$ln_Fmax)
  F_out = Fmax * plogis( Brate * (B - Bmid) / list_bio$K )
  return(F_out)
}

# Project dynamics given biology, policy, process and implementation errors
project <-
function( list_bio,
          list_policy,
          seed = NULL ){

  # Store sequence of updates
  getAll( list_bio, list_policy )
  P_t = F_t = C_t = B_t = rep(0, n_years)
  B_t[1] = list_bio$K

  # Simulate process and implementation errors
  set.seed(seed)
  eps = rnorm(n_years, sd = sigmaB)
  delta = rnorm(n_years, sd = sigmaF)

  # Project each year in sequence
  for( yr in seq_len(n_years) ){
    # Production
    P_t[yr] = r * B_t[yr] * (1 - B_t[yr] / K)
    P_t[yr] = P_t[yr] * exp(eps[yr] - sigmaB^2 / 2)
    # Fishing mortality
    F_t[yr] = policy( B_t[yr] + P_t[yr], list_policy, list_bio )
    F_t[yr] = F_t[yr] * exp(delta[yr] - sigmaF^2 / 2)
    # Catch
    C_t[yr] = (B_t[yr] + P_t[yr]) * (1 - exp(-F_t[yr]))
    # Update biomass at end of year
    if(yr < n_years) B_t[yr+1] = B_t[yr] + P_t[yr] - C_t[yr]
  }
  return(mean(C_t))
}

# Sample realized catches
sample_catch <-
function( list_policy ){

  C_z = sapply( X = seq_len(n_eval) + seed,
                FUN = project,
                list_policy = list_policy,
                list_bio = list_bio )
  return( -1 * log(mean(C_z)) )
}

#####################
# Run model
#####################

n_years = 100
n_eval = 100
seed = 101

# Simulation parameters (fixed)
list_bio = list(
  r = 0.2,  # MSY = rK / 4 = 50
  K = 1000
)

# Policy parameters (to estimate)
list_policy = list(
  ln_Fmax = log(0.05),
  invlogis_Kratio = qlogis(0.5),
  ln_Brate = log(5)
)

# Example run
list_bio$sigmaB = 0.1
list_bio$sigmaF = 0.1
obj <- MakeADFun( sample_catch,
                  list_policy )
opt = nlminb( start = obj$par,
              objective = obj$fn,
              gradient = obj$gr )
opt$par # Optimal policy parameters

####################
# Comparison
####################

# Compare sets of implementation and process error variances
sigmaB_z = c( 0.1, 0.1 )
sigmaF_z = c( 0.1, 0.5 )

#
list_bio$sigmaB = sigmaB_z[1]
list_bio$sigmaF = sigmaF_z[1]
obj <- MakeADFun( sample_catch,
                  list_policy )
opt = nlminb( start = obj$par,
              objective = obj$fn,
              gradient = obj$gr )

# SE doesn't mean anything (objective is not likelihood!), but Hessian should be posdef
H = optimHess( opt$par, obj$fn, obj$gr )
eigen(H)$values

#
list_bio$sigmaB = sigmaB_z[2]
list_bio$sigmaF = sigmaF_z[2]
obj_2 <- MakeADFun( sample_catch,
                    list_policy)
opt_2 = nlminb( start = obj_2$par,
              objective = obj_2$fn,
              gradient = obj_2$gr )

# Visualize production and policy functions
png( file = "utility.png", width=5, height=5, res=200, units="in" )
  par( mfrow=c(2,1), mgp=c(2,0.5,0), mar=c(1,1,0,1), oma=c(2,2,1,2) )
  Pmax = 70
  Fmax = 1
  do_plot = function( obj, sigmaB, sigmaF ){
    parhat = obj$env$parList()
    B = seq(1, 1 * list_bio$K, length=100)
    F = sapply( B, parhat, FUN = policy, list_bio = list_bio )
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

#######################
# Additional exploration
#######################

# Visualize expected catch w.r.t. sigmaF ... not very interesting
J_z2 = rep(NA, 8)
pars_z2z = matrix(NA, nrow=length(J_z2), ncol=3)
sigmaF_z2 = exp( seq(log(0.01),log(1), length=length(J_z2)) )
for(z2 in seq_along(J_z2)){
  message( "Running ", z2, " / ", length(J_z2) )
  list_bio$sigmaF = sigmaF_z2[z2]
  list_bio$sigmaB = 0.1
   obj_z2 <- MakeADFun( sample_catch,
                       list_policy,
                       silent = TRUE )
  opt_z2 = nlminb( start = obj_z2$par,
                objective = obj_z2$fn,
                gradient = obj_z2$gr )
  J_z2[z2] = exp( opt_z2$obj )
  pars_z2z[z2,] = opt_z2$par
}

B = seq(1, 1 * list_bio$K, length=100)
F = apply( pars_z2z, MARGIN=1,
          FUN = \(parvec){
            sapply( B,
                    relist(parvec, list_policy),
                    FUN = policy, list_bio = list_bio )
          } )
png( "utility-2.png", width=4, height=4, res=200, units="in" )
  par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  matplot( x = B, y = F, xlab = "Population biomass", ylab = "Fishing rate policy",
           col = viridis(ncol(F)), xaxs="i", yaxs="i",  ylim = c(0,1),
           type = "l", lty = "solid", lwd = 2 )
  legend( "topleft",
          bty = "n",
          legend = formatC(sigmaF_z2,digits=3,format="f"),
          fill = viridis(ncol(F)),
          title = "sigmaF" )
dev.off()
