

# devtools::install_github( "kaskr/RTMB/RTMBode")

run_dir = R'(C:\Users\James.Thorson\Desktop\Work files\Collaborations\2025 -- Generalized graphical additive mixed model\microcosm)'

#
#Date = Sys.Date()
#  date_dir = file.path(root_dir, paste0("LV_",Date) )
#  dir.create( date_dir )

library(RTMB)
library(readxl)
library(RTMBode)
rk4sys = ecostate::rk4sys

dataset = c( "isle_royale", "lynx_hare", "paramesium" )[3]
integration_method = c( "RK4", "RTMBode" )[1]

if( dataset == "paramesium" ){
  #
  full_dat = read_xls( file.path( run_dir, "VeilleuxMS.exc.xls"), sheet = "fig 11a", skip = 4 )
  # https://jmahaffy.sdsu.edu/courses/f09/math636/lectures/lotka/qualde2.html
  full_dat = cbind( X = as.data.frame(full_dat)[,'paramecium...7'] / 100,
                    Y = as.data.frame(full_dat)[,'didinium...8'] / 100 )

}
if( dataset == "lynx_hare" ){
  # https://jmahaffy.sdsu.edu/courses/f09/math636/lectures/lotka/qualde2.html
  hares = c(30, 47.2, 70.2, 77.4, 36.3, 20.6, 18.1, 21.4, 22, 25.4, 27.1, 40.3, 57, 76.6, 52.3, 19.5, 11.2, 7.6, 14.6, 16.2, 24.7)
  lynx = c(4, 6.1, 9.8, 35.2, 59.4, 41.7, 19, 13, 8.3, 9.1, 7.4, 8, 12.3, 19.5, 45.7, 51.1, 29.7, 15.8, 9.7, 10.1, 8.6)
  full_dat = cbind( X = hares / 10, Y = lynx / 10 )
}

#
f = identity
dzdt <-
function( Time,
          State,
          Pars ){

  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  d1_dt = f(Pars$ln_alpha) * State[1] - f(Pars$ln_beta) * State[1] * State[2]
  d2_dt = -1 * f(Pars$ln_gamma) * State[2] + f(Pars$ln_delta) * State[1] * State[2]
  return( c(d1_dt, d2_dt) )
}

get_nll <-
function( parlist ){
  nll_ti = zhat_ti = array(0, dim = dim(dat) )

  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  # Initialize
  zhat_ti[1,] = ( parlist$ln_z0 )

  for( t_index in 2:nrow(dat) ){
    if( integration_method == "RK4" ){
      proj = rk4sys(
            f = dzdt,
            a = 0,
            b = 1,
            n = 10,
            Pars = parlist,
            y0 = parlist$z_ti[t_index-1,] )
      zhat_ti[t_index,] = proj$y[nrow(proj$y),]
    }else{
      proj = ode(
             y = parlist$z_ti[t_index-1,],
             times = c(0,1),
             func = dzdt,
             parms = parlist )
      zhat_ti[t_index,] = proj[nrow(proj),2:3]
    }
  }
  loglik1 = sum(dnorm( dat, mean = parlist$z_ti, sd = (parlist$ln_sigma), log=TRUE ), na.rm=TRUE)
  loglik2 = sum(dnorm( parlist$z_ti, mean = zhat_ti, sd = parlist$ln_tau, log=TRUE ))
  jnll = -1 * (loglik1 + loglik2)
  ADREPORT( zhat_ti )
  return(jnll)
}

#
dat = full_dat
out_of_sample = nrow(dat) - 0:floor(nrow(dat)/4)
dat[ out_of_sample, ] = NA

parlist = list(
  ln_alpha = (1),
  ln_beta = (1),
  ln_gamma = (1),
  ln_delta = (1),
  ln_z0 = (c(0.1,0.1)),
  ln_sigma = (1),
  ln_tau = (1),
  z_ti = ifelse( is.na(dat), rep(1,nrow(dat))%o%colMeans(dat,na.rm=TRUE), dat )
)

get_nll(parlist)
obj = MakeADFun( func = get_nll,
                 parameters = parlist,
                 random = "z_ti" )
obj$gr(obj$par)

opt = nlminb( obj$par, obj$fn, obj$gr, # hessian = obj$he,
              control = list(iter.max = 1e4, eval.max = 1e4),
              lower = 0.01 )
sdrep = sdreport( obj )

zhat_ti = as.list( sdrep, what = "Estimate", report=TRUE )$zhat_ti
zse_ti = as.list( sdrep, what = "Std. Error", report=TRUE )$zhat_ti

#
png( file = file.path(run_dir,paste0(dataset,".png")), width=5, height = 5, res=200, units="in" )
  par( mar = c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02, xaxs = "i", yaxs = "i" )
  x = seq_len(nrow(full_dat))
  matplot( x = x, y = full_dat, type="p", col = c("black","red"), xlab = "time", ylab = "abundance" )
  matplot( x = x, y = zhat_ti, type="l", lty="solid", col = c("black","red"), add=TRUE )
  polygon( x = c(x,rev(x)), y = c(zhat_ti[,1]+zse_ti[,1],rev(zhat_ti[,1]-zse_ti[,1])),
           type="l", lty="solid", col = rgb(0,0,0,0.2), add=TRUE, border=NA )
  polygon( x = c(x,rev(x)), y = c(zhat_ti[,2]+zse_ti[,2],rev(zhat_ti[,2]-zse_ti[,2])),
           type="l", lty="solid", col = rgb(1,0,0,0.2), add=TRUE, border=NA )
  abline( v = min(out_of_sample) - 0.5, lwd=2, lty="dotted" )
  if(dataset=="paramesium") legend( "topleft", bty = "n", legend = c("paramecium","didinium"), fill = c("black","red"))
dev.off()
