

setwd( R'(C:\Users\James.Thorson\Desktop\Git\autodiff\SDM)' )

# 
library(RTMB)
library(fmesher)
library(Matrix)
library(sf)
library(viridisLite)
library(rnaturalearth)

#
dat = read.csv( "data_real.csv" )
dat_sf = st_as_sf( dat, coords = c("Lon","Lat") )

# Get from shapefile
domain = st_read( file.path(shape_dir,"EBSshelf.shp") )
domain = st_geometry(domain)
domain = st_transform( domain, crs=st_crs("EPSG:4326") )
grid = st_make_grid( domain, cellsize=c(0.25,0.25) )
grid = st_intersection(grid, domain)
grid = st_make_valid(grid)
extrap = st_coordinates(st_centroid(grid))
extrap = cbind( 'Lon'=extrap[,1], 'Lat'=extrap[,2], 'Area_in_survey_km2'=as.numeric(st_area(grid)) / 1e6 )

# Get land map
alaska = ne_states( country = "united states of america" )
alaska = subset( alaska, gn_name == "Alaska" )

#
mesh = fm_mesh_2d( dat[,c('Lon','Lat')], cutoff = 0.5 )
spde = fm_fem( mesh, order=2 )
A_gs = fm_evaluator( mesh, loc=as.matrix(extrap[,c('Lon','Lat')]) )$proj$A
area_g = extrap[,'Area_in_survey_km2']

# Extract
M0 = spde$c0
M1 = spde$g1
M2 = spde$g2

parlist = list(
  mu_c = rep(0, 3),
  omega_sc = array(0, dim=c(mesh$n,3) ),
  ln_kappa = log(1),
  ln_tau = log(1),
  ln_q = log(1),
  ln_phi = log(1),
  invf_p = 0
)

############ START RTMB CODE BLOCK
jnll_spde = function( parlist, what="jnll" ){
  getAll(parlist)
  Q = (exp(4*ln_kappa)*M0 + 2*exp(2*ln_kappa)*M1 + M2) * exp(2*ln_tau)
  omega_ic = A_is %*% omega_sc
  # Likelihood terms
  nll_prior = nll_data = nll_omega = 0
  for( i in seq_along(b_i) ){
    if(Gear[i] == "BT"){
      yhat = exp( ln_q + mu_c[1] + omega_ic[i,1] ) +
             exp( ln_q + mu_c[2] + omega_ic[i,2] )
    }
    if(Gear[i] == "AT2") yhat = exp( mu_c[2] + omega_ic[i,2])
    if(Gear[i] == "AT3") yhat = exp( mu_c[3] + omega_ic[i,3])
    nll_data = nll_data -
      dtweedie( x=b_i[i], mu=yhat, phi=exp(ln_phi),
                p=plogis(invf_p) + 1, log=TRUE)
  }
  for( c_index in 1:3 ){
    nll_omega = nll_omega - dgmrf( omega_sc[,c_index], Q = Q, log=TRUE )
  }
  nll_prior = -1 * dnorm( ln_q, mean=0, sd=0.15, log=TRUE )

  # Report density
  omega_gc = A_gs %*% omega_sc
  D_gc = outer(area_g,rep(1,3)) * exp(outer(rep(1,length(grid)),mu_c) + omega_gc)
  REPORT( D_gc )
  return( nll_data + nll_omega + nll_prior )
}
############ END RTMB CODE BLOCK

#
map = list()
  map$ln_q = factor(NA)
build_obj = function(){
  MakeADFun( 
    func = jnll_spde,
    par = parlist,
    random = "omega_sc",
    silent = TRUE,
    map = map
  )
}
years_to_run = c(2007, 2010)

# Unpack and run
dat_1 = subset( dat, Year == years_to_run[1] )
b_i = dat_1$Catch_KG
Gear = dat_1$Gear
A_is = fm_evaluator( mesh, loc=as.matrix(dat_1[,c('Lon','Lat')]) )$proj$A
obj_1 = build_obj()
opt_1 = nlminb( obj_1$par, obj_1$fn, obj_1$gr,
              control = list(iter.max=1e4, eval.max=1e4, trace=1))
sdrep_1 = sdreport( obj_1 )
D1_gc = obj_1$report()$D_gc

# Unpack and run
dat_2 = subset( dat, Year == years_to_run[2] )
b_i = dat_2$Catch_KG
Gear = dat_2$Gear
A_is = fm_evaluator( mesh, loc=as.matrix(dat_2[,c('Lon','Lat')]) )$proj$A
obj_2 = build_obj()
opt_2 = nlminb( obj_2$par, obj_2$fn, obj_2$gr,
              control = list(iter.max=1e4, eval.max=1e4, trace=1))
sdrep_2 = sdreport( obj_2 )
D2_gc = obj_2$report()$D_gc

#
plotgrid = st_sf( grid,
                  Ptrawl_2007 = rowSums(D1_gc[,1:2]) / rowSums(D1_gc),
                  Paccoustic_2007 = rowSums(D1_gc[,2:3]) / rowSums(D1_gc),
                  Ptrawl_2018 = rowSums(D2_gc[,1:2]) / rowSums(D2_gc),
                  Paccoustic_2018 = rowSums(D2_gc[,2:3]) / rowSums(D2_gc) )

source( "add_legend.R" )
png( file="Proportion_of_total_density.png", width=7.5, height=6, units="in", res=200 )
  par( mar = c(0,0,0,0), mfrow = c(2,2), oma=c(2,2,3,3) )
  for( t in 1:4 ){
    plot( plotgrid[,t], border=NA, breaks = seq(0,1,length=21), pal=viridis, key.pos=NULL, reset=FALSE, main="" )
    plot( st_geometry(alaska), add = TRUE, col = "grey" )
    box()
    if(t==1) add_legend( seq(0,1,length=6), legend_y=c(0.6,1), legend_x=c(0.95,1), col=viridis(10) )
    if(t==1) mtext( side=3, text = "Bottom trawl" )
    if(t==2) mtext( side=3, text = "Accoustics" )
    if(t==2) mtext( side=4, text = years_to_run[1], line=0.5 )
    if(t==4) mtext( side=4, text = years_to_run[2], line=0.5 )
    if(t==1) axis(2)
    if(t==3) axis(2)
    if(t==3) axis(1)
    if(t==4) axis(1)
  }
  mtext( side=c(1,2,3,4), outer=TRUE, text=c("Longitude","Latitude","Proportion available to gear","Year"), line=1, font = 2 )
dev.off()
