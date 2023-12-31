library(INLA); library(raster); library(maptools)
library(gtools); library(sp); library(spdep); library(rgdal)
library(fields); library(mvtnorm); library(gtools); library(geoR)
library(actuar);library(viridisLite);require(grid);require(gridExtra)
require(lattice);require(tidyverse);require(MASS);library(tmap)
library(tmaptools);library(rgdal);library(sf);library(raster)
#install.packages("cartography")
library(cartography) # mapping dedicated package
#install.packages("OpenStreetMap")
#library(OpenStreetMap)

#This script if for 9 areas

path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP000008_UNFPA_PNG/Working/Chris/reviews/sim_study/option2"
out_path <- paste0(path, "/outputs")
file_path <- paste0(path, "/file")


set.seed(500)
model_metrics <- function(obs, pred, upper, lower)
{
  residual = pred - obs
  INACCURACY = mean(abs(residual), na.rm=T)#MAE
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])
  
  output <- list(MAE  = INACCURACY ,
                 RMSE = RMSE,
                 BIAS = abs(BIAS),
                 corr = corr)
  return(output)
}
#model_metrics(obs, pred, upper, lower)


##---Extract province random effects from the best model
prov_ranef <- function(dat, mod)
{
  pro <- rep(1, nrow(dat))
  pr <- mod$summary.random$prov2_ID$mean
  
  uniq <- unique(dat$prov2_ID)
  uniq[1]
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:length(uniq))
    {
      if(dat$prov2_ID[i]==uniq[j]) pro[i] = pr[j]
    }
  }
  pro
}
#ddat$prov_ranef <- prov_ranef(ddat, mod0)

#------------------------------------------------------------------------
####----POSTERIOR SIMULATIONS-----------------------------------
simPop <- function(model, dat, Aprediction, run)
{
  fixedeff  <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  inla.seed = as.integer(runif(1)*.Machine$integer.max)
  #inla.seed =  1657687559
  set.seed(inla.seed)
  print(inla.seed)
  m1.samp <- inla.posterior.sample(run, model, seed = inla.seed ,selection=list(x2=1, x3=1, x4=1,
                                                                                x5=1, x6=1),num.threads="1:1")
  
  sfield_nodes_mean <- model$summary.random$spatial.field['mean']
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  for(i in 1:run)
  {
    #fixedeff[,i] <- model$summary.fixed['Intercept', 'mean'] +
    fixedeff[,i] <-  
      #m1.samp[[i]]$latent[1,] +
      model$summary.fixed['Intercept', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x2'] +
      m1.samp[[i]]$latent[2,] * dat[,'x3'] +
      m1.samp[[i]]$latent[3,] * dat[,'x4'] +
      m1.samp[[i]]$latent[4,] * dat[,'x5'] +
      m1.samp[[i]]$latent[5,] * dat[,'x6'] + 
      dat$prov_ranef +
      model$summary.random$IDsp$mean +
      
      field_mean[,1]
    
    pop_hat[,i]<- exp(fixedeff[,i])
  }
  
  #mean_pop_hat1 <- dat$pop_hat1 #
  mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) #
  median_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.5), na.rm=T) #
  lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) #
  uncert_pop_hat <- (upper_pop_hat - lower_pop_hat)/mean_pop_hat#
  
  
  dat$mean_pop_hat <- mean_pop_hat
  dat$median_pop_hat <- median_pop_hat
  dat$lower_pop_hat <- lower_pop_hat
  dat$upper_pop_hat <- upper_pop_hat
  dat$uncert_pop_hat <- uncert_pop_hat
  dat$sd_pop_hat <- sd_pop_hat
  
  
  
  output <- list(pop_hat = pop_hat,
                 est_data = dat)
  
}
#rm(dat.sim, sim.pops_nb, sim.pops, llg.est, prov.est2)
#run=2000
#run=100
#system.time(str(sim.pop <- simPop(mod0,ddat, A, run)))




#---------------------------------------------------------
####---View a random subset of the Posterior samples
#---------------------------------------------------------

##-----Traceplots 
#dev.off()
traceplot <- function(dat1, dat2, cex.main, lwd, col, file_path)
{
  samp <- sample(nrow(dat1), 6)
  
  jpeg(file_path)
  par(mfrow=c(3,2), mar=c(5,5,2,1))
  for(j in samp)
  {
    plot.ts(dat2[j,], ylab = "pop_hat", col="steel blue",
            cex.main=cex.main,  main=paste0("Pop_hat samples for CU ", dat1$IDsp[j], sep=""))
    abline(h=mean(dat2[j,]), lwd=lwd, col=col)
  }
  dev.off()
}
#trplot <- traceplot(ddat, sim.pop$pop_hat,
#                    cex.main=2, lwd=2, col=2,
#                    file_path = paste0(out_path,"/best_fit_model_traceplots.jpeg"))


##---Histograms
#for(j in samp)
#{
#  hist(sim.pop$pop_hat[j,],xlab="pop_hat", prob=T, ylim=c(0,0.01),
#       cex.axis=1.5, main=paste0("Pop_hat histogram for CU ", 
#                                ddat$IDsp[j], sep=""))
# abline(v=mean(sim.pop$pop_hat[j,]), lwd=2, col="red")
#}


#####----Join the posterior draws to the dataframe------------------------------------------------------


####-----------ADMIN TOTALS AND UNCERTAINTIES----------------------------------------



#----NATIONAL------------------------------
nat_total <- function(dat, thinn)
{
  p_hat <- dat[,thinn]
  tots <- apply(p_hat,2, sum, na.rm=T) #Col sums
  
  tot_sd  <- sd(tots, na.rm=T)
  
  tot_mean  <- mean(tots, na.rm=T)
  
  tot_lower <- quantile(tots, probs=c(0.025))
  tot_median <- quantile(tots, probs=c(0.5))
  tot_upper <- quantile(tots, probs=c(0.975))
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, 
                                                       lower=tot_lower, median=tot_median, upper=tot_upper))))
}

#sum(ddat$pop_hat1 <- fit00)
#dim(ddat)
#names(ddat)
#data.all <- ddat[,c("prov2_ID" ,"pop_hat1")]
#data.all1 <- data.frame(cbind(data.all, sim.pop$pop_hat))#---Contains all the simulated posterior matrix
#dim(data.all1);names(data.all1)
#thin=5 #----Thinning choses every 5th sample for posterior inference 
#thinn <- seq((6+0.2*run),run, by=thin) #--There is a burnin period of the first 20% of the total sample
#(national <- nat_total(data.all1, thinn))

#write.csv(national, file=paste0(out_path, "/best_fit_model_national_estimates.csv"))


##---CU-level 
cu_est <- function(dat, dat2, thinn)
{
  dat$mean  <- apply(dat2[,thinn], 1, mean, na.rm=T)
  dat$lower <- apply(dat2[,thinn], 1, quantile, probs=c(0.025), na.rm=T)
  dat$median <- apply(dat2[,thinn], 1, quantile, probs=c(0.5), na.rm=T)
  dat$upper <- apply(dat2[,thinn], 1, quantile, probs=c(0.975), na.rm=T)
  dat$uncertainty <- (dat$upper-dat$lower)/dat$mean
  return(dat)
}
#cu.est <- cu_est(sim.pop$est_data, data.all1,thinn)
#head(cu.est)
#write.csv(cu.est, file=paste0(out_path, "/best_fit_model_CU_estimates.csv"))


####---- fit metrics
#mod_metrics <- model_metrics(cu.est$pop, cu.est$mean, 
#                             cu.est$upper, cu.est$lower)
#write.csv(mod_metrics, file=paste0(out_path, "/best_fit_model_fit_metrics.csv"))



post_maps <- function(cu.est, shp)
{
  names(cu.est)
  shp.cu <- shp
  shp.cu$obs <- cu.est$popm
  shp.cu$dens_hat <- cu.est$mean_dens_hat
  shp.cu$mean <- cu.est$mean
  shp.cu$lower <- cu.est$lower
  shp.cu$upper <- cu.est$upper
  shp.cu$uncertainty <- cu.est$uncertainty
  
  crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
  shpcu_UTM = spTransform(shp.cu, crs.UTM)
  
  #plot(shpcu_UTM)
  #spplot(shpcu_UTM, "mean")
  
  
  
  tmap_options(check.and.fix = TRUE)
  tmap_mode("plot")
  
  #--Observed counts
  obs <-  tm_shape(shpcu_UTM)+ 
    tm_polygons(col="obs", title="Observed Count",
                style="quantile", palette="viridis")+
    tm_layout(legend.outside = F, legend.text.size=1.5, legend.title.size=2)+
    tm_compass(position = c("right", "bottom"), text.size=1.5)+
    tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 200, 400, 600))+
    tm_layout(main.title = "(a)") 
  
  ##--predicted counts 
  pred <-  tm_shape(shpcu_UTM)+ 
    tm_polygons(col="mean", title="Predicted Count",
                style="quantile", palette="viridis")+
    tm_layout(legend.outside = F, legend.text.size=1.5, legend.title.size=2)+
    tm_compass(position = c("right", "bottom"), text.size=1.5)+
    tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 200, 400, 600))+
    tm_layout(main.title = "(b)") 
  
  
  ##--uncertainties
  uncert <-  tm_shape(shpcu_UTM)+ 
    tm_polygons(col="uncertainty", title="Uncertainty",
                style="quantile", palette="viridis")+
    tm_layout(legend.outside = F, legend.text.size=1.5, legend.title.size=2)+
    tm_compass(position = c("right", "bottom"), text.size=1.5)+
    tm_scale_bar(position = c("left", "bottom"), text.size=1, size=2, breaks=c(0, 200, 400, 600))+
    tm_layout(main.title = "(c)") 
  
  
  t_map1 <-tmap_arrange(obs, pred, uncert, nrow=2)
  return(t_map1)
}
#t_map1 <- post_maps(cu.est, shp)
#tmap_save(t_map1, paste0(out_path,'/best_fit_cu_maps.tiff'), dpi = 300, units="in")


#########---------------------------------------------------------------
###-----Simulate the coordinates based on the lon and lat of PNG
shp <- readOGR(dsn = paste0(path), layer = "PNG_CU_32100_B") #combined shapefile 
plot(shp)

##---Extract the coordinates of the centroids 
lon <- coordinates(shp)[,1]
lat <- coordinates(shp)[,2]
coords <- cbind(lon, lat)
plot(coords)


###
dim(dat.all <- as.data.frame(coords))
head(dat.all)

###----Construct the Mesh, SPDE object and projection matrix
dim(coords)
bnd <- inla.nonconvex.hull(as.matrix(coords),-0.03, -0.05, resolution=c(100,100))
mesh <- inla.mesh.2d(boundary = bnd, max.edge=c(0.6,4), cutoff=0.4)

##
par(mfrow=c(1,1))
plot(mesh)
points(coords, col="red")
mesh$n

####----------------------------------------------

##---spde pARAMETERS
r0 <- 0.3
nu <- 1
sigma0 <- 1
kappa0 <- sqrt(8*nu)/r0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)

spde <- inla.spde2.matern(mesh, 
                          B.tau = matrix(c(log(tau0), -1, +1), nrow=1, ncol=3),
                          B.kappa = matrix(c(log(kappa0), 0, -1), nrow=1, ncol=3),
                          theta.prior.mean = c(0,0),
                          theta.prior.prec = c(0.1, 0.1))


#
Q <- inla.spde2.precision(spde=spde, theta=c(0,0))


#---Simulate the GMRF
sam <- as.vector(inla.qsample(
  n = 1, Q = Q, seed=100))
#length(sam)


###---Build projector matrix A
A <- inla.spde.make.A(mesh=mesh, loc=coords);dim(A)

##----specify the observation indices for estimation 
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)

##---Spatial Random effects
S.pred <- as.vector(A %*% sam)
#hist(S.pred)




#----aDD Province and sETTLEMENT BTYPES UNIFORMLY AT RANDOM 
n.province <- 24
n.set_typ <- 3
n.sample <- nrow(dat.all)

##---Equally likely provinces and settlement types
table(prov1 <- sample(1:n.province,n.sample, prob=rep(1/n.province,n.province), rep=T))
table(set_typ1 <- sample(1:n.set_typ,n.sample, prob=rep(1/n.set_typ,n.set_typ), rep=T))


##---Not Equally likely provinces and settlement types
prop.prov <- c(0.037, 0.031, 0.053, 0.031, 0.049, 0.105, 0.032, 0.018, 0.068,0.036,
               0.055, 0.008, 0.041, 0.088, 0.021, 0.024, 0.014, 0.031, 0.028,0.100, 
               0.023, 0.037, 0.028, 0.042)

prop.set_typ <- c(0.084, 0.8, 0.116)
table(prov2 <- sample(1:n.province,n.sample, prob=prop.prov, rep=T))
table(set_typ2 <- sample(1:n.set_typ,n.sample, prob=prop.set_typ, rep=T))



###---Add to dataset
table(dat.all$clusterID <- as.factor(1:nrow(coord)))
table(dat.all$prov1_ID <- as.factor(prov1))
table(dat.all$prov2_ID <- as.factor(prov2))
table(dat.all$set_typ1_ID <- as.factor(set_typ1))
table(dat.all$set_typ2_ID <- as.factor(set_typ2))

##----------------------------------
###-----Simulate covariates ----------------
nn <- nrow(dat.all)
covs = cbind(1,runif(nn),abs(rnorm(nn)),rpois(nn,2),abs(rnorm(nn)), runif(nn))#--GEOSPATIAL COVARIATES
dim(zpred <- covs)
dim(dat.all)

##--Add to dataset 
dim(ddat <- cbind(dat.all, covs))
ddat <- data.frame(ddat)
dim(ddat)
head(ddat)


names(ddat)[7:12] <- paste0("x", 1:6)
#head(ddat)
#names(ddat)

###############
#Parameters 

##-Nugget effect/iid term
##-Population Count 
#mean_pop <- 417
#var_pop  <- 411^2


##----Nugget effect
sigma_e <- 0.15
eps <- rnorm(nrow(ddat), 0, sigma_e) 


###---Simulate Population count
betaP <- c(5.65, 0.01, 0.12, 0.02, 0.003, 0.012) #--betas - fixed effects

pop <- lambdaP <- numeric(nn)
for (i in 1:nn)
{
  lambdaP[i] <- exp(zpred[i,1]*betaP[1] + zpred[i,2]*betaP[2] + 
                      zpred[i,3]*betaP[3] + zpred[i,4]*betaP[4] + 
                      zpred[i,5]*betaP[5] + zpred[i,6]*betaP[6] + S.pred[i]  + eps[i])
  pop[i] <- rpois(1, lambdaP[i])
}
pop
hist(pop); hist(log(pop))
mean(pop); var(pop)

###--------Add to dataset
ddat$pop <- pop


###---Fit the full model and make predictions

#--settlement - type and province nesting
Zsp <- as(model.matrix( ~ 0 + prov2_ID:set_typ2_ID, data = ddat), "Matrix") 
ddat$IDsp <- 1:nrow(ddat)
ddat$set_prov <- as.factor(apply(Zsp, 1, function(x){names(x)[x == 1]}))#--nesting

#---------------------------------------------------------
####---GLM-based covariates selection--------
#-------------------------------------------------------------------------
library(car) ##--For calculating variance inflation factor (vif)
library(dplyr)
library(tidyverse)

#######-------For Density================================================================
covs_pois<- ddat[,c("ldens", "x2", "x3", "x4", "x5", "x6")]
names(covs_density)
fpois <- glm(formula = ldens ~ x2 + x3 + x4 + x5 + x6, data = covs_pois, family = poisson)
summary(pois)

step_pois <- stepAIC(fpois, scale = 0,
                     direction = c("both"),
                     trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                     k = 2)
step_pois

#--VIF
vif_pois = vif(fpois)
vif_pois #---all vif less than 2, multicollinearity not detected
######---------------------------------------------------------------------------------------- 



####----------------PROCEED WITH Bayesian Hierarchical Models in INLA---------------------------

#######---------------------------------------------------------------------



coverp <- coverb <- c(1, 0.8, 0.6, 0.4, 0.2)# - 1  means not missing
for(i in 1:length(coverp))
{
  print(i)
  result_path <- paste0(out_path,"outputs_for_", coverp[i]*100,"%","_pop_count/")
  if (file.exists(result_path)){
    setwd(file.path(result_path))
  } else {
    dir.create(file.path(result_path))
    setwd(file.path(result_path))
  }
  
  
  pm <- coverp[i]
  set.seed(700)
  if(pm == 1){
    ddat$popm <- ddat$pop
  }
  if(pm < 1)
  {
    ind.mpop <- sample(nrow(ddat), pm*nrow(ddat))#pm*100% missing data
    ddat$popm <- ddat$pop
    ddat$popm[-ind.mpop] = NA#--set all missing pop to NA to be predicted 
  }
  
  print(pm*100)
  png("histogram.png")
  hist(log(ddat$popm), xlab = "Log(population count)", main= "")
  dev.off()
  
  covars <- ddat[,c("x1","x2","x3","x4","x5", "x6","set_prov", 
                    "prov2_ID","set_typ2_ID", "IDsp")]; dim(covars)
  stk_est<- inla.stack(data=list(y=ddat$popm), #the response
                       
                       A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                       
                       effects=list(c(list(Intercept=1), #the Intercept
                                      iset),  #the spatial index
                                    #the covariates
                                    list(covars)),
                       #this is a quick name so you can call upon easily
                       tag='est')
  
  
  form0 <- y ~ -1 + Intercept +  x2 + x3 + x4 + x5 + x6  + f(prov2_ID, model='iid') + 
    f(spatial.field, model=spde) + f(IDsp, model='iid') 
  
  
  ###--Poisson
  mod0 <-inla(form0, #the formula
              data=inla.stack.data(stk_est,spde=spde),  #the data stack
              family= 'poisson',   #which family the data comes from
              control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
              control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
              verbose = FALSE)
  
  ind0 <-inla.stack.index(stk_est, "est")$data #--estimation indices 
  fit0<- exp(mod0$summary.linear.predictor[ind0,"mean"]) #--extract the backtransformed pop_hat
  fit0U <- exp(mod0$summary.linear.predictor[ind0,"0.975quant"]) #--extract the backtransformed pop_hat upper
  fit0L <- exp(mod0$summary.linear.predictor[ind0,"0.025quant"]) #--extract the backtransformed pop_hat lower
  
  fit00 <- fit0
  (POP0 <- cbind(ddat$pop, ddat$popm,fit00))
  apply(POP0, 2, sum, na.rm=T)
  plot(ddat$pop, fit0, col=c("red", "blue"))
  abline(a=0, b=1)
  cor(ddat$pop, fit0)
  
  ##---------------------------------------------------------------------------------------  
  #--Run posterior Simulation-------------------------------------------------------------------
  #--extract posterior random effects
  ddat$prov_ranef <- prov_ranef(ddat, mod0)
  
  #--run posterior simulation
  run=2000
  system.time(str(sim.pop0 <- simPop(mod0,ddat, A, run)))
  ###---------------------------
  trplot <- traceplot(ddat, sim.pop0$pop_hat,
                      cex.main=2, lwd=2, col=2,
                      file_path = "traceplot_0.png")
  ##---------------------------
  
  
  sum(ddat$pop_hat0 <- fit00)
  dim(ddat)
  names(ddat)
  data.all <- ddat[,c("prov2_ID" ,"pop_hat0")]
  data.all0 <- data.frame(cbind(data.all, sim.pop0$pop_hat))#---Contains all the simulated posterior matrix
  dim(data.all0);names(data.all0)
  
  ##--national estimates
  thin=5 #----Thinning choses every 5th sample for posterior inference 
  thinn <- seq((6+0.2*run),run, by=thin) #--There is a burnin period of the first 20% of the total sample
  (national0 <- nat_total(data.all0, thinn))
  write.csv(national0, file= "national_estimates_pois.csv")
  
  ##---cu-level
  cu.est0 <- cu_est(sim.pop0$est_data, data.all0,thinn)
  head(cu.est0)
  write.csv(cu.est0, file="CU_estimates_pois.csv")
  
  ##-Metrics
  mod_metrics0 <- model_metrics(cu.est0$pop, cu.est0$mean, 
                                cu.est0$upper, cu.est0$lower)
  write.csv(mod_metrics0, file= "fit_metrics_pois.csv")
  
  ###---maps
  #t_map0 <- post_maps(cu.est1, shp)
 # tmap_save(t_map0, paste0(out_path,'/posterior_maps_pois.tiff'), dpi = 300, units="in")  
  
  
  
  #####-----------------------------------------
  ###---refit density model densi with predicted building------
  ###------------------------------------------------
  #--Negative binomial
  
  form1 <- y ~ -1 + Intercept +  x2 + x3 + x4 + x5 + x6  + f(prov2_ID, model='iid') + 
    f(spatial.field, model=spde) + f(IDsp, model='iid') 
  
  mod1 <-inla(form1, #the formula
              data=inla.stack.data(stk_est,spde=spde),  #the data stack
              family= 'nbinomial',   #which family the data comes from
              control.predictor=list(A=inla.stack.A(stk_est),compute=TRUE),  #compute gives you the marginals of the linear predictor
              control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
              verbose = FALSE)
  
  ind1 <-inla.stack.index(stk_est, "est_est")$data #--estimation indices 
  fit1<- exp(mod1$summary.linear.predictor[ind0,"mean"]) #--extract the backtransformed pop_hat
  fit1U <- exp(mod1$summary.linear.predictor[ind0,"0.975quant"]) #--extract the backtransformed pop_hat upper
  fit1L <- exp(mod1$summary.linear.predictor[ind0,"0.025quant"]) #--extract the backtransformed pop_hat lower
  
  fit11<- fit1
  fit11U <- fit1U ##---uncertainty based on simulated posterior is required
  fit11L <- fit1L
  (sum_fit <- sum(fit1))
  (POP1 <- cbind(ddat$pop, ddat$popm,fit11))
  apply(POP0, 2, sum, na.rm=T)
  plot(ddat$pop, fit11, col=c("red", "blue"))
  abline(a=0, b=1)
  cor(ddat$pop, fit11)
  ####
  
  
  ###---------------------------------------------------------------------------------------
  #--run posterior simulation----------------------------------------------------------
  #--extract posterior random effects
  ddat$prov_ranef <- prov_ranef(ddat, mod1)
  system.time(str(sim.pop1<- simPop(mod1,ddat, A, run)))
  
  sum(ddat$pop_hat1 <- fit11)
  dim(ddat)
  names(ddat)
  data.all1 <- ddat[,c("prov2_ID" ,"pop_hat1")]
  data.all1 <- data.frame(cbind(data.all1, sim.pop1$pop_hat))#---Contains all the simulated posterior matrix
  dim(data.all1);names(data.all1)
  
  ##--national estimates
  (national1 <- nat_total(data.all1, thinn))
  
  write.csv(national1, file="national_estimates_nb.csv")
  
  ##---cu-level
  cu.est1 <- cu_est(sim.pop1$est_data, data.all1,thinn)
  head(cu.est1)
  write.csv(cu.est1, file="CU_estimates_nb.csv")
  
  png("cu_scatterplot.png")
  plot(cu.est1$mean, cu.est1$mean)
  dev.off()
  ##-Metrics
  mod_metrics1 <- model_metrics(cu.est1$pop, cu.est1$mean, 
                                cu.est1$upper, cu.est1$lower)
  write.csv(mod_metrics1, file="fit_metrics_nb.csv")
  
  
  
  #t_map1 <- post_maps(cu.est1, shp)
  #tmap_save(t_map1, paste0(result_path2,'/posterior_cu_maps_nb.tiff'), dpi = 300, units="in")
  
    
  }







