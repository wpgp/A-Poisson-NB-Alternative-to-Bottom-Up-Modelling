library(INLA); library(raster); library(maptools)
library(gtools); library(sp); library(spdep); library(rgdal)
library(fields); library(mvtnorm);  library(geoR)
library(actuar);library(viridisLite);require(grid);require(gridExtra)
require(lattice);require(tidyverse);require(MASS);library(tmap)
library(tmaptools);library(sf)
#install.packages("cartography")

library(cartography) # mapping dedicated package
#install.packages("OpenStreetMap")
library(OpenStreetMap)

#This script if for 9 areas

path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP000008_UNFPA_PNG/Working/Chris/reviews/sim_study/option2"
out_path <- paste0(path)
#file_path <- paste0(path, "/file")
setwd(out_path)

pop.cover <- c(100, 80, 60, 40, 20)
#bld.cover <- c(100, 95, 90, 85, 80, 75, 70, 65)
metric <- c("mae", "rmse", "bias", "corr")
method <- c("onestep", "twostep")

n.pop.cover <- length(pop.cover)
#n.bld.cover <- length(bld.cover)
n.metric <- length(metric)
n.method <- length(method)

##---build the dataframe for the metrics
n.size <- n.pop.cover*n.metric*n.method
dim(dat.met <- data.frame(expand.grid(method=method,
                                      pop_cover=pop.cover,
                                      metric=metric)))

###

#paste0(out_path, "/outputs_for_80%_pop_count/95%_bldg_count")
dat_met0 <- matrix(nrow=n.pop.cover, ncol=6)
dat_met1 <- matrix(nrow=n.pop.cover, ncol=6)
for(j in 1:n.pop.cover)
{
  j-1
  pathp <- paste0(out_path,"/outputsoutputs_for_", pop.cover[j],"%","_pop_count")
    met0 <- read.csv(paste0(pathp, "/fit_metrics_pois.csv"))
    met1 <- read.csv(paste0(pathp, "/fit_metrics_nb.csv"))
    met0 <- c(pop.cover[j], unlist(met0))
    met1 <- c(pop.cover[j], unlist(met1))
    dat_met0[j,] = met0
    dat_met1[j,] = met1
  }
metrics <- as.data.frame(rbind(dat_met0,dat_met1))

##
names(metrics) <- c("pop_cover", "x", "mae","rmse", "bias", "corr") #-rename columns

var2inc <- c("pop_cover", "mae","rmse", "bias", "corr") #--gets rid of X
metrics <- metrics[,var2inc]
metrics$method <- rep(c("pois", "nb"),each=nrow(metrics)/2)#--add 'method' col

head(metrics)
write.csv(metrics, "combined_fit_metrics.csv", row.names=FALSE)

# Convert to long format for plotting 
require(reshape2)
require(ggpubr)
dim(met_long <- melt(metrics, id.vars=c("pop_cover", "method"),
                     value.name="estimate", variable.name = "metric"))

met_long$method = factor(met_long$method)
met_long$pop_cover = factor(met_long$pop_cover)
head(met_long)
table(met_long$metric)
write.csv(met_long, "combined_fit_metrics_long.csv", row.names=FALSE)

##---Plots
variable_names <- list(
  "mae" = "Mean  Absolute \n Error (MAE)" ,
  "rmse" = "Root Mean Square \n Error (RMSE)",
  "bias" = "Absolute Bias",
  "corr" = "Correlation \n Coefficient"
)

levels(met_long$method) <- c("Poisson","NBinomial")
grp <- levels(met_long$method)

variable_labeller2 <- function(variable,value){
  if (variable=='metric') {
    return(variable_names[value])
  } else {
    return(grp)
  }
}


plot_metrics <- ggplot(met_long, aes(x=pop_cover, y=estimate))+
  geom_point(aes(colour=method, shape=method), size=2)+
  geom_line()+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
  facet_wrap(~metric, scales="free",
             labeller= variable_labeller2)
plot_metrics

plot_met <- ggpar(plot_metrics, xlab="Survey Coverage (%)", ylab="Estimate",
                  legend = "right", legend.title = "Modelling Method",
                  palette = c("dark blue", "orange"),
                  #"#00AFBB", "#E7B800"
                  #colour = "bld_cover",
                  #shape= "bld_cover",
                  font.label = list(size = 15, face = "bold", color ="red"),
                  font.x = c(16),
                  font.y = c(16),
                  xtickslab.rt = 45, ytickslab.rt = 45)

plot_met
##


####============------------------------
#---------------Make Scatter plots for pop counts------------
dat_cu1<- list()
for(j in 1:n.pop.cover)
{
  #j=1
  pathp <- paste0(out_path,"/outputsoutputs_for_", pop.cover[j],"%","_pop_count")
    cu0 <- read.csv(paste0(pathp, "/CU_estimates_pois.csv"))
    cu1 <- read.csv(paste0(pathp, "/CU_estimates_nb.csv"))
    cu0$mean1 <- cu1$mean
    cu0$lower1 <- cu1$lower
    cu0$upper1 <- cu1$upper
    cu0$sd1 <- cu1$sd_pop_hat
    cu0$dens1 <- cu1$mean_dens_hat
    cu0$upper1 <- cu1$upper
    cu0$sd1 <- cu1$sd_pop_hat
    cu0$pop_cover <- rep(pop.cover[j], nrow(cu0))
    #cu0$method <- rep(c("onestep", "twosteps"),nrow(cu0)/2)
    dat_cu1[[j]] = cu0
} 
#dat_cu1  #- very huge and 

##--Extract 
datpop100 <- dat_cu1[[1]]
datpop80 <- dat_cu1[[2]]
datpop60 <- dat_cu1[[3]]
datpop40 <- dat_cu1[[4]]
datpop20 <- dat_cu1[[5]]

dim(datpop100)

##---Select variables to include
Var2Include <- c("lon", "lat", "prov2_ID", "set_typ2_ID",
                 "x2", "x3", "x4", "x5", "x6", "pop",
                 "popm", "prov_ranef",
                 "mean", "lower", "upper", "mean1", "lower1", "upper1",
                 "pop_cover")

##########
##----Extract data for 60% survey coverage
dim(dat_all <- rbind(datpop100[,Var2Include],
                  datpop80[,Var2Include],
                  datpop60[,Var2Include],
                  datpop40[,Var2Include],
                  datpop20[,Var2Include]))


names(dat_all)
dim(dat_pois <- dat_all[,names(dat_all)[c(1:15,19)]])
dat_pois$method <- rep("Poisson", nrow(dat_pois))
dat_pois$method  <- factor(dat_pois$method)


####
dim(dat_nb <- dat_all[,names(dat_all)[c(1:12,16:19)]])
dat_nb$method <- rep("NBinomial", nrow(dat_nb))
dat_nb$method  <- factor(dat_nb$method)

dat_nb2 <- dat_nb
dat_nb2$method <- rep("True", nrow(dat_nb2))
dat_nb2$method  <- factor(dat_nb2$method)


names(dat_nb)  <- names(dat_pois) 
names(dat_nb2) <- names(dat_pois) 
dat_nb2$mean <- dat_nb2$pop

#################################----------------------------------------------
dim(fdat <- rbind(dat_pois, dat_nb, dat_nb2))
names(fdat)
fdat$pop_cover <- factor(fdat$pop_cover)
##------------------------------------------------------------------------------------------------
##-----SCATTER PLOTS---------------------------------------------------------------
#################################----------------------------------------------

##--save
write.csv(fdat, "combined_cu_level_estimates.csv")


require(ggpubr)
#####################
###
var_names <- c(
  "100",
  "80",
  "60",
  "40" ,
  "20"
)

fdat$pop_cover2 <- factor(fdat$pop_cover, levels=var_names)


variable_names <- list(
  "100" = "100% \n Survey Coverage",
  "80" = "80% \n Survey Coverage",
  "60" = "60% \n Survey Coverage",
  "40" = "40% \n Survey Coverage",
  "20" = "20% \n Survey Coverage" 
)

grp <- levels(fdat$method)

variable_labeller <- function(variable,value){
  if (variable=='pop_cover2') {
    return(variable_names[value])
  } else {
    return(grp)
  }
}


#--------------------------plot all survey coverage and all satellite coverage props together


ps <- fdat %>%
  ggplot(aes(x=pop, y=mean))+
  geom_point()+
  geom_smooth(aes(colour=method),method="lm", se=F)+
  xlab(label="Observed Population (count)")+
  ylab(label="Predicted Population (count)")+
  theme_bw()+
  theme(strip.text = element_text(size = 12),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        legend.position = "right",
        #panel.border = element_blank(), 
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  facet_wrap(~pop_cover2, scales="free",
             labeller= variable_labeller, ncol=2)
ps

plots <- ggpar(ps,
               legend = "right", legend.title = "Modelling Method",
               palette = c("#E7B800","#00AFBB","red"),
               #"#999999", "#E69F00", "#56B4E9"
               ##[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
               #[8] "#FFC600" "#FFE200" "#FFFF00"
               #"#E5F5F9" "#99D8C9" "#FFE200"
               #colour = "bld_cover",
               #shape= "bld_cover",
               ##font.label = list(size = 12, face = "bold", color ="red"),
               #font.x = c(12),
               #font.y = c(12),
               ylim= c(0, max(dat_nb2$pop)),
               xtickslab.rt = 45)

plots

##---------------------------------------------------------------------------------

dim(fdat2 <- fdat[fdat$method!="True",])
###-----

###-----------------------------------------------------------------------------------
####============------------------------
#---------------National estimates------------
dat_nat1 <- list()
for(j in 1:n.pop.cover)
{
  pathp <- paste0(out_path,"/outputsoutputs_for_", pop.cover[j],"%","_pop_count")
    nat0 <- read.csv(paste0(pathp, "/national_estimates_pois.csv"))
    nat1 <- read.csv(paste0(pathp, "/national_estimates_nb.csv"))
    nat0$estimates1 <- nat1$estimates
    nat0$pop_cover <- rep(pop.cover[j], nrow(nat0))
    nat0$pop_cover1 <- nat0$pop_cover
  dat_nat1[[j]] = nat0
} 
#dat_nat1  #- very huge and 

##--Extract 
ndatpop100 <- dat_nat1[[1]]
ndatpop80 <- dat_nat1[[2]]
ndatpop60 <- dat_nat1[[3]]
ndatpop40 <- dat_nat1[[4]]
ndatpop20 <- dat_nat1[[5]]

class(ndatpop100)
names(ndatpop100)
##---Select variables to include
Var2Include2 <- c("estimates", "estimates1", "pop_cover", "pop_cover1")


##########
##----Extract data for 100% survey coverage
dim(ndat_all <- rbind(ndatpop100[,Var2Include2],
                   ndatpop80[,Var2Include2],
                   ndatpop60[,Var2Include2],
                   ndatpop40[,Var2Include2],
                   ndatpop20[,Var2Include2]))

ndat_all$measure <- rep(c("mean", "lower", "median", "upper"), nrow(ndat_all)/4)


names(ndat_all)
dim(ndat_pois <- ndat_all[,names(ndat_all)[c(1,3,5)]])
ndat_pois$method <- rep("Poisson", nrow(ndat_pois))
ndat_pois$method  <- factor(ndat_pois$method)

###
dim(ndat_nb <- ndat_all[,names(ndat_all)[c(2,4,5)]])
ndat_nb$method <- rep("NBinomial", nrow(ndat_nb))
ndat_nb$method  <- factor(ndat_nb$method)

####
names(ndat_nb)= names(ndat_pois)


ndat <- rbind(ndat_pois, ndat_nb)

nfdat = ndat
nfdat$pop_cover <- factor(nfdat$pop_cover)
nfdat$method <- factor(nfdat$method)

table(nfdat$pop_cover, nfdat$method)

dim(nfdat_mean <- nfdat[nfdat$measure=="mean",]) #---select only mean
dim(nfdat_lower <- nfdat[nfdat$measure=="lower",]) #---select only lower bound 2.5%
dim(nfdat_median <- nfdat[nfdat$measure=="median",]) #---select only median
dim(nfdat_upper <- nfdat[nfdat$measure=="upper",]) #---select only upper bound 97.5%

##----
names(nfdat_mean) <- c("mean", "pop_cover", "measure", "method")
nfdat_mean$lower <- nfdat_lower$estimates
nfdat_mean$upper <- nfdat_upper$estimates

library(knitr)
kable(nfdat_mean)
###---plots-----------------------
var_names <- c(
  "100",
  "80",
  "60",
  "40" ,
  "20"
)

nfdat_mean$pop_cover2 <- factor(nfdat_mean$pop_cover, levels=var_names)


levels(nfdat_mean$pop_cover2)  <- c(
  "100% \n Survey Coverage",
  "80% \n Survey Coverage",
  "60% \n Survey Coverage",
  "40% \n Survey Coverage",
  "20% \n Survey Coverage" 
)





pp <- ggplot(nfdat_mean, aes(pop_cover, mean, col=method))+ 
  geom_linerange(aes(ymin = lower, ymax = upper), size=1) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size=0.5)+
  geom_errorbar(aes(ymin = lower, ymax = upper), size=1,width = 0.5)+
  geom_line(aes(group = method), size=1) +
  geom_hline(yintercept=11625153 , linetype="dashed", color = "red")+
  theme_bw()

plotp<- ggpar(pp, xlab="Survey Coverage(%)", 
              ylab="Predicted Population(Count)",
              legend = "right", legend.title = "Modelling Method",
              palette = c("dark blue", "orange"),
              #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
              #[8] "#FFC600" "#FFE200" "#FFFF00"
              #"#E5F5F9" "#99D8C9" "#FFE200"
              #colour = "bld_cover",
              #shape= "bld_cover",
              font.label = list(size = 13, face = "bold", color ="red"),
              font.x = c(13),
              font.y = c(13),
              #ylim= c(0, max(p100m$pop)),
              ytickslab.rt = 45)

plotp


#ggsave("mtcars.pdf", width = 20, height = 20, units = "cm")

ggsave(plotp, file="national_all_data_error_bar.tiff", scale=1)

write.csv(nfdat_mean, "national_all_estimates2.csv")

nfdat_mean$per_error <- ((nfdat_mean$mean/11643074)-1)*100
nfdat_mean$per_error2 <- (nfdat_mean$mean/11643074)*100
plote <- ggplot(nfdat_mean, aes(fill=method, y=per_error, x=pop_cover)) + 
  geom_bar(position="stack", stat="identity")+
  geom_hline(yintercept=0 , linetype="dashed", color = "black")+
  theme_bw()
  #s#cale_fill_viridis(discrete=TRUE, name="") 

  plote<- ggpar(plote, xlab="Survey Coverage(%)", 
                ylab="Percentage Error",
                legend = "right", legend.title = "Modelling Method",
                palette = c("dark blue", "orange"),
                font.label = list(size = 13, face = "bold", color ="red"),
                font.x = c(13),
                font.y = c(13),
                #ylim= c(0, max(p100m$pop)),
                ytickslab.rt = 45)

plote
#############
###----OBTAIN OBSERVED TOTAL COUNT FOR EACH SURVEY COVERAGE
p100_summary <- datpop100 %>% 
  summarise(popt=sum(popm, na.rm=T))
p100_summary 

p80_summary <- datpop80  %>% 
  summarise(popt=sum(popm, na.rm=T))
p80_summary 

p60_summary <- datpop60 %>% 
  summarise(popt=sum(popm, na.rm=T))
p60_summary 

p40_summary <- datpop40  %>% 
  summarise(popt=sum(popm, na.rm=T))
p40_summary 

p20_summary <- datpop20  %>%
  summarise(popt=sum(popm, na.rm=T))
p20_summary 

##--save

#################-----MAPS---------------------------------------------------------------
###-----------------------------------------------------------------------------
#####---------------------------------------------------------------------------------
# for loading our data
library(raster);library(readr);library(readxl);library(sf)
# for datasets
library(maps);library(spData)
# for creating animations
library(magick)
# for plotting
library(grid);library(tmap);library(viridis); library(tmaptools)

##=================================================================================================
###-----Simulate the coordinates based on the lon and lat of PNG
shp <- readOGR(dsn = paste0(path), layer = "PNG_CU_32100_B") #combined shapefile 
plot(shp)



#post_maps <- function(cu.est, shp)
#{

shp.cu <- shp
shp.cu$obs1 <- datpop100$popm
shp.cu$mean10 <- datpop100$mean
shp.cu$mean11<- datpop100$mean1
shp.cu$uncert10 <- (datpop100$upper-datpop100$lower)/datpop100$mean
shp.cu$uncert11 <- (datpop100$upper1-datpop100$lower1)/datpop100$mean


shp.cu$obs2 <- datpop80$popm
shp.cu$mean20 <- datpop80$mean
shp.cu$mean21<- datpop80$mean1
shp.cu$uncert20 <- (datpop80$upper-datpop80$lower)/datpop80$mean
shp.cu$uncert21 <- (datpop80$upper1-datpop80$lower1)/datpop80$mean

shp.cu$obs3 <- datpop60$popm
shp.cu$mean30 <- datpop60$mean
shp.cu$mean31<- datpop60$mean1
shp.cu$uncert30 <- (datpop60$upper-datpop60$lower)/datpop60$mean
shp.cu$uncert31 <- (datpop60$upper1-datpop60$lower1)/datpop60$mean

shp.cu$obs4 <- datpop40$popm
shp.cu$mean40 <- datpop40$mean
shp.cu$mean41<- datpop40$mean1
shp.cu$uncert40 <- (datpop40$upper-datpop40$lower)/datpop40$mean
shp.cu$uncert41 <- (datpop40$upper1-datpop40$lower1)/datpop40$mean


shp.cu$obs5 <- datpop20$popm
shp.cu$mean50 <- datpop20$mean
shp.cu$mean51<- datpop20$mean1
shp.cu$uncert50 <- (datpop20$upper-datpop20$lower)/datpop20$mean
shp.cu$uncert51 <- (datpop20$upper1-datpop20$lower1)/datpop20$mean



###
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM = spTransform(shp.cu, crs.UTM)

names(shpcu_UTM)
#plot(shpcu_UTM)
#spplot(shpcu_UTM, "mean")



tmap_options(check.and.fix = TRUE)
tmap_mode("plot")


#
#max(shp100 <- c(shpcu_UTM$obs1, shpcu_UTM$obs2, shpcu_UTM$obs3, shpcu_UTM$obs4, shpcu_UTM$obs5,
#               shpcu_UTM$pred10, shpcu_UTM$pred20, shpcu_UTM$pred30, shpcu_UTM$pred40, shpcu_UTM$pred50,
#                shpcu_UTM$pred11, shpcu_UTM$pred21, shpcu_UTM$pred31, shpcu_UTM$pred41, shpcu_UTM$pred51), na.rm=T)


#max(shp100u <- c(shpcu_UTM$uncert10, shpcu_UTM$uncert20, shpcu_UTM$uncert30, shpcu_UTM$uncert40, shpcu_UTM$uncert50,
#               shpcu_UTM$uncert11, shpcu_UTM$uncert21, shpcu_UTM$uncert31, shpcu_UTM$uncert41, shpcu_UTM$uncert51), na.rm=T)

##=============================================================
#--Observed counts
obs1 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="obs1", title="Simulated Count",
              #style="fixed", palette="viridis")+
              style="fixed", 
              #as.count=T,
              #n=6,
              #legend.hist=T,
              palette=turbo(n=100), 
              legend.show = T#,
              #breaks=c(1,50,150,300,450, 700,1000,3000)
              )+
  tm_layout(legend.outside =T, legend.text.size=1.5, legend.title.size=1.5)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1, breaks=c(0, 200, 400, 600))+
  #tm_layout(main.title = "Simulated Counts at 100% Survey Coverage") 
  tm_layout(main.title = "") 

##--predicted counts 
pred10 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean10", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = T, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert10 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert10", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred11 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean11", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert11<-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert11", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------
#--Observed counts
obs2 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="obs2", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred20 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean20", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert20 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert20", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred21 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean21", title="Predicted Count",
              style="fixed", legend.hist=F,, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert21 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert21", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--------------------------------------------------------------------------------------------------------------
#--Observed counts
obs3 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="obs3", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred30 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean30", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert30 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert30", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred31 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean31", title="Predicted Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert31 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert31", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------------------
#--Observed counts
obs4 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="obs4", title="Observed Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred40 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean40", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert40 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert40", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred41 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean41", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert41<-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert41", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##============================================================================= 
#--Observed counts
obs5 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="obs5", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred50 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean50", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert50 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert50", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred51 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean51", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert51 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert51", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



###########################################################################
#####--------------95% Satellite Observation coverage 
############################################################################
shp.cu2 <- shp
shp.cu2$obs1 <- p100_0_b2$popm
shp.cu2$mean10 <- p100_0_b2$mean
shp.cu2$mean11<- p100_1_b2$mean
shp.cu2$uncert10 <- (p100_0_b2$upper-p100_0_b2$lower)/p100_0_b2$mean
shp.cu2$uncert11 <- (p100_1_b2$upper-p100_1_b2$lower)/p100_1_b2$mean


shp.cu2$obs2 <- p80_0_b2$popm
shp.cu2$mean20 <- p80_0_b2$mean
shp.cu2$mean21 <- p80_1_b2$mean
shp.cu2$uncert20 <- (p80_0_b2$upper-p80_0_b2$lower)/p80_0_b2$mean
shp.cu2$uncert21 <- (p80_1_b2$upper-p80_1_b2$lower)/p80_1_b2$mean


shp.cu2$obs3 <- p60_0_b2$popm
shp.cu2$mean30 <- p60_0_b2$mean
shp.cu2$mean31 <- p60_1_b2$mean
shp.cu2$uncert30 <- (p60_0_b2$upper-p60_0_b2$lower)/p60_0_b2$mean
shp.cu2$uncert31 <- (p60_1_b2$upper-p60_1_b2$lower)/p60_1_b2$mean


shp.cu2$obs4 <- p40_0_b2$popm
shp.cu2$mean40 <- p40_0_b2$mean
shp.cu2$mean41 <- p40_1_b2$mean
shp.cu2$uncert40 <- (p40_0_b2$upper-p40_0_b2$lower)/p40_0_b2$mean
shp.cu2$uncert41 <- (p40_1_b2$upper-p40_1_b2$lower)/p40_1_b2$mean


shp.cu2$obs5 <- p20_0_b2$popm
shp.cu2$mean50 <- p20_0_b2$mean
shp.cu2$mean51 <- p20_1_b2$mean
shp.cu2$uncert50 <- (p20_0_b2$upper-p20_0_b2$lower)/p20_0_b2$mean
shp.cu2$uncert51 <- (p20_1_b2$upper-p20_1_b2$lower)/p20_1_b2$mean



###
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM2 = spTransform(shp.cu2, crs.UTM)

names(shpcu_UTM2)
#plot(shpcu_UTM2)
#spplot(shpcu_UTM2, "mean")



tmap_options(check.and.fix = TRUE)
tmap_mode("plot")


#
max(shp100 <- c(shpcu_UTM2$obs1, shpcu_UTM2$obs2, shpcu_UTM2$obs3, shpcu_UTM2$obs4, shpcu_UTM2$obs5,
                shpcu_UTM2$pred10, shpcu_UTM2$pred20, shpcu_UTM2$pred30, shpcu_UTM2$pred40, shpcu_UTM2$pred50,
                shpcu_UTM2$pred11, shpcu_UTM2$pred21, shpcu_UTM2$pred31, shpcu_UTM2$pred41, shpcu_UTM2$pred51), na.rm=T)


max(shp100u <- c(shpcu_UTM2$uncert10, shpcu_UTM2$uncert20, shpcu_UTM2$uncert30, shpcu_UTM2$uncert40, shpcu_UTM2$uncert50,
                 shpcu_UTM2$uncert11, shpcu_UTM2$uncert21, shpcu_UTM2$uncert31, shpcu_UTM2$uncert41, shpcu_UTM2$uncert51), na.rm=T)

##=============================================================
#--Observed counts
obs1 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="obs1", title="Simulated Count",
              #style="fixed", palette="viridis")+
              style="fixed", 
              #as.count=T,
              #n=6,
              #legend.hist=T,
              palette="viridis", 
              legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside =T, legend.text.size=1.5, legend.title.size=1.5)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1, breaks=c(0, 200, 400, 600))+
  #tm_layout(main.title = "Simulated Counts at 100% Survey Coverage") 
  tm_layout(main.title = "") 

##--predicted counts 
pred10 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean10", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = T, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert10 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert10", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred11 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean11", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert11<-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert11", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------
#--Observed counts
obs2 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="obs2", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred20 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean20", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert20 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert20", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred21 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean21", title="Predicted Count",
              style="fixed", legend.hist=F,, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert21 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert21", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--------------------------------------------------------------------------------------------------------------
#--Observed counts
obs3 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="obs3", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred30 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean30", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert30 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert30", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred31 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean31", title="Predicted Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert31 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert31", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------------------
#--Observed counts
obs4 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="obs4", title="Observed Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred40 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean40", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert40 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert40", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred41 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean41", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert41<-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert41", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##============================================================================= 
#--Observed counts
obs5 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="obs5", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred50 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean50", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert50 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert50", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred51 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean51", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert51 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert51", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


###########################################################################
#####--------------90% Satellite Observation coverage 
############################################################################
shp.cu3 <- shp
shp.cu3$obs1 <- p100_0_b3$popm
shp.cu3$mean10 <- p100_0_b3$mean
shp.cu3$mean11<- p100_1_b3$mean
shp.cu3$uncert10 <- (p100_0_b3$upper-p100_0_b3$lower)/p100_0_b3$mean
shp.cu3$uncert11 <- (p100_1_b3$upper-p100_1_b3$lower)/p100_1_b3$mean


shp.cu3$obs2 <- p80_0_b3$popm
shp.cu3$mean20 <- p80_0_b3$mean
shp.cu3$mean21 <- p80_1_b3$mean
shp.cu3$uncert20 <- (p80_0_b3$upper-p80_0_b3$lower)/p80_0_b3$mean
shp.cu3$uncert21 <- (p80_1_b3$upper-p80_1_b3$lower)/p80_1_b3$mean


shp.cu3$obs3 <- p60_0_b3$popm
shp.cu3$mean30 <- p60_0_b3$mean
shp.cu3$mean31 <- p60_1_b3$mean
shp.cu3$uncert30 <- (p60_0_b3$upper-p60_0_b3$lower)/p60_0_b3$mean
shp.cu3$uncert31 <- (p60_1_b3$upper-p60_1_b3$lower)/p60_1_b3$mean


shp.cu3$obs4 <- p40_0_b3$popm
shp.cu3$mean40 <- p40_0_b3$mean
shp.cu3$mean41 <- p40_1_b3$mean
shp.cu3$uncert40 <- (p40_0_b3$upper-p40_0_b3$lower)/p40_0_b3$mean
shp.cu3$uncert41 <- (p40_1_b3$upper-p40_1_b3$lower)/p40_1_b3$mean


shp.cu3$obs5 <- p20_0_b3$popm
shp.cu3$mean50 <- p20_0_b3$mean
shp.cu3$mean51 <- p20_1_b3$mean
shp.cu3$uncert50 <- (p20_0_b3$upper-p20_0_b3$lower)/p20_0_b3$mean
shp.cu3$uncert51 <- (p20_1_b3$upper-p20_1_b3$lower)/p20_1_b3$mean



###
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM3 = spTransform(shp.cu3, crs.UTM)

names(shpcu_UTM3)
#plot(shpcu_UTM3)
#spplot(shpcu_UTM3, "mean")



tmap_options(check.and.fix = TRUE)
tmap_mode("plot")


#
max(shp100 <- c(shpcu_UTM3$obs1, shpcu_UTM3$obs2, shpcu_UTM3$obs3, shpcu_UTM3$obs4, shpcu_UTM3$obs5,
                shpcu_UTM3$pred10, shpcu_UTM3$pred20, shpcu_UTM3$pred30, shpcu_UTM3$pred40, shpcu_UTM3$pred50,
                shpcu_UTM3$pred11, shpcu_UTM3$pred21, shpcu_UTM3$pred31, shpcu_UTM3$pred41, shpcu_UTM3$pred51), na.rm=T)


max(shp100u <- c(shpcu_UTM3$uncert10, shpcu_UTM3$uncert20, shpcu_UTM3$uncert30, shpcu_UTM3$uncert40, shpcu_UTM3$uncert50,
                 shpcu_UTM3$uncert11, shpcu_UTM3$uncert21, shpcu_UTM3$uncert31, shpcu_UTM3$uncert41, shpcu_UTM3$uncert51), na.rm=T)

##=============================================================
#--Observed counts
obs1 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="obs1", title="Simulated Count",
              #style="fixed", palette="viridis")+
              style="fixed", 
              #as.count=T,
              #n=6,
              #legend.hist=T,
              palette="viridis", 
              legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside =T, legend.text.size=1.5, legend.title.size=1.5)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1, breaks=c(0, 200, 400, 600))+
  #tm_layout(main.title = "Simulated Counts at 100% Survey Coverage") 
  tm_layout(main.title = "") 

##--predicted counts 
pred10 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean10", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = T, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert10 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert10", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20)
  )+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred11 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean11", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert11<-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert11", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------
#--Observed counts
obs2 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="obs2", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred20 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean20", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert20 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert20", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred21 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean21", title="Predicted Count",
              style="fixed", legend.hist=F,, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert21 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert21", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--------------------------------------------------------------------------------------------------------------
#--Observed counts
obs3 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="obs3", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred30 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean30", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert30 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert30", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred31 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean31", title="Predicted Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert31 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert31", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------------------
#--Observed counts
obs4 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="obs4", title="Observed Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred40 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean40", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert40 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert40", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred41 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean41", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert41<-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert41", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##============================================================================= 
#--Observed counts
obs5 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="obs5", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred50 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean50", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert50 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert50", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred51 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean51", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert51 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert51", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 




###########################################################################
#####--------------85% Satellite Observation coverage 
############################################################################
shp.cu4 <- shp
shp.cu4$obs1 <- p100_0_b4$popm
shp.cu4$mean10 <- p100_0_b4$mean
shp.cu4$mean11<- p100_1_b4$mean
shp.cu4$uncert10 <- (p100_0_b4$upper-p100_0_b4$lower)/p100_0_b4$mean
shp.cu4$uncert11 <- (p100_1_b4$upper-p100_1_b4$lower)/p100_1_b4$mean


shp.cu4$obs2 <- p80_0_b4$popm
shp.cu4$mean20 <- p80_0_b4$mean
shp.cu4$mean21 <- p80_1_b4$mean
shp.cu4$uncert20 <- (p80_0_b4$upper-p80_0_b4$lower)/p80_0_b4$mean
shp.cu4$uncert21 <- (p80_1_b4$upper-p80_1_b4$lower)/p80_1_b4$mean


shp.cu4$obs3 <- p60_0_b4$popm
shp.cu4$mean30 <- p60_0_b4$mean
shp.cu4$mean31 <- p60_1_b4$mean
shp.cu4$uncert30 <- (p60_0_b4$upper-p60_0_b4$lower)/p60_0_b4$mean
shp.cu4$uncert31 <- (p60_1_b4$upper-p60_1_b4$lower)/p60_1_b4$mean


shp.cu4$obs4 <- p40_0_b4$popm
shp.cu4$mean40 <- p40_0_b4$mean
shp.cu4$mean41 <- p40_1_b4$mean
shp.cu4$uncert40 <- (p40_0_b4$upper-p40_0_b4$lower)/p40_0_b4$mean
shp.cu4$uncert41 <- (p40_1_b4$upper-p40_1_b4$lower)/p40_1_b4$mean


shp.cu4$obs5 <- p20_0_b4$popm
shp.cu4$mean50 <- p20_0_b4$mean
shp.cu4$mean51 <- p20_1_b4$mean
shp.cu4$uncert50 <- (p20_0_b4$upper-p20_0_b4$lower)/p20_0_b4$mean
shp.cu4$uncert51 <- (p20_1_b4$upper-p20_1_b4$lower)/p20_1_b4$mean



###
require(rgdal)
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM4 = spTransform(shp.cu4, crs.UTM)

names(shpcu_UTM4)
#plot(shpcu_UTM4)
#spplot(shpcu_UTM4, "mean")



tmap_options(check.and.fix = TRUE)
tmap_mode("plot")


#
max(shp100 <- c(shpcu_UTM4$obs1, shpcu_UTM4$obs2, shpcu_UTM4$obs3, shpcu_UTM4$obs4, shpcu_UTM4$obs5,
                shpcu_UTM4$pred10, shpcu_UTM4$pred20, shpcu_UTM4$pred30, shpcu_UTM4$pred40, shpcu_UTM4$pred50,
                shpcu_UTM4$pred11, shpcu_UTM4$pred21, shpcu_UTM4$pred31, shpcu_UTM4$pred41, shpcu_UTM4$pred51), na.rm=T)


max(shp100u <- c(shpcu_UTM4$uncert10, shpcu_UTM4$uncert20, shpcu_UTM4$uncert30, shpcu_UTM4$uncert40, shpcu_UTM4$uncert50,
                 shpcu_UTM4$uncert11, shpcu_UTM4$uncert21, shpcu_UTM4$uncert31, shpcu_UTM4$uncert41, shpcu_UTM4$uncert51), na.rm=T)

##=============================================================
#--Observed counts
obs1 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="obs1", title="Simulated Count",
              #style="fixed", palette="viridis")+
              style="fixed", 
              #as.count=T,
              #n=6,
              #legend.hist=T,
              palette="viridis", 
              legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside =T, legend.text.size=1.5, legend.title.size=1.5)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1, breaks=c(0, 200, 400, 600))+
  #tm_layout(main.title = "Simulated Counts at 100% Survey Coverage") 
  tm_layout(main.title = "") 

##--predicted counts 
pred10 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean10", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = T, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert10 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert10", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25)
  )+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred11 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean11", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert11<-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert11", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------
#--Observed counts
obs2 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="obs2", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred20 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean20", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert20 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert20", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred21 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean21", title="Predicted Count",
              style="fixed", legend.hist=F,, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert21 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert21", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--------------------------------------------------------------------------------------------------------------
#--Observed counts
obs3 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="obs3", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred30 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean30", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert30 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert30", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred31 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean31", title="Predicted Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert31 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert31", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------------------
#--Observed counts
obs4 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="obs4", title="Observed Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred40 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean40", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert40 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert40", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred41 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean41", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert41<-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert41", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##============================================================================= 
#--Observed counts
obs5 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="obs5", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred50 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean50", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert50 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert50", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred51 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean51", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert51 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert51", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 




###########################################################################
#####--------------80% Satellite Observation coverage 
############################################################################
shp.cu5 <- shp
shp.cu5$obs1 <- p100_0_b5$popm
shp.cu5$mean10 <- p100_0_b5$mean
shp.cu5$mean11<- p100_1_b5$mean
shp.cu5$uncert10 <- (p100_0_b5$upper-p100_0_b5$lower)/p100_0_b5$mean
shp.cu5$uncert11 <- (p100_1_b5$upper-p100_1_b5$lower)/p100_1_b5$mean


shp.cu5$obs2 <- p80_0_b5$popm
shp.cu5$mean20 <- p80_0_b5$mean
shp.cu5$mean21 <- p80_1_b5$mean
shp.cu5$uncert20 <- (p80_0_b5$upper-p80_0_b5$lower)/p80_0_b5$mean
shp.cu5$uncert21 <- (p80_1_b5$upper-p80_1_b5$lower)/p80_1_b5$mean


shp.cu5$obs3 <- p60_0_b5$popm
shp.cu5$mean30 <- p60_0_b5$mean
shp.cu5$mean31 <- p60_1_b5$mean
shp.cu5$uncert30 <- (p60_0_b5$upper-p60_0_b5$lower)/p60_0_b5$mean
shp.cu5$uncert31 <- (p60_1_b5$upper-p60_1_b5$lower)/p60_1_b5$mean


shp.cu5$obs4 <- p40_0_b5$popm
shp.cu5$mean40 <- p40_0_b5$mean
shp.cu5$mean41 <- p40_1_b5$mean
shp.cu5$uncert40 <- (p40_0_b5$upper-p40_0_b5$lower)/p40_0_b5$mean
shp.cu5$uncert41 <- (p40_1_b5$upper-p40_1_b5$lower)/p40_1_b5$mean


shp.cu5$obs5 <- p20_0_b5$popm
shp.cu5$mean50 <- p20_0_b5$mean
shp.cu5$mean51 <- p20_1_b5$mean
shp.cu5$uncert50 <- (p20_0_b5$upper-p20_0_b5$lower)/p20_0_b5$mean
shp.cu5$uncert51 <- (p20_1_b5$upper-p20_1_b5$lower)/p20_1_b5$mean



###
require(rgdal)
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM5 = spTransform(shp.cu5, crs.UTM)

names(shpcu_UTM5)
#plot(shpcu_UTM5)
#spplot(shpcu_UTM5, "mean")



tmap_options(check.and.fix = TRUE)
tmap_mode("plot")


#
max(shp100 <- c(shpcu_UTM5$obs1, shpcu_UTM5$obs2, shpcu_UTM5$obs3, shpcu_UTM5$obs4, shpcu_UTM5$obs5,
                shpcu_UTM5$pred10, shpcu_UTM5$pred20, shpcu_UTM5$pred30, shpcu_UTM5$pred40, shpcu_UTM5$pred50,
                shpcu_UTM5$pred11, shpcu_UTM5$pred21, shpcu_UTM5$pred31, shpcu_UTM5$pred41, shpcu_UTM5$pred51), na.rm=T)


max(shp100u <- c(shpcu_UTM5$uncert10, shpcu_UTM5$uncert20, shpcu_UTM5$uncert30, shpcu_UTM5$uncert40, shpcu_UTM5$uncert50,
                 shpcu_UTM5$uncert11, shpcu_UTM5$uncert21, shpcu_UTM5$uncert31, shpcu_UTM5$uncert41, shpcu_UTM5$uncert51), na.rm=T)

##=============================================================
#--Observed counts
obs1 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="obs1", title="Simulated Count",
              #style="fixed", palette="viridis")+
              style="fixed", 
              #as.count=T,
              #n=6,
              #legend.hist=T,
              palette="viridis", 
              legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside =T, legend.text.size=1.5, legend.title.size=1.5)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1, breaks=c(0, 200, 400, 600))+
  #tm_layout(main.title = "Simulated Counts at 100% Survey Coverage") 
  tm_layout(main.title = "") 

##--predicted counts 
pred10 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean10", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = T, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert10 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert10", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4)
  )+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred11 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean11", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert11<-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert11", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------
#--Observed counts
obs2 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="obs2", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred20 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean20", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert20 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert20", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred21 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean21", title="Predicted Count",
              style="fixed", legend.hist=F,, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert21 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert21", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--------------------------------------------------------------------------------------------------------------
#--Observed counts
obs3 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="obs3", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred30 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean30", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert30 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert30", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred31 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean31", title="Predicted Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert31 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert31", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------------------
#--Observed counts
obs4 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="obs4", title="Observed Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred40 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean40", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert40 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert40", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred41 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean41", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert41<-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert41", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##============================================================================= 
#--Observed counts
obs5 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="obs5", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred50 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean50", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert50 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert50", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred51 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean51", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert51 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert51", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



####################################################################################################
######--------------------Perecentage error---------------------------------------------------------
dim(p100mm <- p100m[p100m$method!="FullData",])
p100mm$p_error <- (p100mm$mean/p100mm$pop)*100

dim(p80mm <- p80m[p80m$method!="FullData",])
p80mm$p_error <- (p80mm$mean/p80mm$pop)*100

dim(p60mm <- p60m[p60m$method!="FullData",])
p60mm$p_error <- (p60mm$mean/p60mm$pop)*100

dim(p40mm <- p40m[p40m$method!="FullData",])
p40mm$p_error <- (p40mm$mean/p40mm$pop)*100

dim(p20mm <- p20m[p20m$method!="FullData",])
p20mm$p_error <- (p20mm$mean/p20mm$pop)*100


table(fdat$method)
dim(fdatm <- fdat[fdat$method!="FullData",])
fdatm$R_error <- (fdatm$mean/fdatm$pop)/(sum(fdatm$mean/fdatm$pop))

# plot
bar_plot <- ggplot(fdatm, aes(fill=method, y=R_error, x=pop_cover)) + 
  geom_bar(position="stack", stat="identity")+
  #geom_hline(yintercept=1 , linetype="dashed", color = "black")+
  #scale_fill_viridis(discrete=TRUE, name="") +
  theme_bw()+
  #theme_ipsum() +
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14),
        #panel.border = element_blank(), 
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  facet_wrap(~bld_cover2, ncol=2)
bar_plot<- ggpar(bar_plot, xlab="Survey Coverage(%)", ylab="Normalized Relative Error Rate",
                 legend = "right", legend.title = "Modelling Method",
                 #palette = c("#00AFBB","blue", "red"),
                 #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
                 #[8] "#FFC600" "#FFE200" "#FFFF00"
                 #"#E5F5F9" "#99D8C9" "#FFE200"
                 #colour = "bld_cover",
                 #shape= "bld_cover",
                 font.label = list(size = 15, face = "bold", color ="red"),
                 font.x = c(16),
                 font.y = c(16),
                 #ylim= c(0, max(p100m$pop)),
                 xtickslab.rt = 45, ytickslab.rt = 45)

bar_plot

ggsave(plot3a, file="all_data_scatter_3a.tiff", scale=1)




















