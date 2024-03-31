library(tidyverse)
library(locfit)
library(rdhs)
library(surveyPrev)
library(SUMMER)
library(labelled)
library(sf)
#library(ggpubr)

load(file ='/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Info/Nigeria_general_info.Rdata') # load the country info

# load shapefiles -----
poly.path <- "/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_NGA_shp"
poly.adm1 <- st_read(dsn = poly.path, layer = poly.layer.adm1, options = "ENCODING=UTF-8")
  
geo <- getDHSgeo(country = 'Nigeria', year = 2018)
  
admin1.info <- adminInfo(geo = poly.adm1, admin = 1)
admin1.info$admin.info$admin1 <- 1:nrow(admin1.info$admin.info)
admin1.info$admin.info$admin1.char <- paste0('admin1_',admin1.info$admin.info$admin1)
row.names(admin1.info$admin.mat) <- colnames(admin1.info$admin.mat) <- admin1.info$admin.info$admin1
# re-structure admin1 matrix for INLA
admin1.mat <- matrix(0,ncol=nrow(admin1.info$admin.mat),nrow = nrow(admin1.info$admin.mat))
for(i in 1:nrow(admin1.mat)){
  admin1.mat[i,] <- -1*I(admin1.info$admin.mat[i,]!=0)
}
diag(admin1.mat) <- -colSums(admin1.mat)
save(admin1.mat, file='/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Nigeria/Nigeria_Amat.rda')

# load data ------
load(paste0('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Nigeria/Nigeria_cluster_dat.rda')) 
nmr.dat[nmr.dat$admin1.name=='Lagos',]$urban <- 'urban'
nmr.dat[nmr.dat$admin1.name=='Lagos',]$strata <- 'Lagos south west.urban'

adm1_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_nga/admin1_u1_2006_urban_frac.rds"))
adm1_UR_weights <- adm1_UR_weights[order(as.numeric(str_remove(adm1_UR_weights$adm_idx,'admin1_'))),]
#adm1_UR_weights$dhs_frac <- c(20.46,26.15,4.21,83.70,14.39,28.61,10.89,34.66,14.27,51.42,86.86,58.81,80.25,72.64,70.97,23.51,52.66,11.07,
#                              47.31,45.37,20.19,16.68,36.76,69.76,100,22.84,25.80,50.57,47.8,76.71,71.66,28.64,51.63,21.49,16.17,21.89,18.46)/100
#saveRDS(adm1_UR_weights,file = ("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_nga/admin1_u1_2017_urban_frac.rds"))

# admin1 weights 
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Nigeria/worldpop/adm1_weights_u1.rda"))
adm1_weights <- weight.adm1.u1[weight.adm1.u1$years==2016,]

# direct estimates -----
data_for_direct <- nmr.dat
data_for_direct$age <- 0
data_for_direct$years <- 'All'
data_for_direct$died <- data_for_direct$Y
dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                             regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                             Ntrials = "total",weightsVar = "v005")

natl.dir <- dir.est[dir.est$region=='All',]
admin1.dir <- dir.est[dir.est$region!='All',]
admin1.dir <- admin1.dir %>% dplyr::select(region,mean,lower,upper) %>% rename(admin1_mean = mean, admin1_lower = lower,admin1_upper=upper) %>% rename(admin1.char = region)
admin1.dir$admin1 <- as.numeric(sapply(1:nrow(admin1.dir),function(i){str_split(admin1.dir$admin1.char[i],'_')[[1]][2]}))

# admin 1 random effect model-----
admin1RE.fit <- INLA::inla(Y ~ urban -1 +  f(admin1, model="bym2", graph = admin1.mat, scale.model=T, constr=T,
                                             adjust.for.con.comp=T, 
                                             hyper=list(phi=list(prior="pc", param=c(0.5, 0.5)),
                                                        prec=list(prior='pc.prec', param=c(1, 0.01)))),
                           family = 'binomial', 
                           data = nmr.dat, 
                           control.family=list(variant=1,link='log'),
                           control.predictor = list(compute = FALSE, link = 1), 
                           control.compute = list(config=T),
                           Ntrials = nmr.dat$total)

# sample from posterior
cs <- admin1RE.fit$misc$configs$contents$tag
cs <- cs[cs != "Predictor"]
select <- list()
for (i in 1:length(cs)) {
  select[[i]] <- 0
  names(select)[i] <- cs[i]
}
sampFull <- INLA::inla.posterior.sample(n = 1000, result = admin1RE.fit, intern = TRUE, selection = select)
sampFull.draws <- matrix(NA,1000,length(sampFull[[1]]$latent))
for(i in 1:1000){
  sampFull.draws[i,] <- sampFull[[i]]$latent
}
fields <- colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)

admin1RE.inla <- matrix(0,ncol=(2*max(nmr.dat$admin1)+2),nrow=max(nmr.dat$admin1))
for(area in 1:max(nmr.dat$admin1)){
  admin1RE.inla[area,area] <- 1 #choose admin1 area
}
admin1RE.inla[,2*max(nmr.dat$admin1)+1] <- adm1_UR_weights$urb_frac
admin1RE.inla[,2*max(nmr.dat$admin1)+2] <- 1 - adm1_UR_weights$urb_frac

admin1RE.inla <- admin1RE.inla%*%t(sampFull.draws)
res.admin1RE <- data.frame(admin1 = unique(nmr.dat$admin1),admin1.char = unique(nmr.dat$admin1.char))

res.admin1RE$adm1_median <- (apply(admin1RE.inla,1,function(x){median(exp(x))}))
res.admin1RE$adm1_lower95 <- (apply(admin1RE.inla,1,function(x){quantile(exp(x),0.025)}))
res.admin1RE$adm1_upper95 <- (apply(admin1RE.inla,1,function(x){quantile(exp(x),0.975)}))

natl.inla <- rep(0,2*max(nmr.dat$admin1)+2)
for(area in 1:max(nmr.dat$admin1)){
  natl.inla[area] <- weight.adm1.u1$proportion[area] #choose admin1 area
}
natl.inla[2*max(nmr.dat$admin1)+1] <- sum(adm1_UR_weights$urb_frac)/nrow(adm1_UR_weights)
natl.inla[2*max(nmr.dat$admin1)+2] <- 1 - sum(adm1_UR_weights$urb_frac)/nrow(adm1_UR_weights)
natl.inla <- natl.inla%*%t(sampFull.draws)

res.natl <- data.frame(median = median(exp(natl.inla)),
                       lower = quantile(exp(natl.inla),0.025),
                       upper = quantile(exp(natl.inla),0.975))

#plot -----
plot(res.admin1RE$adm1_median)
abline(h=natl.dir$mean,col='red')
abline(h=natl.dir$lower,col='red',lty=2)
abline(h=natl.dir$upper,col='red',lty=2)
abline(h=res.natl$median,col='blue')
abline(h=res.natl$lower,col='blue',lty=2)
abline(h=res.natl$upper,col='blue',lty=2)
