library(tidyverse)

load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Uganda/Uganda_cluster_dat.rda")
load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Uganda/shapeFiles/Uganda_Amat.rda")
load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Uganda/shapeFiles/Uganda_Amat_names.rda")
load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Uganda/worldpop/adm1_weights_u1.rda")
load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Uganda/worldpop/adm2_weights_u1.rda")

adm1_UR_weights <- readRDS("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Uganda/UR/admin1_u1_urban_weights_2014frame.rds")
adm1_UR_weights <- adm1_UR_weights[adm1_UR_weights$years==2012,]
adm2_UR_weights <- readRDS("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Uganda/UR/admin2_u1_urban_weights_2014frame.rds")
adm2_UR_weights <- adm2_UR_weights[adm2_UR_weights$years==2012,]
## need urban rural props

#2010-2014

years_t <- 2010:2014
nmr.dat <- mod.dat %>% filter(years %in% years_t & survey==2016 & age==0)
weight.adm1.u1 <- weight.adm1.u1[weight.adm1.u1$year==2012,]
weight.adm2.u1 <- weight.adm2.u1[weight.adm2.u1$year==2012,]

admin_key <- readxl::read_excel("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Admin1_Admin2_Key.xlsx") %>% filter(country=='Uganda') %>% dplyr::select(-country)
admin_key$A2 <- as.numeric(str_remove(admin_key$admin2.char,'admin2_'))
admin_key$A1 <- as.numeric(str_remove(admin_key$admin1.char,'admin1_'))

## direct estimates

data_for_direct <- nmr.dat
data_for_direct$age <- 0
data_for_direct$years <- 'All'
data_for_direct$died <- data_for_direct$Y

dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                             regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                             Ntrials = "total",weightsVar = "v005")

sd.fit <- SUMMER::smoothDirect(dir.est,Amat = admin1.mat,time.model = NULL)
sd.est <- getSmoothed(sd.fit)


dir.est %>% ggplot() + geom_point(aes(mean,region)) + geom_point(aes(lower,region)) + geom_point(aes(upper,region)) 
  
## INLA: admin2 with admin1 fixed effect
admin2FE.fit <- INLA::inla(Y ~ urban + admin1.char + 
                             f(admin2, model="bym2", graph = admin2.mat, scale.model=T, constr=T,
                                                          adjust.for.con.comp=T, 
                                                          hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1))) -1,
                           family = 'binomial', 
                           data = nmr.dat, 
                           control.family=list(variant=1,link='log'),
                           control.predictor = list(compute = FALSE, link = 1), 
                           control.compute = list(config=T),
                           Ntrials = total)

# sample from posterior
cs <- admin2FE.fit$misc$configs$contents$tag
cs <- cs[cs != "Predictor"]
select <- list()
for (i in 1:length(cs)) {
  select[[i]] <- 0
  names(select)[i] <- cs[i]
}
sampFull <- INLA::inla.posterior.sample(n = 1000, result = admin2FE.fit, intern = TRUE, selection = select)
sampFull.draws <- matrix(NA,1000,length(sampFull[[1]]$latent))
for(i in 1:1000){
  sampFull.draws[i,] <- sampFull[[i]]$latent
}
fields <- colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)

inla_r_samples <- matrix(0,ncol=(4+2*135+1),nrow=135)
for(area in 1:nrow(inla_r_samples)){
  inla_r_samples[area,area] <- 1 #choose admin2 area
  if(admin_key[admin_key$A2==area,]$A1>1){
    inla_r_samples[area,2*max(admin_key$A2)+1+admin_key[admin_key$A2==area,]$A1] <- 1 #choose admin1 area
  }
}
inla_r_samples[,2*max(admin_key$A2)+1] <- adm2_UR_weights$urban
inla_r_samples[,2*max(admin_key$A2)+2] <- adm2_UR_weights$rural

inla_r_samples <- inla_r_samples%*%t(sampFull.draws)
admin2.inla <- data.frame(admin2 = 1:nrow(admin_key))
admin2.inla$admin2.char <- paste0('admin2_', admin2.inla$admin2)

admin2.inla$adm2_median <- (apply(inla_r_samples,1,function(x){median(exp(x))}))
admin2.inla$adm2_lower95 <- (apply(inla_r_samples,1,function(x){quantile(exp(x),0.025)}))
admin2.inla$adm2_upper95 <- (apply(inla_r_samples,1,function(x){quantile(exp(x),0.975)}))
admin2.inla <- merge(admin2.inla,admin_key)

admin1.inla <- merge(admin2.inla,weight.adm2.u1)
admin1.inla <- admin1.inla %>% group_by(admin1.char) %>% summarise(median = sum(proportion*adm2_median)/sum(proportion))

admin2.res <- merge(admin_key,admin2.inla,by=c('admin2.char'))

dir.est %>% ggplot() + geom_point(aes(mean,region),pch=17,size=3) + geom_point(aes(lower,region),pch=3,size=3) + 
  geom_point(aes(upper,region),pch=3,size=3) + geom_point(data=admin2.inla,aes(adm2_median,admin1.char),col='blue',pch=1,size=3) +
                                               geom_point(data =admin1.inla,aes(median,admin1.char),pch=17,col='red',size=3)

merge(admin1.inla,merge(weight.adm1.u1,admin1.names,by.x = 'admin1.name', by.y='GADM'),by.x = 'admin1.char', by.y = 'Internal') %>%
  summarise(sum(proportion*median)/sum(proportion))
