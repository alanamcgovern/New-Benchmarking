library(INLA)
library(SUMMER)
library(tidyverse)

# load data ------
#data
load("/Users/alanamcgovern/Desktop/UN_Estimates/UN-Subnational-Estimates/Data/Malawi/Malawi_cluster_dat.rda")
dat <- mod.dat[mod.dat$age==0 & mod.dat$survey==2015 & mod.dat$years==2010,]
dat$years <- 2010
dat$died <- dat$Y
dat$years <- as.numeric(as.character(dat$years))
dat$v005 <- dat$v005/1e6
admin_key <- dat %>% dplyr::select(admin1.char,admin2.char) %>% unique()
admin_key_num <- dat %>% dplyr::select(admin1,admin2) %>% unique()

#adjacency matrix
load("/Users/alanamcgovern/Desktop/UN_Estimates/UN-Subnational-Estimates/Data/Malawi/shapeFiles_gadm/Malawi_Amat.rda")
# make separate islands for each admin1 area
admin2.mat.nested <- matrix(0,nrow(admin2.mat),nrow(admin2.mat))
for(area in 1:nrow(admin1.mat)){
  adm2_areas_tmp <- sort(admin_key_num[admin_key_num$admin1==area,]$admin2)
  for(area2 in adm2_areas_tmp){
    if(sum(admin2.mat[area2,adm2_areas_tmp])==0){
      admin2.mat.nested[area2,adm2_areas_tmp] <- admin2.mat[area2,adm2_areas_tmp]
      }else{
       admin2.mat.nested[area2,adm2_areas_tmp] <- admin2.mat[area2,adm2_areas_tmp]/sum(admin2.mat[area2,adm2_areas_tmp])}
  }
}

# get weights for aggregating from admin2 to admin 1
load("/Users/alanamcgovern/Desktop/UN_Estimates/UN-Subnational-Estimates/Data/Malawi/worldpop/adm1_weights_u1.rda")
adm1_weights <- (weight.adm1.u1 %>% filter(years==2010))[,1:2] #not that different across period of years
load("/Users/alanamcgovern/Desktop/UN_Estimates/UN-Subnational-Estimates/Data/Malawi/worldpop/adm2_weights_u1.rda")
adm2_weights <- (weight.adm2.u1 %>% filter(years==2010))[,1:2] #not that different across period of years
adm2_weights <- merge(adm2_weights,admin_key,by.x='region',by.y='admin2.char')
adm1_props <- adm2_weights %>% group_by(admin1.char) %>% summarise(adm1_prop = sum(proportion))
adm2_to_adm1_weights <- left_join(adm2_weights,adm1_props,by='admin1.char')
adm2_to_adm1_weights$adm2_prop <- adm2_to_adm1_weights$proportion/adm2_to_adm1_weights$adm1_prop
adm2_to_adm1_weights <- adm2_to_adm1_weights[,c(1,3,5)]

# get urban/rural weights for admin1 and admin2
adm1_UR_weights <- readRDS("/Users/alanamcgovern/Desktop/UN_Estimates/UN-Subnational-Estimates/Results/Malawi/UR/U1_fraction/admin1_u1_urban_weights.rds")
adm1_UR_weights <- (adm1_UR_weights %>% filter(years==2010))[,c(1,3:4)] #not that different across period of years
adm1_UR_weights <- gather(adm1_UR_weights,strata,proportion,urban:rural)
adm2_UR_weights <- readRDS("/Users/alanamcgovern/Desktop/UN_Estimates/UN-Subnational-Estimates/Results/Malawi/UR/U1_fraction/admin2_u1_urban_weights.rds")
adm2_UR_weights <- (adm2_UR_weights %>% filter(years==2010))[,c(1,3:4)] #not that different across period of years
adm2_UR_weights <- gather(adm2_UR_weights,strata,proportion,urban:rural)

# fit models -----

# admin 1 level models ---------
#direct estimates and smoothed
adm1_dir_est <- getDirect(births = dat, years = unique(dat$years), regionVar = 'admin1.char',
                          timeVar = 'years', clusterVar = '~cluster', Ntrials = 'total')

natl_dir_est <- adm1_dir_est[adm1_dir_est$region=='All',]
adm1_dir_est <- adm1_dir_est[adm1_dir_est$region!='All',]

adm1_sd_fit <- smoothDirect(adm1_dir_est, Amat = admin1.mat, time.model = NULL)

adm1_dir_est <- adm1_dir_est %>% rename(logit.var.est = var.est)
adm1_dir_est$var.est <- (adm1_dir_est$upper-adm1_dir_est$lower)/(2*1.96)
adm1_sd_est <- getSmoothed(adm1_sd_fit, Amat = Q_admin1, CI=0.95, control.compute=list(cpo=T))
adm1_sd_est$mean <- (adm1_sd_est$upper - adm1_sd_est$lower)/2 + adm1_sd_est$lower
adm1_sd_est <- data.frame(model='SD',region = adm1_sd_est$region, 
                          mean = adm1_sd_est$mean, lower95 = adm1_sd_est$lower, upper95 = adm1_sd_est$upper)

#mod1: just admin1 (fixed effect) -- may not have enough data if only 1 year of data is included/many admin1 areas
inla_adm1fe <- INLA::inla(Y ~ urban + factor(admin1) -1,
                          data=dat, family='nbinomial', E=total,
                          control.predictor = list(link = 1),
                          control.compute = list(cpo=T))
mod1_summary <- cbind(data.frame(dat %>% dplyr::select(admin1.char,urban) %>% unique()),unique(inla_adm1fe$summary.fitted.values)[,c(1,3,5)])
mod1_summary$model <- 'mod1'

#mod2: just admin1 (random effect)
hyper.inla <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                   phi = list(prior = "pc",  param = c(0.5, 2/3)))

inla_adm1re <- INLA::inla(Y ~ urban -1 + f(admin1, graph = admin1.mat, 
                                           model = "bym2",hyper=hyper.inla, 
                                           scale.model = T, adjust.for.con.comp = T), 
                        data=dat, family='nbinomial', E=total,
                        control.predictor = list(link = 1),
                        control.compute = list(cpo=T))
mod2_summary <- cbind(data.frame(dat %>% dplyr::select(admin1.char,urban) %>% unique()),unique(inla_adm1re$summary.fitted.values)[,c(1,3,5)])
mod2_summary$model <- 'mod2'

inla_adm1_ests <- rbind(mod1_summary,mod2_summary)
inla_adm1_ests <- mod2_summary
inla_adm1_ests <- merge(adm1_UR_weights,inla_adm1_ests,by.x=c('region','strata'),by.y=c('admin1.char','urban'))
inla_adm1_ests <- inla_adm1_ests %>% group_by(model,region) %>% summarise(mean =(sum(mean*proportion)), 
                                                                          lower95=(sum(`0.025quant`*proportion)),
                                                                          upper95 = (sum(`0.975quant`*proportion)))

# admin 2 level models ---------
# mod3a: both admin1 (fixed effect) and admin2 
# inla_bothfe <- INLA::inla(Y ~ urban + factor(admin1) -1 +
#                             f(admin2, graph = abs(admin2.mat),model = "bym2",hyper=hyper.inla, 
#                               scale.model = T, adjust.for.con.comp = T),
#                           data=dat, family='nbinomial', E=total,
#                           control.predictor = list(compute = FALSE, link = 1),
#                           control.compute = list(cpo=T))
# mod3a_summary <- cbind(data.frame(dat %>% dplyr::select(admin1.char,admin2.char,urban) %>% unique()),unique(inla_bothfe$summary.fitted.values)[,c(1,3,5)])
# mod3a_summary$model <- 'mod3a'

# mod3: both admin1 (fixed effect) and admin2 (nested)
inla_bothfe <- INLA::inla(Y ~ urban + factor(admin1) -1 +
                            f(admin2, graph = (admin2.mat.nested),model = "bym2",hyper=hyper.inla, 
                              scale.model = T, adjust.for.con.comp = T),
                          data=dat, family='nbinomial', E=total,
                          control.predictor = list(compute = F, link = 1),
                          control.compute = list(cpo=T))
mod3_summary <- cbind(data.frame(dat %>% dplyr::select(admin1.char,admin2.char,urban) %>% unique()),unique(inla_bothfe$summary.fitted.values)[,c(1,3,5)])
mod3_summary$model <- 'mod3'

# # mod4a: both admin1 (random effect) and admin2
# inla_bothre <- INLA::inla(Y ~ urban + -1 +
#                             f(admin1, graph = abs(admin1.mat),model = "bym2",hyper=hyper.inla, scale.model = T, adjust.for.con.comp = T) +
#                             f(admin2, graph = abs(admin2.mat),model = "bym2",hyper=hyper.inla, scale.model = T, adjust.for.con.comp = T),
#                           data=dat, family='nbinomial', E=total,
#                           control.predictor = list(compute = FALSE, link = 1))
# mod4a_summary <- cbind(data.frame(dat %>% dplyr::select(admin1.char,admin2.char,urban) %>% unique()),unique(inla_bothre$summary.fitted.values)[,c(1,3,5)])
# mod4a_summary$model <- 'mod4a'

# mod4: both admin1 (random effect) and admin2 (nested)
inla_bothre <- INLA::inla(Y ~ urban + -1 +
                            f(admin1, graph = (admin1.mat),model = "bym2",hyper=hyper.inla, scale.model = T, adjust.for.con.comp = T) +
                            f(admin2, graph = (admin2.mat.nested),model = "bym2",hyper=hyper.inla, scale.model = T, adjust.for.con.comp = T),
                          data=dat, family='nbinomial', E=total,
                          control.predictor = list(compute = FALSE, link = 1),
                          control.compute = list(cpo=T))
mod4_summary <- cbind(data.frame(dat %>% dplyr::select(admin1.char,admin2.char,urban) %>% unique()),unique(inla_bothre$summary.fitted.values)[,c(1,3,5)])
mod4_summary$model <- 'mod4'

# mod5a: just admin 2 
# inla_adm2re <- INLA::inla(Y ~ urban - 1 + f(admin2, graph = abs(admin2.mat),model = "bym2",hyper=hyper.inla, scale.model = T, adjust.for.con.comp = T),
#                           data=dat, family='nbinomial', E=total,
#                           control.predictor = list(compute = FALSE, link = 1))
# mod5a_summary <- cbind(data.frame(dat %>% dplyr::select(admin1.char,admin2.char,urban) %>% unique()),unique(inla_adm2re$summary.fitted.values)[,c(1,3,5)])
# mod5a_summary$model <- 'mod5a'

# mod5: just admin 2 (nested)
inla_adm2re <- INLA::inla(Y ~ urban - 1 + f(admin2, graph = (admin2.mat.nested),model = "bym2",hyper=hyper.inla, scale.model = T, adjust.for.con.comp = T),
                          data=dat, family='nbinomial', E=total,
                          control.predictor = list(compute = FALSE, link = 1),
                          control.compute = list(cpo=T))
mod5_summary <- cbind(data.frame(dat %>% dplyr::select(admin1.char,admin2.char,urban) %>% unique()),unique(inla_adm2re$summary.fitted.values)[,c(1,3,5)])
mod5_summary$model <- 'mod5'

inla_adm2_ests <- rbind(mod3_summary,mod4_summary,mod5_summary)
#merge from U/R to admin2
inla_adm2_ests <- merge(adm2_UR_weights,inla_adm2_ests,by.x=c('region','strata'),by.y=c('admin2.char','urban'))
inla_adm2_ests <- inla_adm2_ests %>% group_by(model,region) %>% summarise(mean = (sum(mean*proportion)), 
                                                                          lower95=(sum(`0.025quant`*proportion)),
                                                                          upper95 = (sum(`0.975quant`*proportion)))
#merge from admin2 to admin1
inla_adm2_ests <- merge(adm2_to_adm1_weights,inla_adm2_ests,by='region')
inla_adm2_to_adm1_ests <- inla_adm2_ests %>% group_by(model,admin1.char) %>% summarise(mean = (sum(mean*adm2_prop)), 
                                                                                       lower95=(sum(lower95*adm2_prop)),
                                                                                       upper95 = (sum(upper95*adm2_prop))) %>% rename(region=admin1.char)

# combine results for admin ------
adm1_ests_inla_all <- rbind(adm1_sd_est,inla_adm1_ests,inla_adm2_to_adm1_ests)

# aggregate to national ------
adm1_ests <- merge(adm1_ests_inla_all,adm1_weights,by='region')
natl_ests <- cbind(adm1_ests %>% group_by(model) %>% summarise(mean=sum(mean*proportion),lower95=sum(lower95*proportion), upper95=sum(upper95*proportion)))

# compare CV diagnostics --------
## first check that they are reliable
sum(inla_adm1fe$cpo$failure + inla_adm1re$cpo$failure + 
      inla_bothfe$cpo$failure + inla_bothre$cpo$failure + inla_adm2re$cpo$failure)==0
sum(inla_adm1re$cpo$failure) + 
      sum(inla_bothfe$cpo$failure) + sum(inla_bothre$cpo$failure) + sum(inla_adm2re$cpo$failure)==0

# very similar, but are in agreement with what we see from comparing aggregated national estimates
cpo_results <- data.frame(model=c('Adm1 FE','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2'),
           CPO=c(-mean(log(inla_adm1fe$cpo$cpo)), -mean(log(inla_adm1re$cpo$cpo)),
                 -mean(log(inla_bothfe$cpo$cpo)), -mean(log(inla_bothre$cpo$cpo)),
                 -mean(log(inla_adm2re$cpo$cpo))))

cpo_results <- data.frame(model=c('Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2'),
                          CPO=c( -mean(log(inla_adm1re$cpo$cpo)),
                                -mean(log(inla_bothfe$cpo$cpo)), -mean(log(inla_bothre$cpo$cpo)),
                                -mean(log(inla_adm2re$cpo$cpo))))

pdf('Nested Model Approach/Malawi 2010 survey=2015 230814.pdf')
{
  plot(seq(-4,-2.8,.01),dnorm(seq(-4,-2.8,.01),natl_dir_est$logit.est,sqrt(1/natl_dir_est$logit.prec)),type='l',main='National Estimates')
  for(k in 1:nrow(natl_ests)){
    abline(v=logit(natl_ests$mean[k]),col=k+1)
  }
  legend('topright',legend=c('SD','Adm1 FE','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2'),
  #legend('topright',legend=c('SD','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2'),
         col = 2:(nrow(natl_ests)+1),lwd=1.5)
  
  g1 <- ggplot() + geom_point(data=adm1_ests_inla_all, aes(x=region,y=logit(mean),col=model),pch=3) +
    geom_point(data=adm1_dir_est,aes(x=region,y=logit.est),pch=2) +
    scale_color_discrete(labels=c('Adm1 FE','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2','SD')) +
    #scale_color_discrete(labels=c('Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2','SD')) +
    ggtitle('Admin1 Estimates')
  print(g1)
  
  #plot distribution of direct estimates
  for(area in 1:nrow(admin1.mat)){
    x_seq <- seq(adm1_dir_est$logit.est[area] - 3*sqrt(adm1_dir_est$logit.var.est[area]),
                 adm1_dir_est$logit.est[area] + 3*sqrt(adm1_dir_est$logit.var.est[area]), 0.01)
    plot(x_seq,dnorm(x_seq,adm1_dir_est$logit.est[area],sqrt(adm1_dir_est$logit.var.est[area])),type='l',
         main=paste0('Distribution of direct estimate for', adm1_dir_est$region[area]),ylab='',xlab='')
    for(k in 1:length(unique(adm1_ests_inla_all$model))){
      dat_tmp <- adm1_ests_inla_all[adm1_ests_inla_all$model == unique(adm1_ests_inla_all$model)[k],]
      abline(v=logit(dat_tmp[area,]$mean),col=k+1)
    }
    abline(v=natl_dir_est$logit.est,lty=2)
    legend('topright',legend=c('SD','Adm1 FE','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2',"Natl direct"),
    #legend('topright',legend=c('SD','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2',"Natl direct"),
           col = c(2:(length(unique(adm1_ests_inla_all$model))+1),1),lwd=1.5,
           lty=c(rep(1,length(unique(adm1_ests_inla_all$model))),2))
    Sys.sleep(0.5)
  }

}
dev.off()
