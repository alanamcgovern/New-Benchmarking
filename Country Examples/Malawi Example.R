country_t <- 'Malawi'
year_range <- 2011:2015
survey_year <- 2015

## load data/organize survey for algorithm -----
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/",country_t,"_cluster_dat_1frame.rda"))
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/shapeFiles_gadm/",country_t,"_Amat.rda"))

load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/worldpop/adm2_weights_u1.rda"))
admin_key <- readxl::read_excel("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Admin1_Admin2_Key.xlsx") %>% filter(country==country_t) %>% dplyr::select(-country)
admin_key$A1 <- as.numeric(sapply(1:nrow(admin_key),function(i){str_split(admin_key$admin1.char[i],'_')[[1]][2]}))
admin_key <- admin_key[order(admin_key$A1),]
admin_key$A2 <-  as.numeric(sapply(1:nrow(admin_key),function(i){str_split(admin_key$admin2.char[i],'_')[[1]][2]}))

adm1_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_mwi/admin1_u1_2013_urban_frac.rds"))
adm1_UR_weights <- adm1_UR_weights[order(as.numeric(str_remove(adm1_UR_weights$adm_idx,'admin1_'))),]

adm2_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_mwi/admin2_u1_2013_urban_frac.rds"))
adm2_UR_weights <- adm2_UR_weights[order(as.numeric(str_remove(adm2_UR_weights$adm_idx,'admin2_'))),]

# get generalized inverse of scaled Q for BYM2
{
  Q.admin2 <- -admin2.mat
  Q.admin2 <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})*Q.admin2
  diag(Q.admin2) <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})
  diag(Q.admin2)[diag(Q.admin2)==0] <- 1
  Q_scaled <- INLA::inla.scale.model(Q.admin2, constr=list(A=t(eigen(Q.admin2)$vectors[,eigen(Q.admin2)$values<1e-10]), e=rep(0,sum(eigen(Q.admin2)$values<1e-10))))
  Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))
}

# prepare direct estimates
{
  data_for_direct <- mod.dat %>% filter(years %in% year_range, survey==survey_year,age==0)
  data_for_direct$years <- 'All'
  data_for_direct$died <- data_for_direct$Y
  dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                               regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                               Ntrials = "total",weightsVar = "v005")
  sd.fit <- smoothDirect(dir.est,Amat = admin1.mat,time.model = NULL)
  sd.est <- getSmoothed(sd.fit)
  
  natl.dir <- dir.est[dir.est$region=='All',]
  dir.est <- dir.est[dir.est$region!='All',]
  dir.est <- dir.est[order(as.numeric(sapply(1:nrow(dir.est),function(i){str_split(dir.est$region[i],'_')[[1]][2]}))),]
}

# input strata level DHS survey info 
{
  dhs_dat <- rbind(admin_key,admin_key)
  dhs_dat$U <- c(rep(1,nrow(admin_key)),rep(0,nrow(admin_key)))
  dhs_dat$A1 <- as.numeric(sapply(1:nrow(dhs_dat),function(i){str_split(dhs_dat$admin1.char[i],'_')[[1]][2]}))
  dhs_dat <- dhs_dat[order(dhs_dat$A1,dhs_dat$A2),]
   # total number of HH
    dhs_dat$M <- c(2924,34856,8574,49234,299,1721,31061,138777,2276,39993,3847,32190,4489,141389,4479,117405,8964,118301,153717,275194,3570,93639,
                   5010,57458,3306,110485,1555,45873,6089,71442,5037,70619,153578,80879,2830,95205,592,70968,5303,109833,8473,177442,3243,124174,
                   3445,18573,366,25049,4227,48373,1117,75562,2405,139634,19041,142394)
    
    # total number of clusters in strata
    dhs_dat$n_clusters <- c(11,205,37,370,4,18,122,825,12,229,12,156,15,486,18,450,29,486,458,1173,12,374,16,177,11,468,6,204,22,416,17,275,412,381,
                            16,380,2,334,19,436,25,614,17,658,9,80,3,157,14,241,3,316,12,674,79,584)
    # average number of HH per cluster
    dhs_dat$M_bar <- round(dhs_dat$M/dhs_dat$n_clusters)
    dhs_dat <- dhs_dat %>% mutate(n_HH_samp=ifelse(U==1,30,33))
}

# organize survey data
{
  alg_dat <- mod.dat %>% filter(years %in% year_range, survey==survey_year,age==0) %>% mutate(U=ifelse(urban=='urban',1,0)) %>% 
    dplyr::select(cluster,total,Y,U,admin1,admin1.char,admin1.name,admin2.char,admin2.name) %>% 
    group_by(cluster,U,admin1,admin1.char,admin1.name,admin2.char,admin2.name) %>% 
    reframe(Z=sum(Y),n=sum(total)) %>% rename(A1=admin1)
  alg_dat <- merge(alg_dat,admin_key)
}

# estimate average number of HH in SAMPLED cluster
{
    frame.info <- merge(dhs_dat,alg_dat %>% group_by(admin2.char,U) %>% summarise(n_clusters_samp = n())) %>% 
      dplyr::select(A2,admin2.char,U,M_bar,n_clusters,n_clusters_samp,n_HH_samp)
    sd_M <- 0.1*frame.info$M_bar
    
    n_iter <- 1000
    M_obs_bar <- rep(NA,n_iter)
    frame.info$M_bar_est <- NA
    
    for(i in 1:nrow(frame.info)){
      for(k in 1:n_iter){
        M = round(rnorm(frame.info$n_clusters[i],frame.info$M_bar[i],sd_M[i]))
        M_obs <- M[sample(1:frame.info$n_clusters[i],frame.info$n_clusters_samp[i],prob = M)]
        M_obs_bar[k] <- mean(M_obs)
      }
      frame.info$M_bar_est[i] <- mean(M_obs_bar)
    }
}

# add estimated number of total births per cluster
alg_dat <- merge(alg_dat,frame.info) %>% dplyr::select(-M_bar,-n_clusters,-n_clusters_samp) %>%
  mutate(N=round(M_bar_est*n/n_HH_samp))

# organize population and survey info
pop_strata_dat <- merge(dhs_dat %>% dplyr::select(A1,A2,admin1.name,admin2.name,U,M),
                            alg_dat %>% group_by(A2,U) %>% summarise(birth_rate = sum(n)/sum(n_HH_samp)))
pop_strata_dat$N <- round(pop_strata_dat$M*pop_strata_dat$birth_rate)
pop_strata_dat <- pop_strata_dat[order(pop_strata_dat$A2,1-pop_strata_dat$U),]
  
data_list <- list(obs_dat = alg_dat,
                    pop_strata_dat = as.data.frame(pop_strata_dat),
                    Q_scaled_inv = Q_scaled_inv)

## INLA model: Admin2 BYM2 -----
admin2.fit <- INLA::inla(Z ~ factor(U) + f(A2, model="bym2", graph = admin2.mat, scale.model=T, constr=T,
                                       hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1))) -1,
                         family = 'binomial',
                         data = data_list$obs_dat,
                         control.family=list(variant=1,link='log'),
                         control.predictor = list(compute = FALSE, link = 1),
                         control.compute = list(config=T),
                         Ntrials = data_list$obs_dat$n)

# sample from posterior
cs <- admin2.fit$misc$configs$contents$tag
cs <- cs[cs != "Predictor"]
select <- list()
for (i in 1:length(cs)) {
  select[[i]] <- 0
  names(select)[i] <- cs[i]
}
sampFull <- INLA::inla.posterior.sample(n = 1000, result = admin2.fit, intern = TRUE, selection = select)
sampFull.draws <- matrix(NA,1000,length(sampFull[[1]]$latent))
for(i in 1:1000){
  sampFull.draws[i,] <- sampFull[[i]]$latent
}
fields <- colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)

inla_r_samples <- exp(sampFull.draws[,1:nrow(admin_key)] +
  sampFull.draws[,2*nrow(admin_key)+1]%*%t(adm2_UR_weights$urb_frac) +
  sampFull.draws[,2*nrow(admin_key)+2]%*%t(1-adm2_UR_weights$urb_frac))

admin2.inla <- data.frame(admin2 = 1:nrow(admin_key))
admin2.inla$admin2.char <- paste0('admin2_', admin2.inla$admin2)

admin2.inla$adm2_median <- (apply(inla_r_samples,2,function(x){median(exp(x))}))
admin2.inla$adm2_lower95 <- (apply(inla_r_samples,2,function(x){quantile(exp(x),0.025)}))
admin2.inla$adm2_upper95 <- (apply(inla_r_samples,2,function(x){quantile(exp(x),0.975)}))


## run NUTS Alg 5.2a ------
data_list$obs_dat$A <- data_list$obs_dat$A2
data_list$pop_strata_dat$A <- data_list$pop_strata_dat$A2

start.time <- Sys.time()
NUTS_5.2a <- postsamp_Alg5.2a_NUTS(eps0=0.01,
                                   M_diag = 1/c(1,2,1,rep(0.5,30),500),
                                   logit_r_hat=natl.dir$logit.est,logit_r_hat_var=natl.dir$var.est,
                                   data_list=data_list, n_iter=10000,warmup=50,chains=4)
Sys.time() - start.time

#save(NUTS_5.2a,file="Handcoded MCMC Simulations/Alg 5.2a BYM2 NUTS 240223, MWI 2011-2015, 4 chains of 10k iterations posterior dist.rda")
#load(file="/Users/alanamcgovern/Desktop/Research/New Benchmarking/Handcoded MCMC Simulations/Alg 5.2a BYM2 NUTS 240223, MWI 2011-2015, 4 chains of 10k iterations posterior dist.rda")

## sort Alg results
r_samples <- adm2_UR_weights$urb_frac*NUTS_5.2a$r[,1:28] + (1-adm2_UR_weights$urb_frac)*NUTS_5.2a$r[,29:56]
#take our burnin
r_samples <- r_samples[c(2001:10001,12002:20002,22003:30003,32004:40004),]

res_r <- data.frame(A2=1:28)
res_r$admin2.char <- paste0('admin2_',res_r$A2)

res_r$median <- apply(r_samples,2,median)
res_r$lower <- apply(r_samples,2,quantile,probs=0.025)
res_r$upper <-  apply(r_samples,2,quantile,probs=0.975)

## Compare Alg5.2a with corresponding INLA model ------
plot(res_r$median,admin2.inla$adm2_median,xlab='Alg 2.3a',ylab = 'INLA')
par(mfrow=c(2,2))
for(i in 1:28){
  hist(inla_r_samples[,i],probability = T)
  hist(r_samples[,i],add=T,probability = T,col='skyblue')
  
}

## INLA model: Admin2 BYM2 with Admin1 fixed effect ------

admin2FE.fit <- INLA::inla(Z ~ factor(U) + factor(A1) + f(A2, model="bym2", graph = admin2.mat, scale.model=T, constr=T,
                                                          adjust.for.con.comp=T, 
                                                          hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1))) -1,
                           family = 'binomial', 
                           data = data_list$obs_dat, 
                           control.family=list(variant=1,link='log'),
                           control.predictor = list(compute = FALSE, link = 1), 
                           control.compute = list(config=T),
                           Ntrials = data_list$obs_dat$n)

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

inla2_r_samples <- matrix(0,ncol=(max(admin_key$A1)+2*max(admin_key$A2)+1),nrow=max(admin_key$A2))
for(area in 1:max(admin_key$A2)){
  inla2_r_samples[area,area] <- 1 #choose admin1 area
  if(admin_key[admin_key$A2==area,]$A1>1){
    inla2_r_samples[area,2*max(admin_key$A2)+1+admin_key[admin_key$A2==area,]$A1] <- 1 #choose admin1 area
  }
}
inla2_r_samples[,2*max(admin_key$A2)+1] <- adm2_UR_weights$urb_frac
inla2_r_samples[,2*max(admin_key$A2)+2] <- 1 - adm2_UR_weights$urb_frac

inla2_r_samples <- inla2_r_samples%*%t(sampFull.draws)
admin2.inla2 <- data.frame(admin2 = 1:nrow(admin_key))
admin2.inla2$admin2.char <- paste0('admin2_', admin2.inla2$admin2)

admin2.inla2$adm2_median <- (apply(inla2_r_samples,1,function(x){median(exp(x))}))
admin2.inla2$adm2_lower95 <- (apply(inla2_r_samples,1,function(x){quantile(exp(x),0.025)}))
admin2.inla2$adm2_upper95 <- (apply(inla2_r_samples,1,function(x){quantile(exp(x),0.975)}))

admin1.inla2 <- data.frame(admin1 = rep(sort(unique(admin_key$A1)),each=2), U = rep(c(0,1),max(admin_key$A1)))
admin1.inla2$est <- rep(c(0,admin2FE.fit$summary.fixed$`0.5quant`[3:nrow(admin2FE.fit$summary.fixed)]),each=2) + 
  rep(admin2FE.fit$summary.fixed$`0.5quant`[c(2,1)],max(admin_key$A1))
admin1.inla2$admin1.char <- paste0('admin1_',admin1.inla2$admin1)
admin1.inla2 <- merge(admin1.inla2,adm1_UR_weights,by.x='admin1.char',by.y='adm_idx')
admin1.inla2[admin1.inla2$U==0,]$urb_frac <- 1-admin1.inla2[admin1.inla2$U==0,]$urb_frac
admin1.inla2 <- admin1.inla2 %>% group_by(admin1.char,admin1) %>% 
  summarise(adm1_median=exp(sum(est*urb_frac)))

all_res <- merge(admin_key,admin2.inla2,by=c('admin2.char'))
all_res <- merge(all_res,admin1.inla2,by=c('admin1.char'))
all_res <- merge(all_res,dir.est,by.x='admin1.char',by.y='region')

g <- all_res %>% ggplot() + geom_point(aes(mean,admin1.char,col='1',pch='1'),size=3) +
  geom_point(aes(lower,admin1.char,col='1'),pch=3,size=3) +
  geom_point(aes(upper,admin1.char,col='1'),pch=3,size=3) +
  geom_point(aes(adm2_median,admin1.char,col='2',pch='2'),size=3) + 
  geom_point(aes(adm1_median,admin1.char,col='3',pch='3'),size=3) +
  #geom_point(aes(admin1_inla_median,admin1.char,col='3',pch='3'),size=3) +
  ggtitle('Malawi 2011-2015') + theme_minimal() +
  ylab('') + xlab('NMR') + theme(legend.position = 'bottom') +
  scale_colour_manual(name = '', values =c('1'='blue','2'='orange','3'='darkgreen'),
                      labels = c('Admin 1 Direct','Admin 2 BYM2 (with Admin 1 FE)','Admin 1 FE')) +
  scale_shape_manual(name = '', values =c('1'=17,'2'=1,'3'=2),
                     labels = c('Admin 1 Direct','Admin 2 BYM2 (with Admin 1 FE)','Admin 1 FE'))

## run NUTS Alg 5.3a -------

#run model to get estimate of covariance matrix
inla.fit <- INLA::inla(Z ~ factor(U) + f(A2, model="bym2", graph = admin2.mat, scale.model=T, constr=T,
                                     adjust.for.con.comp=T, 
                                     hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1))) -1,
                       family = 'nbinomial', data = data_list$obs_dat, 
                       control.family=list(link='log'),
                       control.predictor = list(compute = FALSE, link = 1), 
                       control.compute = list(config=T),
                       E = data_list$obs_dat$n)

# sample from posterior
cs <- inla.fit$misc$configs$contents$tag
cs <- cs[cs != "Predictor"]
select <- list()
for (i in 1:length(cs)) {
  select[[i]] <- 0
  names(select)[i] <- cs[i]
}
sampFull <- INLA::inla.posterior.sample(n = 1000, result = inla.fit, intern = TRUE, selection = select)
sampFull.draws <- matrix(NA,1000,5+max(pop_strata_dat$A2))
for(i in 1:1000){
  sampFull.draws[i,] <- c(sampFull[[i]]$hyperpar[2:3],
                          -sampFull[[i]]$hyperpar[1], #kind of comparable to log(d)
                          tail(sampFull[[i]]$latent,2),sampFull[[i]]$latent[c(1:max(pop_strata_dat$A2))])
}
Sigma <- cov(sampFull.draws)

start.time <- Sys.time()
NUTS_5.3a <- postsamp_Alg5.3a_NUTS_MCMC(eps0 = 0.05, 
                                        M_init = solve(Sigma),
                                        logit_r_hat = dir.est$logit.est,
                                        logit_r_hat_var = dir.est$var.est,
                                        data_list = data_list, 
                                        n_iter = 200, chains=4)
Sys.time() - start.time
save(NUTS_5.3a,file="Handcoded MCMC Simulations/Alg 5.3a BYM2 NUTS 240228, MWI 2011-2015, 4 chains of 2k iterations posterior dist.rda")

iter.include <- c(500:1000,1500:2000,2500:3000,3500:4000)
iter.include <- c(500:2000,2500:4000,4500:6000,6500:8000)
iter.include <- c(50:200,250:400,450:600,650:800)

plot(NUTS_5.3a$logtau[iter.include],type='l')
plot(NUTS_5.3a$logitphi[iter.include],type='l')
plot(NUTS_5.3a$logd[iter.include],type='l')
plot(NUTS_5.3a$alphaU[iter.include],type='l')
plot(NUTS_5.3a$alphaR[iter.include],type='l')
plot(NUTS_5.3a$b[iter.include,1],type='l')
plot(NUTS_5.3a$b[iter.include,17],type='l')
plot(NUTS_5.3a$Yplus[iter.include,1],type='l')
plot(NUTS_5.3a$Yplus[iter.include,2],type='l')
plot(NUTS_5.3a$Yplus[iter.include,3],type='l')


each_chain <- nrow(NUTS_5.3a$reg_params)/4
par(mfrow=c(1,3))
for(k in 1:ncol(NUTS_5.3a$Yplus)){
  plot(cumsum(NUTS_5.3a$Yplus[1:each_chain,k])/(1:each_chain),type='l',col=1,ylab='',
       ylim=c(min(NUTS_5.3a$Yplus[,k]),max(NUTS_5.3a$Yplus[,k])),main = paste0('Running mean of Yplus_',k))
  for(c in 2:chains){
    lines(cumsum(NUTS_5.3a$Yplus[((c-1)*each_chain+1):(c*each_chain),k])/(1:each_chain),col=c)
  }
}


adm2_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_MWI/admin2_u1_2013_urban_frac.rds"))
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Malawi/worldpop/adm2_weights_u1.rda"))
adm2_weights <- weight.adm2.u1[weight.adm2.u1$years==2013,]

strat.res.dat <- pop_dat[order(pop_dat$U,decreasing=T),]
strat.res.dat$mean <- colMeans(NUTS_5.3a$r[iter.include,],na.rm = T)
strat.res.dat$lower <- apply(NUTS_5.3a$r[iter.include,],2,quantile,probs=0.025)
strat.res.dat$upper <- apply(NUTS_5.3a$r[iter.include,],2,quantile,probs=0.975)
strat.res.dat$median <- colMedians(NUTS_5.3a$r[iter.include,],na.rm = T)

# fix later to get bounds
adm2.res.dat <- merge(strat.res.dat,adm2_UR_weights[,2:3],by.x='admin2.char',by.y='adm_idx') %>% mutate(urb_frac = ifelse(U==0,1-urb_frac,urb_frac)) %>% 
  group_by(admin2.char,admin1.char) %>% summarise(mean = sum(urb_frac*mean),median = sum(urb_frac*median))


ggplot() +# geom_point(data=adm2.res.dat,aes(mean,admin1.char,col='2',pch='2'),size=4) +
  geom_point(data=dir.est,aes(mean,region,col='1',pch='1'),size=4) + # admin1 direct
  geom_point(data=dir.est,aes(lower,region,col='1'),pch=3,size=3) + # admin1 direct lower
  geom_point(data=dir.est,aes(upper,region,col='1'),pch=3,size=3) + # admin1 direct upper
  geom_point(data=sd.est,aes(median,region)) +
  geom_point(data=sd.est,aes(lower,region),pch=3) +
  geom_point(data=sd.est,aes(upper,region),pch=3) +
  geom_point(data=all_res,aes(adm1_median,admin1.char,col='4',pch='4'),size=4) +
  geom_point(data=all_res,aes(adm2_median,admin1.char)) +
  ggtitle('Malawi 2011-2015') + theme_minimal() +
  ylab('') + xlab('NMR') + theme(legend.position = 'bottom') 
  scale_colour_manual(name = '', values =c('1'='blue','3'='red','2'='orange','4'='red'),
                      labels = c('Admin 1 Direct','Admin 2','Admin 2, Aggregated','Admin 2 INLA, Aggregated')) +
  scale_shape_manual(name = '', values =c('1'=17,'2'=1,'3'=2,'4'=17),
                     labels = c('Admin 1 Direct','Admin 2','Admin 2, Aggregated','Admin 2 INLA, Aggregated'))
