library(INLA)
library(SUMMER)
library(tidyverse)
library(ggpubr)
library(Rfast)
source('getAggregated.R')
inla.setOption(inla.timeout = 10)


out <- failure.summary <- NULL
countries <- c('Burundi','Ghana','Guinea','Lesotho','Liberia','Malawi','Pakistan','Rwanda','Sierra_Leone')
shapeFileFolders <- c('/shapeFiles','/shapeFiles/gadm40_GHA_shp','/shapeFiles/gadm41_GIN_shp',
                      '/shapeFiles/gadm41_LSO_shp','/shapeFiles/gadm41_LBR_shp','shapeFiles_gadm',
                      '/shapeFiles/gadm41_PAK_shp','/shapeFiles/gadm41_RWA_shp','shapeFiles')

# simulate data and get compare performance of direct and INLA estimate
for(country_t in countries){
  ## load data ------
  message(paste0('Starting ',country_t,'\n'))
  
  data.dir <- (paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,'/'))
  if(file.exists(paste0(data.dir,country_t,"_cluster_dat_1frame.rda"))){
    load(paste0(data.dir,country_t,"_cluster_dat_1frame.rda"))
  }else{
    load(paste0(data.dir,country_t,"_cluster_dat.rda"))
  }
  
  dat <- mod.dat[mod.dat$age==0,]
  dat$died <- dat$Y
  dat$years.int <- as.integer(dat$years)
  dat$years <- as.numeric(as.character(dat$years))
  dat$v005 <- dat$v005/1e6
  
  shapeFilePath <- shapeFileFolders[which(country_t==countries)]
  load(paste0(data.dir,shapeFilePath,'/',country_t,'_Amat.rda'))
  
  survey_year <- max(dat$survey)
  if(country_t=='Ghana')
    survey_year <- 2014
  dat <- dat[dat$survey==survey_year,]
  
  # assign time periods 
  end.year <- survey_year
  if(((end.year-2000+1) %% 3)==0){
    beg.period.years <- seq(2000,end.year,3) 
    end.period.years <- beg.period.years + 2 
  }else if(((end.year-2000+1) %% 3)==1){
    beg.period.years <- c(2000,2000+2,seq(2000+4,end.year,3))
    end.period.years <- c(2000+1,2000+3,seq(2000+6,end.year,3))
  }else if(((end.year-2000+1) %% 3)==2){
    beg.period.years <- c(2000,seq(2000+2,end.year,3))
    end.period.years <- c(2000+1,seq(2000+4,end.year,3))
  }
  
  periods <- paste(beg.period.years, end.period.years, sep = "-") # convert the periods into string
  dat$period <- as.character(cut(dat$years, breaks = c(beg.period.years, beg.period.years[length(beg.period.years)]+5),
                                 include.lowest = T, right = F, labels = periods))
  
  ## simulate population -----
  
  n_areas <- max(dat$admin1)
  dat$g_hi <- 30 # how many HH were sampled
  simdat <- dat %>% filter(period==periods[4]) %>% group_by(period,cluster,strata,urban,age,admin1.char,admin1,admin1.name) %>% 
    summarise(total = sum(total),g_hi=sum(g_hi))
  # num sampled clusters in strata
  strata_dat <- simdat %>% group_by(admin1.char,admin1.name,urban) %>% summarise(a_h = length(unique(cluster))) %>%
    # simulate number of clusters in strata
    mutate(A_h = ifelse(urban=='urban',30*a_h,40*a_h))
  # simulate number of HH in strata
  strata_dat$M_h = sapply(1:nrow(strata_dat),function(i){round(rnorm(1,150*strata_dat$A_h[i],25))})
  
  simdat <- merge(simdat,strata_dat)
  
  #define number of total births in each cluster
  simdat$N <- round(simdat$total/simdat$g_hi*rnorm(1,150,15))
  simdat$v005 <- simdat$M_h/(simdat$a_h*simdat$g_hi)/100
  
  ## run simulations ------
  
  # set parameters
  alphaU <- -3.5
  alphaR <- -3.3
  beta.adm1 <- rnorm(n_areas,0,0.1)
  alpha <- 5 #overdispersion parameter
  
  param_dat <- simdat %>% group_by(admin1,region = admin1.char) %>% summarise(rural = sum(N*I(urban=='rural'))/sum(N))
  param_dat$urban <- 1-param_dat$rural
  param_dat$beta <- beta.adm1
  param_dat$r <- exp(alphaU + param_dat$beta)*param_dat$urban + exp(alphaR + param_dat$beta)*param_dat$rural
  
  nsim <- 50
  coverage.dir.mat <- width.dir.mat <- deathcount.dir.mat <- bias.dir.mat <- coefvar.dir.mat <- matrix(0,nsim,n_areas)
  coverage.sd.mat <- width.sd.mat <- deathcount.sd.mat <- bias.sd.mat <- coefvar.sd.mat <- matrix(0,nsim,n_areas)
  coverage.inla1.mat <- width.inla1.mat <- deathcount.inla1.mat <- bias.inla1.mat <- coefvar.inla1.mat <- matrix(0,nsim,n_areas)
  failure <- nodeaths <- 0
  for(k in 1:nsim){
    if(k%%10==0)
     cat('Starting simulation',k,'\n')
    
    # simulate total deaths in cluster 
    simdat$total_died <- sapply(1:nrow(simdat),function(i){
      rnbinom(1,alpha,alpha/(alpha + simdat$N[i]*exp(alphaU*I(simdat$urban[i]=='urban') + alphaR*I(simdat$urban[i]=='rural') + beta.adm1[simdat$admin1[i]] + rnorm(1,0,0.01))))})
    # simulate observed deaths in cluster
    simdat$died <- sapply(1:nrow(simdat),function(i){
      rhyper(1,simdat$total_died[i],simdat$N[i]-simdat$total_died[i],simdat$total[i])})
    
    ## don't fit models if any area has no observed deaths
    if(sum((simdat %>% group_by(admin1.char) %>% summarise(deaths = sum(died)))$deaths==0)>0){
      inla1.est <- sd.est <- dir.est <- data.frame(admin1 = 1:n_areas,mean=NA,median=NA,lower=NA,upper=NA)
      failure <- failure + 1
      nodeaths <- nodeaths + 1
    }else{
      # INLA estimates
      inla.mod1 <- try({
        INLA::inla(died ~ factor(urban) + factor(admin1) -1,
                   data=simdat, family='nbinomial', E=total,
                   control.predictor = list(compute = T, link = 1),
                   control.family = list(control.link = list(model = 'log'),
                                         hyper = list(theta = list(prior="pc.mgamma", param=0.25))),
                   control.compute = list(config=T))})
      ## stop model if its stalling
      if(is.character(inla.mod1[1])){
        inla1.est <- sd.est <- dir.est <- data.frame(admin1 = 1:n_areas,mean=NA,median=NA,lower=NA,upper=NA)
        failure <- failure + 1
        cat('Iteration ',k,'. Moving on.. \n')
      }else{
        inla.out <- getAggregated(inla.mod1,simdat,admin1.UR.weights = param_dat)
        inla1.est <- inla.out$adm1.est
        
        # direct estimate
        suppressMessages({
          dir.out <- SUMMER::getDirect(simdat, periods[4],regionVar = "admin1.char",
                                       timeVar = "period", clusterVar =  "~cluster",
                                       ageVar = "age", Ntrials = "total",
                                       weightsVar = "v005",national.only = F)
        })
        #put in correct order
        dir.out <- dir.out[dir.out$region!='All',]
        dir.est <- dir.out
        dir.est$region_num <- unlist(lapply(str_split(dir.est$region,'_'),function(x){as.numeric(x[2])}))
        dir.est <- dir.est[order(dir.est$region_num),]
        
        # smoothed direct estimate
        sd.fit <- smoothDirect(dir.est, Amat = admin1.mat,
                               year_label = 'years', time.model = NULL,
                               year_range = periods[4], is.yearly = F)
        sd.est <- getSmoothed(sd.fit, Amat = admin1.mat, year_label = 'years', year_range = periods[4])
        sd.est$region_num <- unlist(lapply(str_split(sd.est$region,'_'),function(x){as.numeric(x[2])}))
        sd.est <- sd.est[order(sd.est$region_num),]
        
      }
    }
    
    # record direct estimate info
    coverage.dir.mat[k,] <- I(dir.est$lower<param_dat$r & dir.est$upper>param_dat$r)
    deathcount.dir.mat[k,] <- (simdat %>% group_by(admin1) %>% summarise(deaths = sum(died)))$deaths
    width.dir.mat[k,] <- dir.est$upper - dir.est$lower
    coefvar.dir.mat[k,] <- (dir.est$upper - dir.est$lower)/(2*1.96*dir.est$mean)
    bias.dir.mat[k,] <- dir.est$mean - param_dat$r
    
    # record smoothed direct estimate info
    coverage.sd.mat[k,] <- I(sd.est$lower<param_dat$r & sd.est$upper>param_dat$r)
    deathcount.sd.mat[k,] <- (simdat %>% group_by(admin1) %>% summarise(deaths = sum(died)))$deaths
    width.sd.mat[k,] <- sd.est$upper - sd.est$lower
    coefvar.sd.mat[k,] <- (sd.est$upper - sd.est$lower)/(2*1.96*sd.est$median)
    bias.sd.mat[k,] <- sd.est$median - param_dat$r
    
    # record INLA estimate info
    coverage.inla1.mat[k,] <- I(inla1.est$lower<param_dat$r & inla1.est$upper>param_dat$r)
    deathcount.inla1.mat[k,] <- (simdat %>% group_by(admin1) %>% summarise(deaths = sum(died)))$deaths
    width.inla1.mat[k,] <- inla1.est$upper - inla1.est$lower
    coefvar.inla1.mat[k,] <- (inla1.est$upper - inla1.est$lower)/(2*1.96*inla1.est$median)
    bias.inla1.mat[k,] <- inla1.est$median - param_dat$r
  }
  
  area_summary <- simdat %>% group_by(region = admin1) %>% summarise(clusternum = length(unique(cluster)), obs_births = sum(total))
  area_summary$r <- param_dat$r
  area_summary$country <- country_t
  area_summary$country.num <- as.factor(which(country_t==countries))
  
  direct.summary <- sd.summary <- inla1.summary <- area_summary
  
  direct.summary$coverage = colMeans(coverage.dir.mat,na.rm=T)
  direct.summary$width = colMeans(width.dir.mat,na.rm=T)
  direct.summary$bias = colMeans(bias.dir.mat,na.rm = T)
  direct.summary$coefvar = colMeans(coefvar.dir.mat,na.rm = T)
  direct.summary$model = 'Direct'
  
  sd.summary$coverage = colMeans(coverage.sd.mat,na.rm=T)
  sd.summary$width = colMeans(width.sd.mat,na.rm=T)
  sd.summary$bias = colMeans(bias.sd.mat,na.rm = T)
  sd.summary$coefvar = colMeans(coefvar.sd.mat,na.rm = T)
  sd.summary$model = 'SD'
  
  inla1.summary$coverage = colMeans(coverage.inla1.mat,na.rm=T)
  inla1.summary$width = colMeans(width.inla1.mat,na.rm=T)
  inla1.summary$bias = colMeans(bias.inla1.mat^2,na.rm = T)
  inla1.summary$coefvar = colMeans(coefvar.inla1.mat,na.rm = T)
  inla1.summary$model = 'INLA'
 
  failure.dat <- data.frame(country=country_t,failures=failure/nsim)
  failure.summary <- rbind(failure.summary,failure.dat)
  
  out <- rbind(out,direct.summary,sd.summary,inla1.summary)
}

save(out,file = 'Evaluating Direct Estimates/Compare_Direct_Simulations_231121.rda')
#load('Evaluating Direct Estimates/Compare_Direct_Simulations_231121.rda')

pdf('Evaluating Direct Estimates/Compare_Direct_Simulations_Figures_231122.pdf')
{
# look at relationship between num of observed births and model performance (strong association) ------

# coverage
b1 <- out %>% 
    #filter(model!='INLA') %>% 
    ggplot() + geom_point(aes(x=obs_births,y=coverage,pch=model,col=country.num),size=1) + geom_hline(yintercept = 0.95,col='red') +
  xlab('Number of observed births') + ylab('Coverage') + scale_x_log10() +scale_shape_manual(name='Model', values = c(3,1,2)) +
    ggtitle('Coverage of 95% CI of Admin 1 level Estimates')

# bias
b2 <- out %>% 
  #filter(model!='Direct') %>% 
  ggplot() + geom_point(aes(x=obs_births,y=bias,pch=model,col=country.num),size=1) + 
  xlab('Number of observed births') +ylab('Average bias') + scale_x_log10() +scale_shape_manual(name='Model', values = c(3,1,2)) + ggtitle('Average bias of Admin 1 level Estimates')

# coefficient of variation
b3 <- out %>% 
  #filter(model!='Direct') %>%
  ggplot() + geom_point(aes(x=obs_births,y=coefvar,pch=model,col=country.num),size=1) + scale_x_log10() +
  xlab('Number of observed births') + ylab('Average coefficient of variation')+scale_shape_manual(name='Model', values = c(3,1,2))+ ggtitle('Average coefficient of variation of Admin 1 level Estimates')

print(b1)
print(b2)
print(b3)


# look at relationship between true NMR and model performance (strong association) ------
# coverage
r1 <- out %>% ggplot() + geom_point(aes(x=r,y=coverage,pch=model,col=country.num),size=1) + geom_hline(yintercept = 0.95,col='red') +
  xlab('true NMR') + ylab('Coverage') + scale_x_log10() +scale_shape_manual(name='Model', values = c(3,1,2))

# bias
r2 <- out %>% ggplot() + geom_point(aes(x=r,y=bias,pch=model,col=country.num),size=1) + 
  xlab('true NMR') +ylab('Average bias') + scale_x_log10() +scale_shape_manual(name='Model', values = c(3,1,2))

#coefficient of variation
r3 <- out %>% ggplot() + geom_point(aes(x=r,y=coefvar,pch=model,group=country.num,col=country.num),size=1) + scale_x_log10() +
  xlab('true NMR') + ylab('Coefficient of Variation')+scale_shape_manual(name='Model', values = c(3,1,2))

rfig <- ggarrange(plotlist = list(r1,r2,r3),common.legend = T,legend='right')

print(annotate_figure(rfig, top = text_grob("Relationship between Admin 1 level estimates and true NMR",size = 12)))
}
dev.off()
