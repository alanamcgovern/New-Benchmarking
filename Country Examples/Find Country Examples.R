library(tidyverse)
library(locfit)
library(rdhs)
library(surveyPrev)
library(SUMMER)
library(labelled)
library(sf)
library(ggpubr)

set_rdhs_config(email = "amcgov@uw.edu",
                project = "Spatial Modeling for Subnational Administrative Level 2 Small-Area Estimation - Under 5 Mortality Rate")
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]
setwd(paste(code.path.splitted[1: (length(code.path.splitted)-1)], collapse = "/"))


easy_countries <- c('Angola','Cameroon','Ethiopia','Guinea','Haiti','Liberia','Madagascar','Myanmar',
               'Mauritania','Rwanda','Chad','Tanzania','Zambia','Zimbabwe')
my_countries <- c('Burundi','Malawi','Sierra_Leone')
noAdm2_countries <- c('Kenya', 'Lesotho', 'Nigeria', 'Ghana', 'Togo') # no admin2 UR fractions bc analysis was not done at admin2 level

# has admin 1 areas without deaths: Myanmar (3 year range), Kenya (3 year range), Tanzania (2010 and 2022)
# problem with polygon files -- Nepal, 'Benin', 'Bangladesh', 
# geo function/other surveyPrev stuff doesn't work: 'Pakistan', 'Senegal', 'Uganda' (no admin2)
# no UR fractions at all - 'Namibia'
# not all admin 1 areas are showing up in DHS data? Mali

countries <- c(easy_countries,my_countries)
#countries <- countries[18:length(countries)]

for(country_t in countries){
  plot_list <- list()
  # load country info
  load(file = paste0('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Info/',paste0(country_t, "_general_info.Rdata"))) # load the country info
  countryId <- dhs_countries()[dhs_countries()$ISO3_CountryCode==toupper(gadm.abbrev),]
  potential_surveys <- dhs_datasets(countryIds = countryId$DHS_CountryCode, surveyYearStart = 2010, surveyType = 'DHS') %>% 
    dplyr::filter((FileType == 'Births Recode' & FileFormat=='Stata dataset (.dta)'))
  survey_years <- as.numeric(potential_surveys$SurveyYear)
  for(survey_year in survey_years){
    for(year_range in c(3,5)){
      beg.year <- survey_year - year_range + 1
      mid.year <- survey_year - (year_range + 1)/2
      
      # load shapefiles -----
      # if country uses regular shapefiles
      poly.path <- paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_",gadm.abbrev,"_shp")
      # other cases
      if(!dir.exists(poly.path)){
        poly.path <- paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_alt/",country_t)
      }
      
      if(country_t %in% easy_countries){
        poly.adm1 <- st_read(dsn = poly.path, layer = poly.layer.adm1, options = "ENCODING=UTF-8")
        poly.adm2 <- st_read(dsn = poly.path, layer = poly.layer.adm2, options = "ENCODING=UTF-8")
        
        geo <- getDHSgeo(country = country_t, year = survey_year)
        
        cluster.info <- clusterInfo(geo = geo, poly.adm1 = poly.adm1, poly.adm2 = poly.adm2)
        admin1.info <- adminInfo(geo = poly.adm1, admin = 1)
        admin1.info$admin.info$admin1 <- 1:nrow(admin1.info$admin.info)
        admin1.info$admin.info$admin1.char <- paste0('admin1_',admin1.info$admin.info$admin1)
        row.names(admin1.info$admin.mat) <- colnames(admin1.info$admin.mat) <- admin1.info$admin.info$admin1
        
        admin2.info <- adminInfo(geo = poly.adm2, admin = 2)
        admin2.info$admin.info$admin2 <- 1:nrow(admin2.info$admin.info)
        admin2.info$admin.info$admin2.char <- paste0('admin2_',admin2.info$admin.info$admin2)
        row.names(admin2.info$admin.mat) <- colnames(admin2.info$admin.mat) <- admin2.info$admin.info$admin2
        admin.key <- merge(admin1.info$admin.info[,c(1,4,5)],admin2.info$admin.info[,c(1,2,6,7)])
        admin.key <- admin.key[order(admin.key$admin2),]
        cluster.info$cluster.info <- merge(cluster.info$cluster.info[,c(1,5,6)],admin.key)
        
        #make admin2 mat nested
        for(i in 1:nrow(admin2.info$admin.mat)){
          admin2.info$admin.mat[i,which(admin.key$admin1!=admin.key$admin1[i])] <- 0
          if(sum(admin2.info$admin.mat[i,])>0){
            admin2.info$admin.mat[i,] <- admin2.info$admin.mat[i,]/sum(admin2.info$admin.mat[i,])
          }
        }
        
      }else if(country_t %in% my_countries){
        load(paste0(poly.path,'/',country_t,'_Amat.rda'))
        admin2.info <- list()
        admin2.info$admin.mat <- admin2.mat
        admin.key <- readxl::read_excel("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Admin1_Admin2_Key.xlsx") %>% 
          filter(country==country_t) %>% select(-country)
        admin.key$admin1 <- as.numeric(str_remove(admin.key$admin1.char,'admin1_'))
        admin.key$admin2 <- as.numeric(str_remove(admin.key$admin2.char,'admin2_'))
        admin.key <- admin.key[order(admin.key$admin2),]
        
        row.names(admin2.info$admin.mat) <- colnames(admin2.info$admin.mat) <- 1:nrow(admin2.mat)
        #make admin2 mat nested
        for(i in 1:nrow(admin2.info$admin.mat)){
          admin2.info$admin.mat[i,which(admin.key$admin1!=admin.key$admin1[i])] <- 0
          if(sum(admin2.info$admin.mat[i,])>0){
            admin2.info$admin.mat[i,] <- admin2.info$admin.mat[i,]/sum(admin2.info$admin.mat[i,])
          }
        }
      }else if(country %in% noAdm2_countries){
        poly.adm1 <- st_read(dsn = poly.path, layer = poly.layer.adm1, options = "ENCODING=UTF-8")
        
        geo <- getDHSgeo(country = country_t, year = survey_year)
        
        cluster.info <- clusterInfo(geo = geo, poly.adm1 = poly.adm1, poly.adm2 = poly.adm1)
        admin1.info <- adminInfo(geo = poly.adm1, admin = 1)
        admin1.info$admin.info$admin1 <- 1:nrow(admin1.info$admin.info)
        admin1.info$admin.info$admin1.char <- paste0('admin1_',admin1.info$admin.info$admin1)
        row.names(admin1.info$admin.mat) <- colnames(admin1.info$admin.mat) <- admin1.info$admin.info$admin1
        
        admin.key <- admin1.info$admin.info[,c(1,4,5)]
        cluster.info$cluster.info <- merge(cluster.info$cluster.info[,c(1,5,6)],admin.key)
      }
      
      # prepare NMR data ------
      
      if(country_t %in% c(easy_countries,noAdm2_countries)){
        data.paths.tmp <- get_datasets(potential_surveys[potential_surveys$SurveyYear==survey_year,]$FileName, clear_cache = T)
        raw.dat <- readRDS(paste0(data.paths.tmp))
        
        # convert some variables to factors
        alive <- attr(raw.dat$b5, which = "labels")
        names(alive) <- tolower(names(alive))
        raw.dat$b5 <- ifelse(raw.dat$b5 == alive["yes"][[1]], "yes", "no")
        raw.dat$b5 <- factor(raw.dat$b5, levels = c("yes", "no"))
        
        strat <- attr(raw.dat$v025,which='labels')
        names(strat) <- tolower(names(strat))
        raw.dat$v025 <- ifelse(raw.dat$v025 == strat["urban"][[1]],'urban','rural')
        raw.dat$v025 <- factor(raw.dat$v025, levels = c('urban','rural'))
        
        if(country=='Ethiopia'){
          cmc.adjust <- 92
        }else if(country=='Nepal'){
          cmc.adjust <- -681
        }else{cmc.adjust <- 0}
        
        if(country=='Pakistan' & survey_year==2017){
          raw.dat[raw.dat$v005==0,]$v005 <- raw.dat[raw.dat$v005==0,]$sv005
        }
        
        # read DHS data
        u5.dat <- getBirths(data=raw.dat,
                            surveyyear = survey_year,
                            year.cut = seq(beg.year, survey_year + 1, 1),
                            cmc.adjust = cmc.adjust,compact = T)
        
        # retrieve the some columns of the full data
        u5.dat <- u5.dat[ ,c("v001", "v024", "time", "total", "age", "v005", "v025", "strata", "died")]
        
        # remove wrong points in the data if any
        u5.dat <- u5.dat[!(u5.dat$v001 %in% cluster.info$wrong.points),]
        
        # label variables and add admin areas
        nmr.dat <- u5.dat %>% filter(age==0) %>% dplyr::select(v001,time,total,v005,v025,strata,died)
        colnames(nmr.dat) <- c("cluster","years", "total","v005", "urban","strata","Y")
        nmr.dat <- merge(nmr.dat,cluster.info$cluster.info)
      }else if(country_t %in% my_countries){
        load(paste0('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/',country_t,'/',country_t, "_cluster_dat.rda"))
        nmr.dat <- mod.dat %>% filter(survey==survey_year,years %in% beg.year:survey_year,age=='0')
      }
      
      if(country_t %in% noAdm2_countries){
        nmr.dat <- nmr.dat %>% group_by(cluster,urban,admin1,admin1.char,admin1.name,strata,v005) %>% 
          reframe(total = sum(total), Y = sum(Y))
      }else{
        nmr.dat <- nmr.dat %>% group_by(cluster,urban,admin1,admin1.char,admin1.name,admin2,admin2.char,admin2.name,strata,v005) %>% 
          reframe(total = sum(total), Y = sum(Y))
      }
      nmr.dat <- as.data.frame(nmr.dat)
      
      save(nmr.dat, file = paste0('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/',country_t,'/',beg.year,'_', survey_year,'_',country_t,'_cluster_dat.rda')) 
      
      # load urban rural weights ----
      adm1_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_",gadm.abbrev,"/admin1_u1_",mid.year,"_urban_frac.rds"))
      adm1_UR_weights <- adm1_UR_weights[order(as.numeric(str_remove(adm1_UR_weights$adm_idx,'admin1_'))),]
      if(!(country_t %in% noAdm2_countries)){
        adm2_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_",gadm.abbrev,"/admin2_u1_",mid.year,"_urban_frac.rds"))
        adm2_UR_weights <- adm2_UR_weights[order(as.numeric(str_remove(adm2_UR_weights$adm_idx,'admin2_'))),]
        
      }
      
      # admin2 weights (don't have for every country)
      # load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/worldpop/adm2_weights_u1.rda"))
      # adm2_weights <- weight.adm2.u1[weight.adm2.u1$years==2014,]
      
      # admin1 weights (don't have for every country)
      #load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/worldpop/adm1_weights_u1.rda"))
      # adm1_weights <- weight.adm1.u1[weight.adm1.u1$years==mid.year,]
      
      # direct estimates -----
      data_for_direct <- nmr.dat
      data_for_direct$age <- 0
      data_for_direct$years <- 'All'
      data_for_direct$died <- data_for_direct$Y
      data_for_direct$v005 <- data_for_direct$v005
      dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                                   regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                                   Ntrials = "total",weightsVar = "v005")
      
      natl.dir <- dir.est[dir.est$region=='All',]
      admin1.dir <- dir.est[dir.est$region!='All',]
      admin1.dir <- admin1.dir %>% dplyr::select(region,mean,lower,upper) %>% rename(admin1_mean = mean, admin1_lower = lower,admin1_upper=upper) %>% rename(admin1.char = region)
      admin1.dir$admin1 <- as.numeric(sapply(1:nrow(admin1.dir),function(i){str_split(admin1.dir$admin1.char[i],'_')[[1]][2]}))
      
      # admin 1 fixed effect model-----
      admin1.fit <- INLA::inla(Y ~ urban +  factor(admin1) -1,
                               family = 'binomial',
                               data = nmr.dat,
                               control.family=list(variant=1,link='log'),
                               control.predictor = list(compute = FALSE, link = 1),
                               control.compute = list(config=T),
                               Ntrials = nmr.dat$total)

      admin1.inlafe <- cbind(nmr.dat,est = exp(admin1.fit$summary.linear.predictor$`0.5quant`)) %>% dplyr::select(admin1.char,urban,est) %>% unique()
      admin1.inlafe <- merge(admin1.inlafe,adm1_UR_weights[,2:3],by.x='admin1.char',by.y='adm_idx') %>%
        mutate(urb_frac = ifelse(urban=='rural',1-urb_frac,urb_frac)) %>% group_by(admin1.char) %>% summarise(admin1_inla_median = sum(est*urb_frac))

      
      # admin 1 random effect model-----
      # admin1RE.fit <- INLA::inla(Y ~ urban -1 +  f(admin1, model="bym2", graph = admin1.info$admin.mat, scale.model=T, constr=T,
      #                                            adjust.for.con.comp=T, 
      #                                            hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1))) -1,
      #                          family = 'binomial', 
      #                          data = nmr.dat, 
      #                          control.family=list(variant=1,link='log'),
      #                          control.predictor = list(compute = FALSE, link = 1), 
      #                          control.compute = list(config=T),
      #                          Ntrials = nmr.dat$total)
      # 
      # # sample from posterior
      # cs <- admin1RE.fit$misc$configs$contents$tag
      # cs <- cs[cs != "Predictor"]
      # select <- list()
      # for (i in 1:length(cs)) {
      #   select[[i]] <- 0
      #   names(select)[i] <- cs[i]
      # }
      # sampFull <- INLA::inla.posterior.sample(n = 1000, result = admin1RE.fit, intern = TRUE, selection = select)
      # sampFull.draws <- matrix(NA,1000,length(sampFull[[1]]$latent))
      # for(i in 1:1000){
      #   sampFull.draws[i,] <- sampFull[[i]]$latent
      # }
      # fields <- colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)
      # 
      # admin1RE.inla <- matrix(0,ncol=(2*max(admin.key$admin1)+2),nrow=max(admin.key$admin1))
      # for(area in 1:max(admin.key$admin1)){
      #   admin1RE.inla[area,area] <- 1 #choose admin1 area
      # }
      # admin1RE.inla[,2*max(admin.key$admin1)+1] <- adm1_UR_weights$urb_frac
      # admin1RE.inla[,2*max(admin.key$admin1)+2] <- 1 - adm1_UR_weights$urb_frac
      # 
      # admin1RE.inla <- admin1RE.inla%*%t(sampFull.draws)
      # res.admin1RE <- data.frame(admin1 = unique(admin.key$admin1),admin1.char = unique(admin.key$admin1.char))
      # 
      # res.admin1RE$adm1_median <- (apply(admin1RE.inla,1,function(x){median(exp(x))}))
      # res.admin1RE$adm1_lower95 <- (apply(admin1RE.inla,1,function(x){quantile(exp(x),0.025)}))
      # res.admin1RE$adm1_upper95 <- (apply(admin1RE.inla,1,function(x){quantile(exp(x),0.975)}))
      # 
      # natl.inla <- rep(0,2*max(admin.key$admin1)+2)
      # for(area in 1:max(admin.key$admin1)){
      #   natl.inla[area] <- weight.adm1.u1$proportion[area] #choose admin1 area
      # }
      # natl.inla[2*max(admin.key$admin1)+1] <- sum(adm1_UR_weights$urb_frac)/nrow(adm1_UR_weights)
      # natl.inla[2*max(admin.key$admin1)+2] <- 1 - sum(adm1_UR_weights$urb_frac)/nrow(adm1_UR_weights)
      # natl.inla <- natl.inla%*%t(sampFull.draws)
      # 
      # res.natl <- data.frame(median = median(exp(natl.inla)),
      #                        lower = quantile(exp(natl.inla),0.025),
      #                        upper = quantile(exp(natl.inla),0.975))
      
      # admin 1 FE + admin 2 BYM2 model -----
      if(!(country_t %in% noAdm2_countries)){
        admin2FE.fit <- INLA::inla(Y ~ urban + factor(admin1) + f(admin2, model="bym2", graph = admin2.info$admin.mat, scale.model=T, constr=T,
                                                                  adjust.for.con.comp=T, 
                                                                  hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1))) -1,
                                   family = 'binomial', 
                                   data = nmr.dat, 
                                   control.family=list(variant=1,link='log'),
                                   control.predictor = list(compute = FALSE, link = 1), 
                                   control.compute = list(config=T),
                                   Ntrials = nmr.dat$total)
        
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
        
        admin2FE.inla <- matrix(0,ncol=(max(admin.key$admin1)+2*max(admin.key$admin2)+1),nrow=max(admin.key$admin2))
        for(area in 1:max(admin.key$admin2)){
          admin2FE.inla[area,area] <- 1 #choose admin1 area
          if(admin.key[admin.key$admin2==area,]$admin1>1){
            admin2FE.inla[area,2*max(admin.key$admin2)+1+admin.key[admin.key$admin2==area,]$admin1] <- 1 #choose admin1 area
          }
        }
        admin2FE.inla[,2*max(admin.key$admin2)+1] <- adm2_UR_weights$urb_frac
        admin2FE.inla[,2*max(admin.key$admin2)+2] <- 1 - adm2_UR_weights$urb_frac
        
        admin2FE.inla <- admin2FE.inla%*%t(sampFull.draws)
        admin2.inla2 <- data.frame(admin2 = 1:nrow(admin.key))
        admin2.inla2$admin2.char <- paste0('admin2_', admin2.inla2$admin2)
        
        # take out region in Tanzania (should really do this before modeling -- fix later)
        if(country_t=='Tanzania'){
          admin2FE.inla <- admin2FE.inla[-c(135),]
          admin2.inla2 <- admin2.inla2[-c(135),]
        }
        
        admin2.inla2$adm2_median <- (apply(admin2FE.inla,1,function(x){median(exp(x))}))
        admin2.inla2$adm2_lower95 <- (apply(admin2FE.inla,1,function(x){quantile(exp(x),0.025)}))
        admin2.inla2$adm2_upper95 <- (apply(admin2FE.inla,1,function(x){quantile(exp(x),0.975)}))
        
        admin1.inla2 <- data.frame(admin1 = rep(sort(unique(admin.key$admin1)),each=2),
                                   U = rep(c(0,1),max(admin.key$admin1)))
        admin1.inla2$est <- rep(c(0,admin2FE.fit$summary.fixed$`0.5quant`[3:nrow(admin2FE.fit$summary.fixed)]),each=2) + 
          rep(admin2FE.fit$summary.fixed$`0.5quant`[c(2,1)],max(admin.key$admin1))
        admin1.inla2$admin1.char <- paste0('admin1_',admin1.inla2$admin1)
        admin1.inla2 <- merge(admin1.inla2,adm1_UR_weights,by.x='admin1.char',by.y='adm_idx')
        admin1.inla2[admin1.inla2$U==0,]$urb_frac <- 1-admin1.inla2[admin1.inla2$U==0,]$urb_frac
        admin1.inla2 <- admin1.inla2 %>% group_by(admin1.char,admin1) %>% 
          summarise(adm1_median=exp(sum(est*urb_frac)))
      }
      
      # admin 2 BYM2 model -----
      admin2.fit <- INLA::inla(Y ~ urban + f(admin2, model="bym2", graph = admin2.info$admin.mat, scale.model=T, constr=T,
                                              hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1))) -1,
                              family = 'binomial',
                        data = nmr.dat,
                        control.family=list(variant=1,link='log'),
                        control.predictor = list(compute = FALSE, link = 1),
                        control.compute = list(config=T),
                        Ntrials = nmr.dat$total)

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

      admin2_res_draws <- sampFull.draws[,1:nrow(admin.key)] +
        sampFull.draws[,2*nrow(admin.key)+1]%*%t(adm2_UR_weights$urb_frac) +
        sampFull.draws[,2*nrow(admin.key)+2]%*%t(1-adm2_UR_weights$urb_frac)

      admin2.inla <- data.frame(admin2 = 1:nrow(admin.key))
      admin2.inla$admin2.char <- paste0('admin2_', admin2.inla$admin2)

      admin2.inla$adm2_median <- (apply(admin2_res_draws,2,function(x){median(exp(x))}))
      admin2.inla$adm2_lower95 <- (apply(admin2_res_draws,2,function(x){quantile(exp(x),0.025)}))
      admin2.inla$adm2_upper95 <- (apply(admin2_res_draws,2,function(x){quantile(exp(x),0.975)}))
      
      # compare ----
      
      if(country_t %in% noAdm2_countries){
        all_res <- merge(admin1.dir,admin1.inlafe)
        
        g <- all_res %>% ggplot() + geom_point(aes(admin1_mean,admin1.char,col='1',pch='1'),size=4) +
          geom_point(aes(admin1_lower,admin1.char,col='1'),pch=3,size=4) +
          geom_point(aes(admin1_upper,admin1.char,col='1'),pch=3,size=4) +
          geom_point(aes(admin1_inla_median,admin1.char,col='3',pch='3'),size=4) +
          ggtitle(paste0(country_t,' ',beg.year,'-',survey_year)) + theme_minimal() +
          ylab('') + xlab('NMR') + theme(legend.position = 'bottom') +
          scale_colour_manual(name = '', values =c('1'='blue','3'='darkgreen'),
                              labels = c('Admin 1 Direct','Admin 1 FE')) +
          scale_shape_manual(name = '', values =c('1'=17,'3'=2),
                             labels = c('Admin 1 Direct','Admin 1 FE'))
      }else{
        all_res <- merge(admin.key,admin2.inla2,by=c('admin2.char','admin2'))
        all_res <- merge(all_res,admin2.inla,by=c('admin2.char','admin2'))
        all_res <- merge(all_res,admin1.inla2,by=c('admin1','admin1.char'))
        all_res <- merge(all_res,admin1.inlafe)
        all_res <- merge(all_res,admin1.dir)
        #all_res <- merge(all_res,adm2_weights[,1:2],by.x='admin2.char',by.y='region')
        #all_res <- merge(all_res,all_res %>% group_by(admin1.char,admin1,admin1_mean,admin1_inla_median) %>% summarise(adm1_median_agg = sum(proportion*adm2_median)/sum(proportion)))
        
        g <- all_res %>% ggplot() + geom_point(aes(admin1_mean,admin1.char,col='1',pch='1'),size=3) +
          geom_point(aes(admin1_lower,admin1.char,col='1'),pch=3,size=3) +
          geom_point(aes(admin1_upper,admin1.char,col='1'),pch=3,size=3) +
          geom_point(aes(adm2_median,admin1.char,col='2',pch='2'),size=3) + 
          geom_point(aes(adm1_median,admin1.char,col='3',pch='3'),size=3) +
          #geom_point(aes(admin1_inla_median,admin1.char,col='3',pch='3'),size=3) +
          ggtitle(paste0(country_t,' ',beg.year,'-',survey_year)) + theme_minimal() +
          ylab('') + xlab('NMR') + theme(legend.position = 'bottom') +
          scale_colour_manual(name = '', values =c('1'='blue','2'='orange','3'='darkgreen'),
                              labels = c('Admin 1 Direct','Admin 2 BYM2 (with Admin 1 FE)','Admin 1 FE')) +
          scale_shape_manual(name = '', values =c('1'=17,'2'=1,'3'=2),
                             labels = c('Admin 1 Direct','Admin 2 BYM2 (with Admin 1 FE)','Admin 1 FE'))
      }
      
      plot_list[[length(plot_list)+1]] <- g
      
    }
  }
  pdf(file = paste0(country_t,' surveys.pdf'))
  if(length(unique(admin.key$admin1.name))<=10){
    print(ggarrange(plotlist = plot_list, common.legend = T, nrow=3,ncol=1))
  }else if(length(unique(admin.key$admin1.name))<=20){
    print(ggarrange(plotlist = plot_list, common.legend = T, nrow=2,ncol=1))
  }else{
    print(ggarrange(plotlist = plot_list, common.legend = T, nrow=1,ncol=1))
  }
  dev.off()
}

