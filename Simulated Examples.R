

# choose some country to use as population surface: Angola ----
poly.adm1 <- read_sf(dsn = file.path("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_AGO_shp/gadm41_AGO_1.shp"))
poly.adm1$admin1 <-1:nrow(poly.adm1)

poly.adm2 <- read_sf(dsn = file.path("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_AGO_shp/gadm41_AGO_2.shp"))
admin2.mat <- nb2mat(poly2nb(poly.adm2), zero.policy = TRUE)
colnames(admin2.mat) <- rownames(admin2.mat) <- paste0("admin2_", 1:dim(admin2.mat)[1])
poly.adm2$admin2 <- 1:nrow(admin2.mat)

admin.key <- merge(as.data.frame(poly.adm1[,c('NAME_1','admin1')])[,1:2],
                        as.data.frame(poly.adm2[,c('NAME_1','NAME_2','admin2')])[,1:3])

# adjust adjacency matrix (can try collapsing admin1 areas together, or only using some)-----

{
  admin.key$admin1 <- ifelse(admin.key$NAME_1 %in% c('Zaire','Uíge'),1,
                                    ifelse(admin.key$NAME_1 %in% c('Bengo','Luanda'),2,
                                           ifelse(admin.key$NAME_1 %in% c('Cuanza Norte','Cuanza Sul'),3,
                                                  ifelse(admin.key$NAME_1 %in% c('Benguela','Namibe','Huíla'),4,
                                                         ifelse(admin.key$NAME_1 %in% c('Cunene','Cuando Cubango'),5,
                                                                ifelse(admin.key$NAME_1 %in% c('Huambo','Bié'),6,
                                                                       ifelse(admin.key$NAME_1 %in% c('Lunda Norte','Malanje'),7,
                                                                              ifelse(admin.key$NAME_1 %in% c('Lunda Sul','Moxico'),8, NA))))))))
  admin.key <- admin.key[!is.na(admin.key$admin1),]
  n_admin1 <- length(unique(admin.key$admin1))
}

{
  n_admin1 <- 8
  admin.key <- admin.key[admin.key$admin1 %in% 1:n_admin1,]
}
admin2.mat <- admin2.mat[admin.key$admin2,admin.key$admin2]
admin.key$admin1.char <- paste0('admin1_',admin.key$admin1)
n_admin2 <- length(unique(admin.key$admin2))
admin.key$admin2 <- 1:nrow(admin.key)
admin.key$admin2.char <- paste0('admin2_',admin.key$admin2)
colnames(admin2.mat) <- rownames(admin2.mat) <- admin.key$admin2.char

admin2.mat.nested <- admin2.mat
for(i in 1:nrow(admin2.mat.nested)){
  admin2.mat.nested[i,which(admin.key$admin1!=admin.key$admin1[i])] <- 0
  if(sum(admin2.mat.nested[i,])>0){
    admin2.mat.nested[i,] <- admin2.mat.nested[i,]/sum(admin2.mat.nested[i,])
  }
}

Q.admin2 <- -admin2.mat
Q.admin2 <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})*Q.admin2
diag(Q.admin2) <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})
diag(Q.admin2)[diag(Q.admin2)==0] <- 1
Q_scaled <- INLA::inla.scale.model(Q.admin2, constr=list(A=t(eigen(Q.admin2)$vectors[,eigen(Q.admin2)$values<1e-10]), e=rep(0,sum(eigen(Q.admin2)$values<1e-10))))
Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))

# define parameters ------
d <- 1
alphaU <- -3.7
alphaR <- -4
tau <- 50
phi <- 0.3
b <- as.vector(Rfast::rmvnorm(1,rep(0,n_admin2),1/tau*(diag((1-phi),n_admin2) + phi*Q_scaled_inv)))
#make so that smaller areas have higher rates
# for(i in 1:n_admin1){
#   b[admin.key[admin.key$admin1==i,]$admin2] <- sort(b[admin.key[admin.key$admin1==i,]$admin2])
# }

# simulate survey data -----

# parameters other than directly above and below: 
#number of admin1 areas
#number of admin2 areas
#size of admin2 areas within an admin1 areas

n_clusters_urban <- round(rnorm(n_admin1,2000,100)) #number of urban clusters in each admin1 area
n_clusters_rural <- round(rnorm(n_admin1,5000,100)) #number of rural clusters in each admin1 area
n_births_urban <- 100 #average number of births per urban cluster 
sd_n_births_urban <- 10
n_births_rural <-  150 #average number of births per rural cluster 
sd_n_births_rural <- 15
n_clusters_urban_samp <-  round(0.05*n_clusters_urban) #number of urban clusters sampled in each admin1 area
n_clusters_rural_samp <-  round(0.02*n_clusters_rural) #number of rural clusters sampled in each admin1 area
n_births_urban_samp <- 20 # average number of births sampled from each urban cluster 
sd_births_urban_samp <- 2
n_births_rural_samp <- 25 # average number of births sampled from each rural cluster
sd_n_births_rural_samp <- 2

all_dat <- data.frame(
  # (true) number of births in each cluster
  N = round(c(rnorm(sum(n_clusters_urban),n_births_urban,sd_n_births_urban),
              rnorm(sum(n_clusters_rural),n_births_rural,sd_n_births_rural))),
  # urban or rural strata (U=1 if urban, otherwise 0)
  U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
  # admin area of each cluster
  admin1 = c(unlist(sapply(1:n_admin1,function(i){
    rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_admin1,function(i){rep(i,n_clusters_rural[i])}))))
#assign clusters
all_dat$cluster <- 1:nrow(all_dat)
#assign admin2 -- CAN CHANGE TO VARY SIZES OF ADMIN2 AREAS
all_dat$admin2 <- sapply(1:nrow(all_dat),function(x){
  sample(admin.key[admin.key$admin1==all_dat$admin1[x],]$admin2,1,
         prob = admin.key[admin.key$admin1==all_dat$admin1[x],]$admin2 - min(admin.key[admin.key$admin1==all_dat$admin1[x],]$admin2) + 1)})
#draw deaths
all_dat$Y <- sapply(1:nrow(all_dat),function(i){
  rnbinom(1,
          size = 1/d*all_dat$N[i]*exp(alphaU*all_dat$U[i] + alphaR*(1-all_dat$U[i]) + b[all_dat$admin2[i]] + rnorm(1,0,0.05)),
          prob=1/(1+d))})
### sample clusters
obs_dat <- NULL
for(area in 1:n_admin1){
  # randomly select clusters from urban strata based on size of cluster
  obs_dat <- rbind(obs_dat,all_dat[all_dat$admin1==area & all_dat$U==1,][sample(1:n_clusters_urban[area],n_clusters_urban_samp[area],prob = all_dat[all_dat$U==1 & all_dat$admin1==area,]$N),],
                   # randomly select cluster from rural strata  based on size of cluster
                   all_dat[all_dat$admin1==area & all_dat$U==0,][sample(1:n_clusters_rural[area],n_clusters_rural_samp[area],prob = all_dat[all_dat$U==0& all_dat$admin1==area,]$N),])
}

# sample births
obs_dat$n <- obs_dat$Z <- NA
obs_dat[obs_dat$U==1,]$n <- round(min(rnorm(sum(obs_dat$U),n_births_urban_samp,sd_births_urban_samp),obs_dat[obs_dat$U==1,]$N))
obs_dat[obs_dat$U==0,]$n <- round(min(rnorm(sum(1-obs_dat$U),n_births_rural_samp,sd_n_births_rural_samp),obs_dat[obs_dat$U==0,]$N))
for(i in 1:nrow(obs_dat)){
  obs_dat$Z[i] <- sum(sample(c(rep(1,obs_dat$Y[i]),rep(0,(obs_dat$N[i]-obs_dat$Y[i]))), obs_dat$n[i]))
}

# assign sampling  and population weights -----
#births in strata/(# cluster sampled * #births sampled in cluster)
pop_dat <- all_dat %>% group_by(admin1,U) %>% summarise(N=sum(N))
pop_dat$num_clust_samp <- pop_dat$num_births_samp <- NA
pop_dat$num_clust_samp[pop_dat$U==0] <- n_clusters_rural_samp
pop_dat$num_clust_samp[pop_dat$U==1] <- n_clusters_urban_samp
pop_dat$num_births_samp[pop_dat$U==0] <- n_births_rural_samp
pop_dat$num_births_samp[pop_dat$U==1] <- n_births_urban_samp
pop_dat$wt <- pop_dat$N/(pop_dat$num_births_samp*pop_dat$num_clust_samp)
pop_dat$wt <- pop_dat$wt/100
obs_dat <- merge(obs_dat,pop_dat[c('admin1','U','wt')])
obs_dat <- obs_dat[order(obs_dat$admin1),]
#calculate urban/rural weights
admin1_UR_weights <- all_dat %>% group_by(admin1) %>% summarise(urban_prop = sum(U*N)/sum(N))
admin2_UR_weights <- all_dat %>% group_by(admin2) %>% summarise(urban_prop = sum(U*N)/sum(N))
admin2_weights <- all_dat %>% group_by(admin2,admin1) %>% summarise(prop = sum(N)/sum(all_dat$N))
admin2_to_admin1_weights <- admin2_weights %>% group_by(admin1) %>% mutate(prop = prop/sum(prop))

obs_dat$admin1.char <- paste0('admin1_',obs_dat$admin1)
obs_dat$years <- 'All'
obs_dat$age <- '0'
obs_dat$died <- obs_dat$Z
obs_dat$strata <- paste0(obs_dat$admin1,'.',c('rural','urban')[obs_dat$U+1])

# get direct estimates -----
dir.est <- SUMMER::getDirect(obs_dat, years='All',
                               regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                               Ntrials = "n", weightsVar = "wt", CI=0.9)
dir.est <- dir.est[dir.est$region!='All',]
dir.est <- dir.est[order(as.numeric(str_remove(dir.est$region,'admin1_'))),]

# get admin1 fixed effect + admin2 BYM2 estimates ------
admin2.inla.fit <- INLA::inla(Z ~ factor(U) + factor(admin1) + 
                             f(admin2, model="bym2", graph = admin2.mat.nested, scale.model=T, constr=T,
                               adjust.for.con.comp=T,
                               hyper=list(
                                 phi=list(prior="pc", param=c(0.5, 2/3), initial=-3),
                                 prec=list(prior="pc.prec",param=c(0.2/0.31,0.01),initial = 5))) -1,
                           family = 'nbinomial', 
                           data = obs_dat, 
                           control.family=list(link='log'),
                           control.fixed = list(mean =list("factor(U)0"=-3.5, "factor(U)1"=-3.5, default = 0),
                                                prec=list(default = 1/9)),
                           control.predictor = list(compute = FALSE, link = 1), 
                           control.compute = list(config=T),
                           E = n)

# sample from posterior
cs <- admin2.inla.fit$misc$configs$contents$tag
cs <- cs[cs != "Predictor"]
select <- list()
for (i in 1:length(cs)) {
  select[[i]] <- 0
  names(select)[i] <- cs[i]
}
sampFull <- INLA::inla.posterior.sample(n = 1000, result = admin2.inla.fit, intern = TRUE, selection = select)
sampFull.draws <- matrix(NA,1000,length(sampFull[[1]]$latent))
for(i in 1:1000){
  sampFull.draws[i,] <- sampFull[[i]]$latent
}
fields <- colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)

inla_r_samples <- matrix(0,ncol=length(fields),nrow=n_admin2)
for(area in 1:nrow(inla_r_samples)){
  inla_r_samples[area,area] <- 1 #choose admin2 area
  if(admin.key[admin.key$admin2==area,]$admin1>1){
    inla_r_samples[area,2*n_admin2+1+admin.key[admin.key$admin2==area,]$admin1] <- 1 #choose admin1 area
  }
}
inla_r_samples[,2*n_admin2+1] <- 1 - admin2_UR_weights$urban_prop
inla_r_samples[,2*n_admin2+2] <- admin2_UR_weights$urban_prop

inla_r_samples <- inla_r_samples%*%t(sampFull.draws)
admin2.inla <- data.frame(admin2 = 1:nrow(admin.key))
admin2.inla$admin2.char <- paste0('admin2_', admin2.inla$admin2)

admin2.inla$median <- (apply(inla_r_samples,1,function(x){median(exp(x))}))
admin2.inla$lower95 <- (apply(inla_r_samples,1,function(x){quantile(exp(x),0.025)}))
admin2.inla$upper95 <- (apply(inla_r_samples,1,function(x){quantile(exp(x),0.975)}))
admin2.inla <- merge(admin2.inla,admin.key)

admin1.agg <- merge(admin2.inla,admin2_to_admin1_weights) %>% group_by(admin1.char) %>% summarise(median = sum(prop*median))

dir.est %>% ggplot() + geom_point(aes(mean,region),pch=17,size=3) + 
  geom_point(aes(lower,region),pch=3,size=3) + geom_point(aes(upper,region),pch=3,size=3) +
  geom_point(data=admin1.agg,aes(median,admin1.char),pch=17,col='red',size=3) +
  geom_point(data=admin2.inla,aes(median,admin1.char),pch=1,col='blue',size=3)



