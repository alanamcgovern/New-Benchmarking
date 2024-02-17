
# urban Dodoma, Tanzania 2022
{
  n_clusters <- 621 # total number of clusters *GIVEN*
  M_bar <- 117 # average number of births in each cluster *GIVEN*
  n_clusters_samp <-  5 # number of clusters sampled *GIVEN*
}

# urban Oti, Ghana 2022
{
  n_clusters <- 401 # total number of clusters *GIVEN*
  M_bar <- 220 # average number of births in each cluster *GIVEN*
  n_clusters_samp <-  14 # number of clusters sampled *GIVEN*
}

# rural Volta, Ghana 2022
{
  n_clusters <- 1885 # total number of clusters *GIVEN*
  M_bar <- 147 # average number of births in each cluster *GIVEN*
  n_clusters_samp <-  20 # number of clusters sampled *GIVEN*
}

sd_seq <- c(.05,0.1,.15,.20)

n_iter <- 1000
M_vec <- matrix(NA,n_iter,length(sd_seq))
M_obs_bar <- rep(NA,n_iter)

for(i in 1:length(sd_seq)){
  n_births_sd <- M_bar*sd_seq[i] # between cluster variation of cluster size
  
  for(k in 1:n_iter){
    M = round(rnorm(n_clusters,M_bar,n_births_sd))
    M_obs <- M[sample(1:n_clusters,n_clusters_samp,prob = M)]
    M_obs_bar[k] <- mean(M_obs)
  }
  
  M_vec[,i] <- M_obs_bar
}

par(mfrow=c(3,2))
plot(sd_seq*M_bar,colmeans(M_vec), 
     main = paste0('N = ', n_clusters,', n = ', n_clusters_samp, ', Average cluster size = ', M_bar),ylab='SD')

sapply(1:length(sd_seq), function(x){hist(M_vec[,x],main=paste0('SD = ', sd_seq[x]*M_bar),xlab='Average size of observed clusters',
                                          xlim=c(130,170)
                                          )})

