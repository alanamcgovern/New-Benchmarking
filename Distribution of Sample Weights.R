

cluster_num <- 2500 # number of clusters in strata
cluster_samp_num <- 15 # number of clusters sampled
HH_samp_num <- 30 # number of HH sampled per cluster

clusters <- data.frame(id=1:cluster_num)
clusters$M <- rpois(cluster_num,90) # number of households according to original sampling frame
clusters$L <- clusters$M + round(rnorm(cluster_num,5,5)) # actual number of households after listing (assumes some population growth on average)
clusters$select_prob <- clusters$M/sum(clusters$M)
clusters$wt <- clusters$M*sum(clusters$M)/(HH_samp_num*cluster_samp_num*clusters$M)
#ignore normalization for now

# choose clusters to sample
clusters_samp <- clusters[sample(1:cluster_num,cluster_samp_num,prob=clusters$M),]
