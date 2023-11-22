# function to get admin1 and national aggregated estimates from INLA model
getAggregated <- function(inla_mod,data,nsim=1000,
                          admin1.weights=NULL,admin1.UR.weights=NULL,
                          admin2.to.admin1.weights = NULL, admin2.weights = NULL){
  
  #put weights in right order
  if(!is.null(admin1.weights)){
    admin1.weights <- admin1.weights[order(as.numeric(str_remove(admin1.weights$region,'admin1_'))),]
    }else{admin1.weights = data.frame(region = unique(data$admin1.char), proportion = 0)}
  if(!is.null(admin1.UR.weights)){
    admin1.UR.weights <- admin1.UR.weights[order(as.numeric(str_remove(admin1.UR.weights$region,'admin1_'))),]
  }else{admin1.UR.weights = data.frame(region = unique(data$admin1.char),urban = 0,rural=0)}
  if(!is.null(admin2.weights))
    admin2.weights <- admin2.weights[order(as.numeric(str_remove(admin2.weights$region,'admin2_'))),]
  if(!is.null(admin2.to.admin1.weights))
    admin2.to.admin1.weights <- admin2.to.admin1.weights[order(as.numeric(str_remove(admin2.to.admin1.weights$region,'admin2_'))),]
  
  # data frame for admin 1 overall estimates
  outFullOverall <- data.frame(admin1 = 1:max(data$admin1))
  outFullOverall$admin1.char <- paste0('admin1_', outFullOverall$admin1)
  
  # sample from posterior
  cs <- inla_mod$misc$configs$contents$tag
  cs <- cs[cs != "Predictor"]
  select <- list()
  for (i in 1:length(cs)) {
    select[[i]] <- 0
    names(select)[i] <- cs[i]
  }
  
  sampFull <- INLA::inla.posterior.sample(n = nsim, result = inla_mod, intern = TRUE, selection = select)
  sampFull.draws <- matrix(NA,nsim,length(sampFull[[1]]$latent))
  for(i in 1:nsim){
    sampFull.draws[i,] <- sampFull[[i]]$latent
  }
  fields <- colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)
  
  ## matrix which specifies which columns to include and how to weight them for admin1 
  AA.wt.adm1 <- matrix(0,nrow = length(fields),ncol = max(data$admin1))
  # matrix which specifies which columns to include and how to weight them for national
  AA.wt.natl <- rep(0,ncol(sampFull.draws))
  
  for(k in 1:max(data$admin1)){
    # add admin1 effects (if in model)
    col.id <- c(which(paste0("admin1:",k) == fields),which(paste0("factor(admin1)",k,":1") == fields))
    if(length(col.id)==1){
      AA.wt.adm1[col.id,] <- I(outFullOverall$admin1==k)
      AA.wt.natl[col.id] <- admin1.weights$proportion[k]
    }
    
    # add weights for admin1.strata effects (if in model)
    urban.id <- which(paste0("factor(admin1.strata)",k,":1") == fields)
    if(length(urban.id)==1){
      AA.wt.adm1[urban.id,] <- I(outFullOverall$admin1==k)*admin1.UR.weights$urban[k]
      AA.wt.natl[urban.id] <- admin1.weights$proportion[k]*admin1.UR.weights$urban[k]
    }
    rural.id <- which(paste0("factor(admin1.strata)",k+max(data$admin1),":1") == fields)
    if(length(rural.id)==1){
      AA.wt.adm1[rural.id,] <- I(outFullOverall$admin1==k)*admin1.UR.weights$rural[k]
      AA.wt.natl[rural.id] <- admin1.weights$proportion[k]*admin1.UR.weights$rural[k]
    }
  }
  # add weights for urban/rural (if in model)
  rural.id <- c(which("factor(strata)0:1" == fields),which("factor(urban)rural:1" == fields))
  if(length(rural.id)==1){
    AA.wt.adm1[rural.id,] <- admin1.UR.weights$rural
    AA.wt.natl[rural.id] <- sum(admin1.UR.weights$rural*admin1.weights$proportion)
  }
  urban.id <- c(which("factor(strata)1:1" == fields),which("factor(urban)urban:1" == fields))
  if(length(urban.id)==1){
    AA.wt.adm1[urban.id,] <- admin1.UR.weights$urban
    AA.wt.natl[urban.id] <- sum(admin1.UR.weights$urban*admin1.weights$proportion)
  }
  # add weights for any admin2 effects (if in model)
  if('admin2' %in% colnames(data)){
    for(k in 1:max(data$admin2)){
      col.id <- c(which(paste0("admin2:",k) == fields),which(paste0("factor(admin2)",k,":1") == fields))
      if(length(col.id)==1){
        adm1.area.t <- as.numeric(str_sub(admin2.to.admin1.weights$admin1.char[k],-1))
        AA.wt.adm1[col.id,adm1.area.t] <- admin2.to.admin1.weights$proportion[k]
        AA.wt.natl[col.id] <- admin2.weights$proportion[k]
      }
    }
  }
  
  ## admin1 estimates
  AA.adm1 <- sampFull.draws %*% AA.wt.adm1
  outFullOverall$log.mean <- colMeans(AA.adm1)
  outFullOverall$mean <- colMeans(exp(AA.adm1))
  outFullOverall$median <- colMedians(exp(AA.adm1))
  outFullOverall$variance <- colVars(exp(AA.adm1))
  outFullOverall$lower <- apply(exp(AA.adm1),2,quantile,0.05)
  outFullOverall$upper <- apply(exp(AA.adm1),2,quantile,0.95)
  
  ## national estimates
  AA.natl <- sampFull.draws %*% AA.wt.natl
  outFullNatl <- data.frame(mean = colMeans(exp(AA.natl)))
  outFullNatl$log.mean <- colMeans(AA.natl)
  outFullNatl$median <- colMedians(exp(AA.natl))
  outFullNatl$variance <- colVars(exp(AA.natl))
  outFullNatl$lower <- apply(exp(AA.natl),2,quantile,0.05)
  outFullNatl$upper <- apply(exp(AA.natl),2,quantile,0.95)
  
  return(list(adm1.est = outFullOverall, natl.est = outFullNatl))
}
