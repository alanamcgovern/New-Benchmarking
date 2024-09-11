## process results from cluster
suppressMessages({
  library(openxlsx)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly=TRUE)

# which dataset
i = as.numeric(args[1])
# which setting
setting_num = as.numeric(args[2])

# number of chains for each dataset
nchains <- 4
# number of iterations for each chain
niter <- 2000
# burnin
burnin <- 1000

setwd(paste0("Sim",setting_num,'_Exact'))

raw.res <- NULL
if(sum(str_detect(list.files(),paste0('Results',i,'-')))==nchains){
  for(j in 1:nchains){
    raw.res.tmp <- read.xlsx(paste0('Results',i,'-',j,'.xlsx'))
    
    if(unique(raw.res.tmp$chain)!=j | unique(raw.res.tmp$dataset!=i))
      stop('Problem with indexing.')
    # check efficiency
    if(median(raw.res.tmp$tree_depth,na.rm = T)<=4)
      warning(paste0('For dataset ', i,' step size is too large'))
    if(median(raw.res.tmp$tree_depth,na.rm = T)>=9)
      warning(paste0('For dataset ', i,' step size is too small'))
    
    raw.res.tmp$iter.id <- 1:nrow(raw.res.tmp)
    raw.res <- rbind(raw.res,raw.res.tmp)
  }
  
  ## remove burnin and save
  raw.res <- raw.res %>% filter(iter.id>burnin) %>% dplyr::select(-c(iter.id))
  saveRDS(raw.res, file=paste0('Results',i,'.rds'))
  
  for(j in 1:nchains){
    suppressMessages(file.remove(paste0('Results',i,'-',j,'.xlsx')))
  }
}
  

