library(raster)
library(rgdal)
library(terra)

data.dir <- "/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Sierra_Leone"

#load data
setwd(data.dir)
load(paste0('Sierra_Leone_cluster_dat_1frame.rda'),envir = .GlobalEnv)
mod.dat <- mod.dat[mod.dat$years==2015,]
cluster_list<-mod.dat[!duplicated(mod.dat[c('cluster','survey','LONGNUM','LATNUM')]),]

# download population surface
url <- "https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/2016/SLE/sle_ppp_2016_1km_Aggregated_UNadj.tif"
download.file(url, destfile='Population/sle_ppp_2016_1km_Aggregated_UNadj.tif', method = "libcurl",mode="wb")

# UNadjusted population counts
worldpop <- rast('Population/sle_ppp_2016_1km_Aggregated_UNadj.tif')
info.cells <- (1:nrow(coords))[as.vector(!is.nan(values(worldpop)))]
info.coords <- xyFromCell(worldpop,info.cells)
# longitudes in ordered vector
longs <- unique(sort(info.coords[,1]))
lats <- unique(sort(info.coords[,2]))

# for each urban cluster, average population within 2km distance
urban_cluster_list <- cluster_list[cluster_list$urban=='urban',]
urban_cluster_list$pop_area_est <- sapply(1:nrow(urban_cluster_list),function(i){
  #find closest long and lat
  closest_coord_id <- info.coords[which.min((urban_cluster_list[i,]$LONGNUM - info.coords[,1])^2 + (urban_cluster_list[i,]$LATNUM - info.coords[,2])^2),]
  #make circle of 2km radius around closest coordinates -- (ends up being rotated 2x2 square)
  long_neighborhood <- longs[(which(longs==closest_coord_id[1])-2):(which(longs==closest_coord_id[1])+2)]
  lat_neighborhood <- lats[(which(lats==closest_coord_id[2])-2):(which(lats==closest_coord_id[2])+2)]
  closest_coords <- rbind(as.matrix(expand.grid(x=long_neighborhood[2:4],y=lat_neighborhood[2:4])),
                          cbind(long_neighborhood[c(1,5)],lat_neighborhood[3]),cbind(long_neighborhood[3],lat_neighborhood[c(1,5)]))
  #get the population values at these coordinates
  cell_ids <- cellFromXY(worldpop,closest_coords)
  pops <- values(worldpop)[cell_ids]
  #return the average value over the 3x3 area, excluding any cells that don't have a count (b/c they are outside boundary)
  return(mean(pops,na.rm=T))
})

# for each rural cluster, average population within 2km distance
rural_cluster_list <- cluster_list[cluster_list$urban=='rural',]
rural_cluster_list$pop_area_est <- sapply(1:nrow(rural_cluster_list),function(i){
  #find closest long and lat
  closest_coord_id <- info.coords[which.min((rural_cluster_list[i,]$LONGNUM - info.coords[,1])^2 + (rural_cluster_list[i,]$LATNUM - info.coords[,2])^2),]
  #make 10 km radius circle around closest coordinate -- (ends up being rotated 5x5 square) CODE
  long_neighborhood <- longs[(which(longs==closest_coord_id[1])-5):(which(longs==closest_coord_id[1])+5)]
  lat_neighborhood <- lats[(which(lats==closest_coord_id[2])-5):(which(lats==closest_coord_id[2])+5)]
  #make circle of 5km radius around closest coordinates
  closest_coords <- rbind(as.matrix(expand.grid(x=long_neighborhood[3:9],y=lat_neighborhood[3:9])),
                          cbind(long_neighborhood[c(1,11)],lat_neighborhood[6]),
                          cbind(long_neighborhood[6],lat_neighborhood[c(1,11)]),
                          as.matrix(expand.grid(x=long_neighborhood[c(2,10)],y=lat_neighborhood[4:8])),
                          as.matrix(expand.grid(x=long_neighborhood[4:8],y=lat_neighborhood[c(2,10)])))
  
  #get the population values at these 9 coordinated
  cell_ids <- cellFromXY(worldpop,closest_coords)
  pops <- values(worldpop)[cell_ids]
  #return the average value over the 3x3 area, excluding any cells that don't have a count (b/c they are outside boundary)
  return(mean(pops,na.rm=T))
})

# use that proportion and EA size mean to estimate EA size
urban_cluster_list[order(urban_cluster_list$pop_area_est),]
