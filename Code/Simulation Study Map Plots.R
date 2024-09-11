library(sf)
library(tidyverse)
library(ggpubr)

poly.adm2 <- read_sf(dsn = file.path("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_AGO_shp/gadm41_AGO_2.shp"))

#combine admin1 areas
poly.adm2$A1 <- ifelse(poly.adm2$NAME_1 %in% c('Zaire','Uíge'),1,
                         ifelse(poly.adm2$NAME_1 %in% c('Bengo','Luanda'),2,
                                ifelse(poly.adm2$NAME_1 %in% c('Cuanza Norte','Cuanza Sul'),3,
                                       ifelse(poly.adm2$NAME_1 %in% c('Benguela','Namibe','Huíla'),4,
                                              ifelse(poly.adm2$NAME_1 %in% c('Cunene','Cuando Cubango'),5,
                                                     ifelse(poly.adm2$NAME_1 %in% c('Huambo','Bié'),6,
                                                            ifelse(poly.adm2$NAME_1 %in% c('Lunda Norte','Malanje'),7,
                                                                   ifelse(poly.adm2$NAME_1 %in% c('Lunda Sul','Moxico'),8, NA))))))))
poly.adm2 <- poly.adm2[!is.na(poly.adm2$A1),]
poly.adm2$A2 <- 1:nrow(poly.adm2)

ggplot() + geom_sf(data=poly.adm2,aes(fill=factor(A1))) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),legend.position = 'bottom') +
  ggtitle('Administrative areas of generated surface') +
  scale_fill_discrete(name='First administrative area')

maps <- list()
for(i in 1:4){
  load(file=paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim_setting_",i,".rda"))
  urb_prop <- (sim_setting$all_dat %>% group_by(A2) %>% summarise(prop=sum(U*N)/sum(N)))$prop
  poly.adm2$nmr <- exp(urb_prop[poly.adm2$A2]*sim_setting$alphaU + (1-urb_prop[poly.adm2$A2])*sim_setting$alphaR + sim_setting$beta[poly.adm2$A1] + sim_setting$b[poly.adm2$A2])
  #poly.adm2$nmr <- (sim_setting$all_dat %>% group_by(A2) %>% summarise(nmr = sum(Y)/sum(N)))$nmr

  maps[[i]] <- ggplot() + geom_sf(data=poly.adm2,aes(fill=nmr)) +
    geom_sf(fill = "transparent", color = "#ABA300", lwd=0.75, data = poly.adm2 %>% group_by(A1) %>% summarise()) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), axis.title.y=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),legend.position = 'bottom') +
    scico::scale_fill_scico(palette = "vik",name='true NMR')+
    ggtitle(paste0('Setting ',i))
}

ggarrange(plotlist = maps,nrow=2,ncol=2,common.legend = T)
