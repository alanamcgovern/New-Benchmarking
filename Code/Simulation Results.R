library(ggpubr)
library(tidyverse)
library(scales)
library(sf)
library(patchwork)

# with Exact results
admin1.full.res <- admin2.full.res <- NULL
for(setting_num in c(1:4)){
  load(file = paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/Data/sim_setting_",setting_num,".rda"))
  load(file = paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/Results/sim",setting_num,"results_all_admin1.rda"))
  load(file = paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/Results/sim",setting_num,"results_all_admin2.rda"))
  
  n_admin1 <- length(unique(admin1.res$admin1))
  
  admin1_UR_weights <- sim_setting$all_dat %>% group_by(A1) %>% summarise(urban_prop = sum(U*N)/sum(N))
  admin2_UR_weights <- sim_setting$all_dat %>% group_by(A2) %>% summarise(urban_prop = sum(U*N)/sum(N))
  admin2_weights <- sim_setting$all_dat %>% group_by(A2,A1) %>% summarise(prop = sum(N)/sum(sim_setting$all_dat$N))
  admin2_to_admin1_weights <- admin2_weights %>% group_by(A1) %>% mutate(prop = prop/sum(prop))
  
  ## ADD VARIABLES ------
  admin1.res <- admin1.res %>% mutate(stan.disc = abs(stan.median - dir.mean),
                                      exact.disc = abs(exact.median - dir.mean),
                                      bench.disc = abs(bench.median - dir.mean),
                                      stan.coverage = I(pop_nmr>=stan.lower & pop_nmr<=stan.upper),
                                      exact.coverage = I(pop_nmr>=exact.lower & pop_nmr<=exact.upper),
                                      bench.coverage = I(pop_nmr>=bench.lower & pop_nmr<=bench.upper)) 
  admin1.pop.nmr <- admin1.res %>% dplyr::select(admin1,pop_nmr) %>% unique() %>% arrange(pop_nmr)
  admin1.pop.nmr$admin1.ordered <- 1:nrow(admin1.pop.nmr)
  admin1.res <- merge(admin1.res, admin1.pop.nmr)
  
  admin2.res <- admin2.res %>% mutate(stan.width = stan.upper - stan.lower,
                                      bench.width = bench.upper - bench.lower,
                                      exact.width = exact.upper - exact.lower,
                                     stan.50coverage = I(pop_nmr>=stan.lower50 & pop_nmr<=stan.upper50),
                                      stan.coverage = I(pop_nmr>=stan.lower & pop_nmr<=stan.upper),
                                      exact.coverage = I(pop_nmr>=exact.lower & pop_nmr<=exact.upper),
                                     exact.50coverage = I(pop_nmr>=exact.lower50 & pop_nmr<=exact.upper50),
                                     bench.50coverage = I(pop_nmr>=bench.lower50 & pop_nmr<=bench.upper50),
                                      bench.coverage = I(pop_nmr>=bench.lower & pop_nmr<=bench.upper))
  admin2.pop.nmr <- admin2.res %>% dplyr::select(A2,pop_nmr) %>% unique() %>% arrange(pop_nmr)
  admin2.pop.nmr$admin2.ordered <- 1:nrow(admin2.pop.nmr)
  admin2.res <- merge(admin2.res, admin2.pop.nmr)
  
  admin1.full.res <- rbind(admin1.full.res, data.frame(admin1.res,sim=setting_num))
  admin2.full.res <- rbind(admin2.full.res, data.frame(admin2.res,sim=setting_num))

}

admin1.full.res$sim_label <- ifelse(admin1.full.res$sim==4,'1a',admin1.full.res$sim)
admin2.full.res$sim_label <- ifelse(admin2.full.res$sim==4,'1a',admin2.full.res$sim)

#admin2.full.res %>% mutate(n_deaths_t = if_else(is.na(n_deaths),0,n_deaths)) %>% group_by(sim) %>% summarise(max(n_deaths_t))
#admin2.full.res %>% group_by(sim) %>% summarise(val = mean(is.na(n_deaths))) 
## PLOTS -------

colors <- hue_pal()(4)
stan_color <- colors[1]
bench_color <- colors[3]
exact_color <- colors[2]

disc_plot <-  admin1.full.res %>% group_by(admin1.ordered,sim_label,pop_nmr) %>% 
  reframe(stan = 1000*mean(stan.disc),exact = 1000*mean(exact.disc),
          bench = 1000*mean(bench.disc)) %>% 
  ggplot() + 
  geom_point(aes(x=admin1.ordered,y=stan,color='stan',pch='stan'),size=3) + 
  geom_point(aes(x=admin1.ordered,y=bench,color='bench',pch='bench'),size=3) + 
  geom_point(aes(x=admin1.ordered,y=exact,color='exact',pch='exact'),size=3) + 
  ggtitle('') + 
  ylab('Avg. discrepancy with direct estimate') + 
  scale_x_continuous(name=c('First administrative area, in order of population NMR')) +
  scale_color_manual(name='Model', values=c('bench'=bench_color,'stan'=stan_color,'exact'=exact_color),
                     labels=c('bench'='Soft DABUL','stan'='Standard unit-level (nested)','exact'='Hard DABUL')) +
  scale_shape_manual(name='Model', values=c('bench'=16,'stan'=17,'exact'=3),
                     labels=c('bench'='Soft DABUL','stan'='Standard unit-level (nested)','exact'='Hard DABUL')) +
  facet_grid(~sim_label) +theme_bw() + 
  theme(legend.position = 'bottom',axis.title = element_text(size=12), axis.text = element_text(size=10))


perc_disc_plot <- admin1.full.res %>%
  group_by(admin1,sim_label) %>% reframe(val = (mean(stan.disc - bench.disc))/mean(stan.disc)) %>%
  ggplot() + geom_boxplot(aes(factor(sim_label),100*val)) + xlab('Simulation setting') + 
  ylab('% decrease in discrepancy b/w UL and soft DABUL') + 
  theme_bw() +theme(axis.title = element_text(size=12), axis.text = element_text(size=10))

ggarrange(plotlist = list(disc_plot,perc_disc_plot),nrow=2,heights = c(3,2))

#what about when discrepancy between direct and stan is larger?
disc_plot2 <- admin1.full.res %>% filter(stan.disc>0.001) %>% group_by(admin1.ordered,sim_label,pop_nmr) %>% 
  reframe(stan = 1000*mean(stan.disc),exact = 1000*mean(exact.disc),
          bench = 1000*mean(bench.disc)) %>% 
  ggplot() + 
  geom_point(aes(x=admin1.ordered,y=stan,color='stan',pch='stan'),size=3) + 
  geom_point(aes(x=admin1.ordered,y=bench,color='bench',pch='bench'),size=3) + 
  geom_point(aes(x=admin1.ordered,y=exact,color='exact',pch='exact'),size=3) + 
  ggtitle('') + 
  ylab('Avg. discrepancy with direct estimate') + 
  scale_x_continuous(name=c('First administrative area, in order of population NMR')) +
  scale_color_manual(name='Model', values=c('bench'=bench_color,'stan'=stan_color,'exact'=exact_color),
                     labels=c('bench'='Soft DABUL','stan'='Standard unit-level (nested)','exact'='Hard DABUL')) +
  scale_shape_manual(name='Model', values=c('bench'=16,'stan'=17,'exact'=3),
                     labels=c('bench'='Soft DABUL','stan'='Standard unit-level (nested)','exact'='Hard DABUL')) +
  facet_grid(~sim_label) +theme_bw() + ylim(c(0,2.5)) +
  theme(legend.position = 'bottom',
        axis.title = element_text(size=12),axis.text = element_text(size=10))

perc_disc_plot2 <-  admin1.full.res %>% filter(stan.disc>0.001) %>% 
  group_by(admin1,sim_label) %>% reframe(val = (mean(stan.disc - bench.disc))/mean(stan.disc)) %>%
  ggplot() + geom_boxplot(aes(sim_label,100*val)) + xlab('Simulation setting') + 
  ylab('% decrease in discrepancy b/w UL and soft DABUL') + 
  theme_bw() +theme(axis.title = element_text(size=12),axis.text = element_text(size=10))

ggarrange(plotlist = list(disc_plot2,perc_disc_plot2), nrow=2,heights = c(3,2))

mae_plot <-  admin2.full.res %>% group_by(admin2.ordered,sim_label) %>% 
  reframe(stan = (mean(abs((1000*stan.median)-(1000*pop_nmr)))),
          exact = (mean(abs((1000*exact.median)-(1000*pop_nmr)))),
          bench = (mean(abs((1000*bench.median)-(1000*pop_nmr))))) %>%
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + 
  geom_line(aes(x=admin2.ordered,y=value,col=name,lty=name)) +
  scale_x_continuous(name=c('Second administrative area, in order of population NMR')) +
  facet_grid(~sim_label) +theme_bw() +theme(axis.title = element_text(size=12)) +
  scale_color_manual(name='Model', 
                     labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)'),
                     values=c(bench_color,exact_color,stan_color)) +
  scale_linetype_discrete(name='Model', labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)')) +
  #ggtitle('Avg. Squared Error of Admin 2 estimates across 500 simulations') + 
  ylab('Avg. absolute error')

mae_boxplot <- admin2.full.res %>% group_by(A2, sim_label) %>% 
  reframe(stan = (mean(abs((1000*stan.median)-(1000*pop_nmr)))),
          exact = (mean(abs((1000*exact.median)-(1000*pop_nmr)))),
          bench = (mean(abs((1000*bench.median)-(1000*pop_nmr))))) %>%
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + geom_boxplot(aes(x=factor(sim_label),fill=name,y=value)) + 
  scale_fill_manual(name='Model',
                    labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)'),
                    values=c(bench_color,exact_color,stan_color)) +
  xlab('Simulation setting') +theme_bw() +theme(axis.title = element_text(size=12)) +
  ylab('Avg. absolute error per area')

ggarrange(plotlist = list(mae_plot,mae_boxplot),nrow=2,heights = c(2,3),common.legend = T)

coefvar_plot <-  admin2.full.res %>% group_by(admin2.ordered,sim_label,pop_nmr) %>% 
 reframe(stan = mean(stan.sd/stan.median),
         exact=mean(exact.sd/exact.median),
         bench = mean(bench.sd/bench.median)) %>% 
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + 
  geom_line(aes(x=admin2.ordered,y=value,col=name,lty=name)) +
  scale_x_continuous(name=c('Second administrative area, in order of population NMR')) +
  facet_grid(~sim_label) +theme_bw() +
  scale_color_manual(name='Model', 
                     labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)'),
                     values=c(bench_color,exact_color,stan_color)) +
  scale_linetype_discrete(name='Model', labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)')) +
  ylab('Avg. coefficient of variation') +
  theme(legend.position = 'bottom',axis.title = element_text(size=12))

coefvar_boxplot <- admin2.full.res %>% group_by(A2, sim_label) %>% 
  reframe(stan = mean(stan.sd/stan.median),
          exact=mean(exact.sd/exact.median),
          bench = mean(bench.sd/bench.median)) %>% 
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + geom_boxplot(aes(x=factor(sim_label),fill=name,y=value)) + 
  scale_fill_manual(name='Model',
                    labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)'),
                    values=c(bench_color,exact_color,stan_color)) +
  xlab('Simulation setting') +theme_bw() +
  theme(axis.title = element_text(size=12)) +
  ylab('Avg. coefficient of variation per area') 

ggarrange(plotlist = list(coefvar_plot,coefvar_boxplot),nrow=2,heights = c(2,3),common.legend = T)

cov_plot <-  admin2.full.res %>% group_by(admin2.ordered,sim_label,pop_nmr) %>% 
  reframe(stan = mean(stan.coverage),
          exact=mean(exact.coverage),
          bench = mean(bench.coverage)) %>% 
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + 
  geom_line(aes(x=admin2.ordered,y=value,col=name,lty=name)) +
  geom_hline(yintercept = 0.9,color='sienna') +
  scale_x_continuous(name=c('Second administrative area, in order of population NMR')) +
  facet_grid(~sim_label) + theme_bw() +
  scale_color_manual(name='Model', 
                     labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)'),
                     values=c(bench_color,exact_color,stan_color)) +
  scale_linetype_discrete(name='Model', labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)')) + 
  ylab('90% credible interval coverage') + theme(legend.position = 'bottom',axis.title = element_text(size=12))

cov_boxplot <- admin2.full.res %>% 
  group_by(A2, sim_label) %>% 
  reframe(stan = mean(stan.coverage),
          exact=mean(exact.coverage),
          bench= mean(bench.coverage)) %>% 
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + geom_boxplot(aes(x=factor(sim_label),fill=name,y=value)) + 
  geom_hline(yintercept = 0.9,color='sienna') +
  scale_fill_manual(name='Model',
                    labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)'),
                    values=c(bench_color,exact_color,stan_color)) +
  xlab('Simulation setting') +theme_bw() + 
  theme(legend.position = 'none',axis.title = element_text(size=12)) +
  ylab('90% credible interval coverage per area') 

ggarrange(plotlist = list(cov_plot,cov_boxplot), nrow=2,heights = c(2,3),common.legend = T)

ciwidth_plot <-  admin2.full.res %>% group_by(admin2.ordered,sim_label,pop_nmr) %>% 
  reframe(stan = mean(1000*stan.width),
           exact=mean(1000*exact.width),
          bench = mean(1000*bench.width)) %>% 
  pivot_longer(cols=c('stan',
                         'exact',
                      'bench')) %>% 
  ggplot() + 
  geom_line(aes(x=admin2.ordered,y=value,col=name,lty=name)) +
  scale_x_continuous(name=c('Second administrative area, in order of population NMR')) +
  facet_grid(~sim_label) + theme_bw() +
  scale_color_manual(name='Model', 
                     labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)'),
                     values=c(bench_color,exact_color,stan_color)) +
  scale_linetype_discrete(name='Model', labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)')) + 
  ylab('Avg. 90% credible interval width') + theme(legend.position = 'bottom',axis.title = element_text(size=12))

ciwidth_boxplot <- admin2.full.res %>% 
  group_by(A2, sim_label) %>% 
  reframe(stan = mean(1000*stan.width),
           exact=mean(1000*exact.width),
          bench= mean(1000*bench.width)) %>% 
  pivot_longer(cols=c('stan',
                        'exact',
                      'bench')) %>% 
  ggplot() + geom_boxplot(aes(x=factor(sim_label),fill=name,y=value)) + 
  scale_fill_manual(name='Model',
                    labels=c('Soft DABUL','Hard DABUL','Standard unit-level (nested)'),
                    values=c(bench_color,exact_color,stan_color)) +
  xlab('Simulation setting') +theme_bw() + 
  theme(legend.position = 'none',axis.title = element_text(size=12)) +
  ylab('Avg. 90% credible interval width per area') 

ggarrange(plotlist = list(ciwidth_plot,ciwidth_boxplot), nrow=2,heights = c(2,3),common.legend = T)



## 50% coverage

cov50_plot <-  admin2.full.res %>% group_by(admin2.ordered,sim_label,pop_nmr) %>% 
  reframe(stan = mean(stan.50coverage),
           exact=mean(exact.50coverage),
          bench = mean(bench.50coverage)) %>% 
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + 
  geom_line(aes(x=admin2.ordered,y=value,col=name,lty=name)) +
  geom_hline(yintercept = 0.5,color='sienna') +
  scale_x_continuous(name=c('Second administrative area, in order of population NMR')) +
  facet_grid(~sim_label) + theme_bw() +
  scale_color_manual(name='Model', 
                     labels=c('DABUL','Standard unit-level (nested)','Exact DABUL'),
                     values=c(bench_color,stan_color,exact_color)) +
  scale_linetype_discrete(name='Model', labels=c('DABUL','Standard unit-level (nested)','Exact DABUL')) + 
  ylab('50% credible interval coverage') + theme(legend.position = 'bottom',axis.title = element_text(size=12))

cov50_boxplot <- admin2.full.res %>% 
  group_by(A2, sim_label) %>% 
  reframe(stan = mean(stan.50coverage),
          exact=mean(exact.50coverage),
          bench= mean(bench.50coverage)) %>% 
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + geom_boxplot(aes(x=factor(sim_label),fill=name,y=value)) + 
  geom_hline(yintercept = 0.5,color='sienna') +
  scale_fill_manual(name='Model',
                    labels=c('DABUL','Standard unit-level (nested)','Exact DABUL'),
                    values=c(bench_color,stan_color,exact_color)) +
  xlab('Simulation setting') +theme_bw() + 
  theme(legend.position = 'none',axis.title = element_text(size=12)) +
  ylab('50% credible interval coverage per area') 

ggarrange(plotlist = list(cov50_plot,cov50_boxplot), nrow=2,heights = c(2,3),common.legend = T)

## map plots -----

poly.adm1 <- read_sf(dsn = file.path("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_AGO_shp/gadm41_AGO_1.shp"))
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

poly.adm1 <- merge(poly.adm1,unique.array(data.frame(NAME_1=poly.adm2$NAME_1,A1=poly.adm2$A1)))
#only keep one from each combined area
poly.adm1 <- poly.adm1 %>% filter(NAME_1 %in% c('Uíge','Bengo','Cuanza Sul','Huíla','Cuando Cubango','Bié','Lunda Norte','Moxico'))

ggplot() + geom_sf(data=poly.adm2,aes(fill=factor(A1))) + theme_bw() +
  geom_sf_label(data=poly.adm1,aes(label=A1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),legend.position = 'none') +
  ggtitle('Administrative areas of generated surface')


##########################################

admin1.full.res %>% 
  group_by(admin1, sim_label) %>% 
  reframe(stan = mean(stan.coverage),
          exact=mean(exact.coverage),
          bench= mean(bench.coverage)) %>% 
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + geom_boxplot(aes(x=factor(sim_label),fill=name,y=value)) + 
  geom_hline(yintercept = 0.9,color='sienna') +
  scale_fill_manual(name='Model',
                    labels=c('DABUL','Exact DABUL','Standard unit-level (nested)'),
                    values=c(bench_color,exact_color,stan_color)) +
  xlab('Simulation setting') +theme_bw() + 
  theme(legend.position = 'none',axis.title = element_text(size=12)) +
  ylab('90% credible interval coverage per area') 

admin1.full.res %>% group_by(admin1, sim_label) %>% 
  reframe(stan = mean(stan.sd/stan.median),
          exact=mean(exact.sd/exact.median),
          bench = mean(bench.sd/bench.median)) %>% 
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + geom_boxplot(aes(x=factor(sim_label),fill=name,y=value)) + 
  scale_fill_manual(name='Model',
                    labels=c('DABUL','Exact DABUL','Standard unit-level (nested)'),
                    values=c(bench_color,exact_color,stan_color)) +
  xlab('Simulation setting') +theme_bw() +
  theme(axis.title = element_text(size=12)) +
  ylab('Avg. coefficient of variation per area') 

admin1.full.res %>% group_by(admin1, sim_label) %>% 
  reframe(stan = mean(stan.upper - stan.lower),
          exact=mean(exact.upper - exact.lower),
          bench = mean(bench.upper - bench.lower)) %>% 
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + geom_boxplot(aes(x=factor(sim_label),fill=name,y=value)) + 
  scale_fill_manual(name='Model',
                    labels=c('DABUL','Exact DABUL','Standard unit-level (nested)'),
                    values=c(bench_color,exact_color,stan_color)) +
  xlab('Simulation setting') +theme_bw() +
  theme(axis.title = element_text(size=12))


admin1.full.res %>% group_by(admin1, sim_label) %>% 
  reframe(stan = (mean(abs((1000*stan.median)-(1000*pop_nmr)))),
          exact = (mean(abs((1000*exact.median)-(1000*pop_nmr)))),
          bench = (mean(abs((1000*bench.median)-(1000*pop_nmr))))) %>%
  pivot_longer(cols=c('stan',
                      'exact',
                      'bench')) %>% 
  ggplot() + geom_boxplot(aes(x=factor(sim_label),fill=name,y=value)) + 
  scale_fill_manual(name='Model',
                    labels=c('DABUL','Exact DABUL','Standard unit-level (nested)'),
                    values=c(bench_color,exact_color,stan_color)) +
  xlab('Simulation setting') +theme_bw() +theme(axis.title = element_text(size=12)) +
  ylab('Avg. absolute error per area')


