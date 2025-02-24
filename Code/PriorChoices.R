
hist(exp(rnorm(1000,-3.5,sqrt(2))))

# prior on urban/rural intercept
a1 <- ggplot(data=data.frame(x=exp(rnorm(1e4,-3.5,1)))) + 
  geom_histogram(aes(x),fill='grey',color='black') +
  theme_bw() + ggtitle('SD=1') +
  xlim(c(-.1,1)) + 
  xlab('Prevalence') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

a2 <- ggplot(data=data.frame(x=exp(rnorm(1e4,-3.5,3)))) + 
  geom_histogram(aes(x),fill='grey',color='black') +
  theme_bw() + ggtitle('SD=3') +
  xlim(c(-.1,1)) + 
  xlab('Prevalence') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

a3 <- ggplot(data=data.frame(x=exp(rnorm(1e4,-3.5,10)))) + 
  geom_histogram(aes(x),fill='grey',color='black') +
  theme_bw() + ggtitle('SD=10') +
  xlim(c(-.1,1)) + 
  xlab('Prevalence') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggarrange(plotlist = list(a1,a2,a3),ncol=3)

# prior on beta, conditional on strata intercept
  # prior on fixed effect for each area, conditional 
  # assuming rare events so we have to choose prior variance so that we still recover that
b1 <- ggplot(data=data.frame(x=exp(rnorm(1000,log(0.015),3)))) + 
  geom_histogram(aes(x),fill='grey',color='black') +
  theme_bw() + ggtitle('Overall prevalence = 0.015') +
  xlab('Prevalence') +
  xlim(c(-0.05,1)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

b2 <- ggplot(data=data.frame(x=exp(rnorm(1000,log(0.03),3)))) + 
  geom_histogram(aes(x),fill='grey',color='black') +
  xlab('Prevalence') +
  xlim(c(-0.05,1)) +
  theme_bw() +ggtitle('Overall prevalence = 0.03') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

b3 <- ggplot(data=data.frame(x=exp(rnorm(1000,log(0.045),3)))) + 
  geom_histogram(aes(x),fill='grey',color='black') +
  xlab('Prevalence') +
  xlim(c(-0.05,1)) +
  theme_bw() +ggtitle('Overall prevalence = 0.05') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggarrange(plotlist = list(b1,b2,b3),nrow=3)

# prior on d
ggplot(data=data.frame(x=rexp(1000,1))) + 
  geom_histogram(aes(x),fill='grey',color='black') +
  ggtitle('Prior on overdispersion parameter') +
  theme_bw() + xlab('d') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# cluster variation given rate and overdispersion parameter
d_seq <- c(1,2,5)
r_seq <- c(0.02,0.03,0.04)
samples4plot <- NULL
for(d_t in d_seq){
  for(r_t in r_seq){
    samples4plot <- rbind(samples4plot,
                     data.frame(deaths = rnbinom(1000, size = 1/d_t*100*r_t, prob = 1/(1+d_t)),
                                r=paste0('area prevalence = ', r_t),d=paste0('d=',d_t)))
  }
}

samples4plot %>% ggplot() + 
    geom_histogram(aes(deaths/100),fill='grey',color='black') +
  facet_grid(d~r) +
    theme_bw() + xlab('Prevalence') +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())


par(mfrow=c(2,1))

hist(rexp(1000,1))




