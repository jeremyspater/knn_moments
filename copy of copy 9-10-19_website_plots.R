#9-10-19 based on 5-24-19_accuracyplots and 5-22-19_sim
#Point is just to make and plot a binary and a clustered neighborhoods, using colors

#Preamble
library(ggplot2)
library(tidyr)
library(plyr);library(dplyr, warn.conflicts = F)
library(tictoc)
library(stargazer)
rm(list=ls())
s = function(x){summary(factor(x))}
####################################################################################################

####################################################################################################
#Set global parameters
xlim = 100 #neighborhood size x. call this meters. or km if this is city scale??
ylim = 100 #neighborhood size y. call this meters. or km if this is city scale??
ngrid = 10 #number of grid spaces
####################################################################################################


####################################################################################################
#Clustered process

library(spatstat)

#https://rdrr.io/cran/spatstat/man/rPoissonCluster.html
# multitype Neyman-Scott process (each cluster is a multitype process)
#[This is same as Matern cluster process, only with marks?]
nclust2 <- function(x0, y0, radius, n, types=c("yes", "no"), p) {
  X <- runifdisc(n, radius, centre=c(x0, y0))
  #M <- sample(types, n, replace=TRUE, prob = c(p, 1-p)) #this just makes the marks random
  M <- sample(types, 1, replace=TRUE, prob = c(p, 1-p)) %>% rep(n)#this gives everybody in the cluster the same mark
  marks(X) <- M
  return(X)
}

#Nonsegregated version: Type is NOT by cluster
nclust1 <- function(x0, y0, radius, n) {
  X <- runifdisc(n, radius, centre=c(x0, y0))
  return(X)
}

ClusterSeg = function(n_samp = 100, p = 0.3, nclus = 13, ka = 0.05, r = 3){
  kids = rPoissonCluster(kappa = ka, #determines number of clusters
                         expand = 10, #expansion window for parents (so cluster center can be offscreen)
                         rcluster = nclust2, #function for generating clusters, above
                         radius = r,#3, #tried 10, but the clusters overlap too much. Size of cluster
                         n = nclus, #tried 5, but it didn't look very segregated. Number in cluster
                         p = p, #prob of minority
                         win = owin(c(0,100),c(0,100))) #limits of window
  
  kids_df = data.frame(long = kids$x,
                       lat = kids$y,
                       minor = kids$marks) 
  ns = min(nrow(kids_df), n_samp)
  return(kids_df[sample(nrow(kids_df), ns), ]) #return sample of n rows
}

ClusterNonSeg = function(n_samp = 100, p = 0.3, nclus = 13, ka = 0.05, r = 3){
  kids = rPoissonCluster(kappa = ka, #determines number of clusters
                         expand = 10, #expansion window for parents (so cluster center can be offscreen)
                         rcluster = nclust1, #function for generating clusters, above
                         radius = r,#3, Size of cluster
                         n = nclus, #Number in cluster
                         #p = p,
                         win = owin(c(0,100),c(0,100))) #limits of window
  
  marks(kids) = sample( c('yes', 'no'), kids$n, replace=TRUE, prob = c(p, 1-p)) #select marks; not by cluster (because nonseg). Determines majority/minority
  
  kids_df = data.frame(long = kids$x,
                       lat = kids$y,
                       minor = kids$marks) 
  ns = min(nrow(kids_df), n_samp)
  return(kids_df[sample(nrow(kids_df), ns), ]) #return sample of n rows
}

#make one and plot
d = ClusterSeg(ka = 0.01, nclus = 100, r = 10, p = input$pminority, n_samp = 1E6) %>% 
  rename(Ethnicity = minor) %>%
  mutate(Ethnicity = mapvalues(Ethnicity, from = c('yes', 'no'), to = c('Minority','Majority')))

#ggplot(data = d) + geom_point(aes(x = long, y = lat, shape = Ethnicity), size = 2) + theme_bw() +
ggplot(data = d) + geom_point(aes(x = long, y = lat, color = Ethnicity), size = 2) + theme_bw() +
  ggtitle('Cluster segregation example') + 
  theme(text = element_text(size=24), plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(values = c(1, 2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_rect(fill = "grey")) +
  scale_color_manual(values=c("darkmagenta", "green4")) +
  xlim(0,100) + ylim(0,100)

path0 = paste0('/Users/jeremyspater/Dropbox/duke/political economy core/prospectus/methods paper/results_output/',Sys.Date(),'/')
dir.create(path0)
#ggsave(filename = paste0(path0, 'cluster_seg_example_website.png'), height = 200, width = 200, units = 'mm'); rm(path0)


####################################################################################################
#Binary segregation

#Nonsegregated neighborhoods: just a uniform distribution
create_Nonseg = function(n, xlim, ylim, pminority){
  fr = data.frame(long = runif(n,0,xlim),
                  lat = runif(n,0,xlim),
                  minor = rbinom(n,1,pminority) == 1)
  return(fr)
}

#Seg neighborhoods: Assume that majority is uniform across whole space

create_Seg = function(n, xlim, ylim, pminority, conc){ #concentration parameter
  fr = data.frame(minor = rbinom(n,1,pminority) == 1,
                  long = NA,
                  lat = NA,
                  ghe = rbinom(n,1,conc)) #"ghetto" parameter: proportion of each group on the "correct" side
  #lat longs for majority (not concentrated), ie ghe = F
  fr$long[fr$minor == F & fr$ghe == F] = runif(sum(fr$minor == F & fr$ghe == F), 0, xlim) #x: runif from 0 to xlim
  fr$lat[fr$minor == F & fr$ghe == F] = runif(sum(fr$minor == F & fr$ghe == F), 0, ylim) #y: runif from 0 to ylim
  #lat longs for majority (concentrated)
  fr$long[fr$minor == F & fr$ghe == T] = runif(sum(fr$minor == F & fr$ghe == T), min = 0, #x: runif from 0 to ghetto limit
                                               max = (1 - pminority)*xlim) #ghetto limit: proportional to size of majority
  fr$lat[fr$minor == F & fr$ghe == T] = runif(sum(fr$minor == F & fr$ghe == T), 0, ylim) #y: runif from 0 to xlim
  #lat longs for minority (not concentrated)
  fr$long[fr$minor == T & fr$ghe == F] = runif(sum(fr$minor == T & fr$ghe == F), 0, xlim) #x: runif from 0 to xlim
  fr$lat[fr$minor == T & fr$ghe == F] = runif(sum(fr$minor == T & fr$ghe == F), 0, ylim) #y: runif from 0 to ylim
  #lat longs for minority (concentrated)
  fr$long[fr$minor == T & fr$ghe == T] = runif(sum(fr$minor == T & fr$ghe == T),  #x: runif from ghetto edge to xlim
                                               min = (1-pminority)*xlim, max = xlim) #ghetto limit: prop. to size of minority
  fr$lat[fr$minor == T & fr$ghe == T] = runif(sum(fr$minor == T & fr$ghe == T), 0, ylim) #y: runif from 0 to ylim
  return(fr)
}

#make one and plot
d = create_Seg(n = 10000, xlim = 100, ylim = 100, pminority = 0.3, conc = 0.8) %>% rename(Ethnicity = minor) %>%
  mutate(Ethnicity = mapvalues(Ethnicity, from = c('TRUE', 'FALSE'), to = c('Minority','Majority')))

ggplot(data = d) + geom_point(aes(x = long, y = lat, color = Ethnicity), size = 2) + theme_bw() +
  ggtitle('Binary segregation example') + theme(text = element_text(size=24), plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(values = c(1, 2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        panel.background = element_rect(fill = "grey")) +
  scale_color_manual(values=c("darkmagenta", "green4")) +
  xlim(0,100) + ylim(0,100)
path0 = paste0('/Users/jeremyspater/Dropbox/duke/political economy core/prospectus/methods paper/results_output/',Sys.Date(),'/')
dir.create(path0)
#ggsave(filename = paste0(path0, 'binary_seg_example_website.png'), height = 200, width = 200, units = 'mm'); rm(path0)
