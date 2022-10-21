setwd("src/comparison_study/")
source("simulation.R")
dat<-simulate(m = 4, map.length = 100, n.ind = 200, n.mrk = 300, 
              n.chr = 4, seed.for.config = 10, seed.for.pop = 20)
dat<-dat$dat.mp
plot(dat)
