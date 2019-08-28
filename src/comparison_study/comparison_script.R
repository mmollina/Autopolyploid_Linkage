#R CMD BATCH --no-save --no-restore '--args n.sim=10 n.mrk=200 prob.dose=c(0.4,0.4,0.1,0.1) LODth=3' simulation3_script.R  res_10_200_4411.Rout &

##Reading arguments
options(echo=TRUE)
arg <- commandArgs(trailingOnly = TRUE)
print(arg)
for(j in 1:length(arg)){
  eval(parse(text=arg[[j]]))
}
print(arg)
rm(arg)

require(mappoly)
require(polymapR)
setwd("~/repos/Autopolyploid_Linkage/src/comparison_study/")
source("~/repos/Autopolyploid_Linkage/src/comparison_study/simulation.R")
source("~/repos/Autopolyploid_Linkage/src/comparison_study/build_using_polymapR.R")
source("~/repos/Autopolyploid_Linkage/src/comparison_study/build_using_MAPpoly.R")
# n.sim=28
# n.mrk=200
# prob.dose=c(0.4,0.4,0.1,0.1)
# LODth=3

cat(n.sim, "\n")

#creating directories
dir.create("pedsim_files", showWarnings = FALSE)
dir.create("rec_mat", showWarnings = FALSE)

#####
## Data simulation
#####
dt<-simulate(m = 6, 
             map.length = 100, 
             n.ind = 200, 
             n.mrk = n.mrk,
             n.chr = 1, 
             seed.for.config = n.mrk, 
             seed.for.pop = n.sim*n.mrk, 
             prob.dose = prob.dose)
#####
## Building map using polymapR
#####
time.pm<-system.time(res.pm<-try(build_map_polymapR(dat = dt$dat.pm, 
                                                    n.chr = 1, 
                                                    LOD_bridge = LODth, 
                                                    LOD_assign_LG = LODth, 
                                                    LOD_assign_HM = LODth)))
#####
## Building map using HMM-based procedure
#####
dat = dt$dat.mp
time.mp<-system.time(res.mp<-build_map_MAPpoly(dat = dt$dat.mp, 
                                               start.set = 4, 
                                               thres.twopt = LODth, 
                                               thres.hmm = LODth, 
                                               extend.tail = 50, 
                                               phase.number.limit = 20, 
                                               sub.map.size.diff.limit = Inf))
save.image(paste0("res_nsim:", n.sim,
                  "_nmrk:", n.mrk,
                  "_prob:", paste0(prob.dose, collapse = "-"),
                  "_LODs:", LODth,
                  ".RData"))

try(system("mv *.Rout Rout/"))
try(system("mv *.RData RData/"))

