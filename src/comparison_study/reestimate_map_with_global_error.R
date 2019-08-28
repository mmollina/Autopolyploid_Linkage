#R CMD BATCH --no-save --no-restore '--args n.sim=1' reestimate_map.R  reetimate_1.Rout &

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
fl<-list.files("~/repos/Autopolyploid_Linkage/src/simulation3/RData", full.names = TRUE)
(fl<-grep("LODs:5", fl, value = TRUE))
(fl<-grep(paste0("nsim:",nsim,"_"), fl, value = TRUE))
prefix<-"~/repos/Autopolyploid_Linkage/src/simulation3/reestimate_with_globar_error/RData/"
for(i1 in 1:length(fl)){
  current.file<-fl[i1]
  load(current.file)
  posfix<-paste0("globerr_",substring(current.file, first = 69))
  map.globerr<-est_full_hmm_with_global_error(input.map = res.mp$res.mds, error = 0.1, tol = 10e-4, verbose = TRUE)
  save(map.globerr, current.file, file = paste0(prefix, posfix))
}