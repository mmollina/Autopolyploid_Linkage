#R CMD BATCH --no-save --no-restore '--args m=6 rec.frac=0.01 n.ind=200 n.mrk=10 n.sim.start=1 n.sim.stop=10 seed.for.config=140 thres.twopt=1 thres.hmm=1 extend.tail=10 tol=10e-3 tol.final=0.01 case=1' simulation_and_analisys.R  res_1.Rout &
require(mappoly)

##Reading arguments
options(echo=TRUE)
arg <- commandArgs(trailingOnly = TRUE)
print(arg)
for(i in 1:length(arg)){
  eval(parse(text=arg[[i]]))
}
rm(arg)

##Simulation Arguments
 # m<-6                    # ploidy level
 # rec.frac<-0.01          # recombination fraction
 # n.ind<-200              # number of individuals
 # n.mrk<-10               # number of markers
 # n.sim.start<-1          # simulation starting point
 # n.sim.stop<-3          # simulation stop point
 # seed.for.config<-140    # seed for simulating the linkage phase configuration
 # thres.twopt<-3          # threshold for two-point analysis
 # extend.tail<-10         # number of markers to be considered when using HMM
 #                         # at the end of a linkage group
 # thres.hmm<-1            # threshold for the HMM analysis
 # tol<-0.1                # tolerance for the HMM
 # tol.final<-0.001        # tolerance for the final map
 # case<-1                 # Case simulated

##Which scenrario should be simulated
if(case==1)
{
  max.ph=0;
  c.h.a.r<-FALSE
}
if(case==2)
{
  max.ph=m/2;
  c.h.a.r<-FALSE
}
if(case==3)
{
  max.ph=m/2;
  c.h.a.r<-TRUE
}

##Simulating homologous chromosomes
ph.temp<-sim_homologous(m=m,
                        n.mrk=n.mrk,
                        max.d=m/2,    #max dosage number
                        max.ph=max.ph,
                        choose.hom.at.random=c.h.a.r,
                        seed=seed.for.config)

#Recombination fraction vector
rf.vec<-rep(rec.frac, (n.mrk-1))

#dummy simulation just to cache genotype counts
dat <- poly_cross_simulate(m = m,
                           rf.vec = rf.vec,
                           n.mrk = n.mrk,
                           n.ind = n.ind,
                           hom.allele = ph.temp)
dummy_sim <- make_seq_mappoly(dat, "all")
counts <- cache_counts_twopt(dummy_sim, get.from.web=TRUE)

results<-vector("list", length(n.sim.start:n.sim.stop))
##simulation loop
for(i in n.sim.start:n.sim.stop)
{
  ##Simulating full-sib population
  dat<-poly_cross_simulate(m = m,
                           rf.vec = rf.vec,
                           n.mrk = n.mrk,
                           n.ind = n.ind,
                           hom.allele = ph.temp,
                           seed = i,
                           draw = TRUE,
                           file = paste0("config_ploidy_",m , "_config_", seed.for.config,"_case_", case, ".pdf"))
  start<-proc.time()
  ##Making a sequence with all markers
  all.dat<-make_seq_mappoly(dat, arg="all")
  ##Calculating the pairwise recombination fraction for all markers
  all.pairs<-est_pairwise_rf(input.seq=all.dat,
                             counts,
                             n.clusters=1,
                             verbose=FALSE)
  ##Map reconstriction including estimation of linkage phases
  res<-est_rf_hmm_sequential(input.seq=all.dat, thres.twopt=thres.twopt,
                             thres.hmm=thres.hmm, extend.tail=extend.tail,
                             twopt=all.pairs, verbose=TRUE, tol=tol,
                             tol.final=tol.final)
  end<-proc.time()
  ##loglike for each map
  l<-sapply(res$maps, function(x) x$loglike)
  LOD<-max(l)-l
  best<-which(abs(LOD) < 0.1)
  ##Comparing the simulated and estimated linakge phase configuration
  ph.sim.P<-ph.temp$hom.allele.p
  ph.sim.Q<-ph.temp$hom.allele.q
  Q.result<-P.result<-vector("list", length(best))
  for(j in 1:length(best))
  {
    ph.est.P<-res$maps[[best[j]]]$seq.ph$P
    P.result[[j]]<-compare_haplotypes(m, ph.sim.P, ph.est.P)
    ph.est.Q<-res$maps[[best[j]]]$seq.ph$Q
    Q.result[[j]]<-compare_haplotypes(m, ph.sim.Q, ph.est.Q)
  }

  P<-sapply(P.result, function(x) sum(as.logical(x$haplo.ord), na.rm=TRUE))
  Q<-sapply(Q.result, function(x) sum(as.logical(x$haplo.ord), na.rm=TRUE))
  best<-which.max(P+Q)
  ##saving results
  results[[i]]<-list(start.proc.time = start,
                     end.proc.time = end,
                     P.result = P.result[[best]],
                     Q.result = Q.result[[best]],
                     rf.vec = res$maps[[best]]$seq.rf,
                     maps = res,
                     dat = dat)
  save(list=c("results", "dat"), file=paste0("result_ploidy_",m , "_config_", seed.for.config, "_case_", case, "_from_", n.sim.start, "_to_", n.sim.stop, ".RData"))
}

## Check linkage phases
sapply(results, function(x) x$P.result$is.same.haplo)
sapply(results, function(x) x$Q.result$is.same.haplo)
## Linkage phase OK
plot(results[[1]]$maps)
## Linkage phase not OK
plot(results[[2]]$maps)

