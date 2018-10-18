#R CMD BATCH --no-save --no-restore '--args m=6 dist.cm=0.5 n.ind=200 n.mrk=50 n.sim.start=1 n.sim.stop=1 seed.for.config=140 thres.twopt=10 extend.tail=100 thres.hmm=10 tol=0.1 tol.final=0.001 case=2 ' sequential.R  res_1.Rout &
require(mappoly)

##Reading arguments
 options(echo=TRUE)
 arg <- commandArgs(trailingOnly = TRUE)
 print(arg)
 for(i in 1:length(arg)){
   eval(parse(text=arg[[i]]))
 }
print(arg)
rm(arg)

##Simulation Arguments
 # m<-4                            # ploidy level
 # map.length<-100                 # distance in cM
 # n.ind<-200                      # number of individuals
 # n.mrk<-200                      # number of markers
 # n.sim.start<-1                  # simulation starting point
 # n.sim.stop<-3                   # simulation stop point
 # seed.for.config<-40             # seed for simulating the linkage phase configuration
 # thres.twopt<-3                  # threshold for two-point analysis
 # extend.tail<-50                 # number of markers to be considered when using HMM at the end of a linkage group
 # thres.hmm<-10                   # threshold for the HMM analysis
 # tol<-10e-2                      # tolerance for the intermediate HMMs
 # tol.final<-10e-4                # tolerance for the final map
 # prefPairing <- 0                # preferential pairing
 # quadrivalents <- 0              # quadrivalent formation rate    

 centromere <- 20 # long arm 4 times as long as the short arm
 phase.number.limit = +Inf      # Limit of number of phases to test in each marker addition 
 sub.map.size.diff.limit = +Inf # Limit of map increment in each marker addition
 natural.pairing <- 0 #probabilities of bivalents and quadrivalents are determined by the “quadrivalents” argument

 #Scenario C
 max.ph=m/2;
 c.h.a.r<-TRUE #choose homologous at random

 #####
 ##Simulating homologous chromosomes
 ph.temp<-sim_homologous(m=m,
                         n.mrk=n.mrk,
                         max.d=m/2,    #max dosage number
                         max.ph=max.ph,
                         choose.hom.at.random=c.h.a.r,
                         seed=seed.for.config)
 php<-ph_list_to_matrix(ph.temp$hom.allele.p, m)
 phq<-ph_list_to_matrix(ph.temp$hom.allele.q, m)
 w<-cbind(php,phq)

 #dummy simulation just to cache genotype counts
 dat <- poly_cross_simulate(m = m,
                            rf.vec = rep(.001, (n.mrk-1)),
                            n.mrk = n.mrk,
                            n.ind = n.ind,
                            hom.allele = ph.temp)
 dummy_sim <- make_seq_mappoly(dat, "all")
 counts <- cache_counts_twopt(dummy_sim, get.from.web=TRUE)
 
 results<-vector("list", length(n.sim.start:n.sim.stop))
 #####
 
 set.seed(m*seed.for.config*case*prefPairing*quadrivalents)
 id<-sample(1:10e6)[1:1000]
 ##simulation loop
 for(i in n.sim.start:n.sim.stop)
 {
   ##Simulating full-sib population
   file.prefix<-paste0("sim_m_", m,
                       "seed_", seed.for.config,
                       "case_", case,
                       "_pp_", prefPairing,
                       "_qua_", quadrivalents,
                       "number_", i)
   
   #########PedigreeSim simulation
   
   #Simulation parameters
   write(x = paste("PLOIDY =", m), file = paste0(file.prefix, ".par"))
   write(x = "MAPFUNCTION = HALDANE", file = paste0(file.prefix, ".par"), append = TRUE)
   write(x = "MISSING = NA", file = paste0(file.prefix, ".par"), append = TRUE)
   write(x = "POPTYPE = F1", file = paste0(file.prefix, ".par"), append = TRUE)
   write(x = paste("POPSIZE =", n.ind), file = paste0(file.prefix, ".par"), append = TRUE)
   write(x = paste("NATURALPAIRING =", natural.pairing), file = paste0(file.prefix, ".par"), append = TRUE)
   write(x = paste0("FOUNDERFILE = ", file.prefix, ".gen"), file = paste0(file.prefix, ".par"), append = TRUE)
   write(x = paste0("CHROMFILE = ", file.prefix, ".chrom"), file = paste0(file.prefix, ".par"), append = TRUE)
   write(x = paste0("MAPFILE = ", file.prefix, ".map"), file = paste0(file.prefix, ".par"), append = TRUE)
   write(x = paste0("OUTPUT = ", file.prefix, "_out"), file = paste0(file.prefix, ".par"), append = TRUE)
   write(x = paste("SEED =", id[i]), file = paste0(file.prefix, ".par"), append = TRUE)
   
   #chromosome file
   a<-c("Ch1", map.length,	20.0,		prefPairing,		quadrivalents)
   a<-matrix(a, nrow = 1)
   write.table(a, file =paste0(file.prefix, ".chrom"), quote = FALSE, row.names = FALSE, col.names = c("chromosome",	"length",	"centromere",	"prefPairing", "quadrivalents"))
   
   ## map
   x<-data.frame(paste0("M",1:n.mrk), "Ch1", seq(0, map.length,length.out = n.mrk))
   write.table(x, file =paste0(file.prefix, ".map"), quote = FALSE, row.names = FALSE, col.names = c("marker", "chromosome", "position"))
   
   ## parent genotypes
   y<-data.frame(paste0("M",1:n.mrk), w)
   write.table(y, file =paste0(file.prefix, ".gen"), quote = FALSE, row.names = FALSE, col.names = c("marker", paste0("P1_", 1:m), paste0("P2_", 1:m)))
   
   ## PedigreeSim
   system(paste0("java -jar PedigreeSim.jar ", file.prefix, ".par"))
   
   ############End of simulation
   
   ## Reading data from pedigreesim
   dat.g<-read.table(file = paste0(file.prefix, "_out_genotypes.dat"), header = TRUE, row.names = 1)
   dat.f<-read.table(file = paste0(file.prefix, "_out.hsa"))
   
   ## Adequating data to mappoly format
   data("hexafake")
   dat<-hexafake
   rm(hexafake)
   
   dat$n.ind<-n.ind
   dat$n.mrk<-n.mrk
   dat$m<-m
   dat$ind.names<-paste0("Ind_", 1:n.ind)
   dat$mrk.names<-paste0("M",1:n.mrk)
   dat$dosage.p<-apply(w[,1:m], 1, function(x) sum(x==1))
   dat$dosage.q<-apply(w[,(m+1):(2*m)], 1, function(x) sum(x==1))
   dat$sequence<-NA
   dat$sequence.pos<-NA
   dat$nphen<-0
   dat$phen<-NA
   dat$geno.dose[1:10, 1:10]
   
   G<-NULL
   ct<-1
   for(j in 1:(n.ind+2))
   {
     dat.g[,ct:(ct+m-1)]
     G<-cbind(G, apply(dat.g[,ct:(ct+m-1)], 1, function(x) sum(x==1)))
     ct<-ct+m
   }
   G<-G[,-c(1:2)]
   colnames(G)<-paste0("Ind_", 1:n.ind)
   
   ## Attribute NA to impossible classes
   ## in this case, instead of NA, we use ploidy level + 1 for missing values (mappoly codification)
   for(j in 1:n.mrk)
   {
     expect<-segreg_poly(m = m, dP = dat$dosage.p[j], dQ = dat$dosage.q[j])
     names(expect)<-c(0:m)
     dr<-is.na(match(G[j,], names(expect[expect!=0])))
     if(any(dr)){
       cat(expect, "--->", table(G[j,]) ,"\n")
       G[j,dr] <- m + 1   
     }
   }
   dat$geno.dose<-G
   dat
   
   # move PedrigreeSim files to adequate dir 
   system(paste0("mv ", file.prefix, "* pedsim_files/"))
   
   start<-proc.time()
   ##Making a sequence with all markers
   all.dat<-make_seq_mappoly(dat, arg="all")
   ##Calculating the pairwise recombination fraction for all markers
   all.pairs<-est_pairwise_rf(input.seq=all.dat,
                              counts,
                              n.clusters=1,
                              verbose=FALSE)
   
   mat<-rf_list_to_matrix(all.pairs)
   
   pdf(paste0("mat", file.prefix, ".pdf"))
   plot(mat)
   dev.off()
   
   # move recombination fraction matrix pdf's to adequate dir
   system(paste0("mv *", file.prefix, "*.pdf rec_mat/"))
   
   ##Map reconstruction including estimation of linkage phases
   res<-est_rf_hmm_sequential(input.seq = all.dat,
                              thres.twopt = thres.twopt,
                              thres.hmm = thres.hmm,
                              extend.tail = extend.tail,
                              twopt = all.pairs,
                              verbose = TRUE,
                              tol = tol,
                              tol.final = tol.final,
                              phase.number.limit = phase.number.limit,
                              sub.map.size.diff.limit =  sub.map.size.diff.limit,
                              info.tail = TRUE,
                              reestimate.single.ph.configuration = FALSE)
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
   fname<-paste0("simu_1_m-", m,
                 "_len-", map.length,
                 "_n.ind-", n.ind,
                 "_n.mrk-", n.mrk,
                 "_n.start-", n.sim.start,
                 "_n.stop-",n.sim.stop,
                 "_seed-",seed.for.config,
                 "_ttw-", thres.twopt,
                 "_etail-", extend.tail,
                 "_thmm-", thres.hmm,
                 "_t-",tol,
                 "_pp_", prefPairing,
                 "_qua_", quadrivalents,
                 "_tf-", tol.final, ".RData")
   save(list=c("ph.temp", "results", "dat"), file=fname)
 }
 