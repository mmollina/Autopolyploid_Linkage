## Auxiliary codes for the Mollinari and Garcia 2018
## ------------------------------------------------------
## EXAMPLE in Supporting Information 3
##
## System information
###############################################
## Sys.info()
##                                      sysname 
##                                      "Linux" 
##                                      release 
##                          "4.13.0-38-generic" 
##                                      version 
##"#43-Ubuntu SMP Wed Mar 14 15:20:44 UTC 2018" 
##                                     nodename 
##                                    "haldane" 
##                                      machine 
##                                     "x86_64" 

devtools::install_github("mmollina/mappoly")
require(mappoly)

##NOTE: The data could be simulated using the following
##commented code. However, due to operational systems'
##particularities, we are using the same dataset presented
##in the S3 Appendix. You can try to use other seeds and
##perform different simulations.
##Simulation Arguments
## m<-4           # ploidy level
## n.ind<-30      # number of individuals
## n.mrk<-3       # number of markers
## Linkage phase configuration
##
## -o--------      -o---------
## -o------o-      --------o--
## --------o-  X   --------o--
## ----o-----      ----o---o--
##
## h.temp<-list(hom.allele.p=list(c(1,2), 4, c(2,3)),
##              hom.allele.q=list(1, 4, c(2,3,4)),
##              p=c(2,1,2),
##              q=c(1,1,3))
##
## Simulating full-sib population
## dat<-poly_cross_simulate(m=m,
##                          rf.vec=mf_h(1),
##                          n.mrk=n.mrk,
##                          n.ind=n.ind,
##                          hom.allele = h.temp,
##                          seed=8764,
##                          draw=TRUE,
##                          file="test.pdf")

## Creating the dataset used in Supporting Information 3
dat<-list()
dat$m <- 4
dat$n.ind <- 30
dat$n.mrk <- 3
dat$ind.names <- paste0("Ind_", 1:30)
dat$mrk.names <- paste0("M_", 1:3)  
dat$dosage.p <- c(2,1,2)
dat$dosage.q <- c(1,1,3)
dat$sequence <- NA
dat$sequence.pos <- NA
dat$geno.dose <- matrix(c(2,1,3,1,2,1,2,2,2,3,1,2,1,2,1,1,3,2,2,2,2,3,2,2,1,3,2,2,1,1,
                          0,1,0,0,1,2,1,1,0,1,1,0,1,0,2,2,1,0,1,1,0,0,0,1,1,0,1,1,2,2,
                          3,2,2,3,3,3,3,3,3,2,2,2,3,2,2,2,2,3,3,3,3,2,2,1,4,2,2,2,3,2),
                          byrow = TRUE, nrow = 3, ncol = 30,
                          dimnames=list(dat$mrk.names, dat$ind.names)) 
dat$nphen <- 0
dat$phen  <- NULL
class(dat)<- 'mappoly.data'

#Making a sequence with all markers
all_mrk<-make_seq_mappoly(dat, arg="all")
counts<-cache_counts_twopt(input.seq=all_mrk)

##Calculating the pairwise recombination fraction for all markers
all_pairs<-est_pairwise_rf(input.seq=all_mrk,
                           counts)

round(all_pairs$pairwise$`1-2`, 2)
round(all_pairs$pairwise$`1-3`, 2)
round(all_pairs$pairwise$`2-3`, 2)

#Check the most likely linkage phase configurations using the two
#point analysis given a LOD Score threshold
ph <- ls_linkage_phases(input.seq = all_mrk, thres = 1, twopt = all_pairs)
plot(ph)

#Reconstruction a map with a known order and the subset of linkage
#phases previously computed
res<-est_rf_hmm(input.seq=all_mrk, input.ph=ph, verbose=TRUE, tol=10e-7)
get_LOD(res)
plot(res, config = 1)
