load("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
require(mappoly)
for(i in 1:12){
  denovo.dist<-est_full_hmm_with_prior_dist(input.map = all.maps[[i]][[1]], dat.dist = solcap.dat, 
                                            phase.config = 1, tol = 10e-5, verbose = FALSE)
  genomic.dist<-est_full_hmm_with_prior_dist(input.map = all.maps[[i]][[2]], dat.dist = solcap.dat, 
                                             phase.config = 1, tol = 10e-5, verbose = FALSE)
  denovo.error<-est_full_hmm_with_global_error(input.map = all.maps[[i]][[1]], 
                                               error = 0.05, tol = 10e-5)
  genomic.error<-est_full_hmm_with_global_error(input.map = all.maps[[i]][[2]], 
                                                error = 0.05, tol = 10e-5)
  all.maps[[i]][[3]]<-denovo.dist$reestimated.map
  all.maps[[i]][[4]]<-genomic.dist$reestimated.map
  all.maps[[i]][[5]]<-denovo.error
  all.maps[[i]][[6]]<-genomic.error
  save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
}

