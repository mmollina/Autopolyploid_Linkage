#####
## Map Construction using MAPpoly
#####
build_map_MAPpoly<-function(dat,
                            start.set = 4,
                            thres.twopt = 10,
                            thres.hmm = 10,
                            extend.tail = 100,
                            verbose = TRUE,
                            tol = 10e-2,
                            tol.final = 10e-3,
                            phase.number.limit = 20,
                            sub.map.size.diff.limit = Inf){
  counts<-get_cache_two_pts_from_web(m = dat$m)
  ##Making a sequence with all markers
  all.dat<-make_seq_mappoly(dat, arg="all")
  ##Calculating the pairwise recombination fraction for all markers
  w0<-system.time(all.pairs<-est_pairwise_rf(input.seq = all.dat,
                             count.cache = counts,
                             n.clusters = 1,
                             verbose = TRUE))
  mat<-rf_list_to_matrix(all.pairs)
  mds.res<-mds_mappoly(input.mat = mat)
  s<-make_seq_mappoly(mds.res)
  w1<-system.time(res.orig<-est_rf_hmm_sequential(input.seq = all.dat,
                                                  twopt = all.pairs,
                                                  start.set = start.set,
                                                  thres.twopt = thres.twopt,
                                                  thres.hmm = thres.hmm,
                                                  extend.tail = extend.tail,
                                                  verbose = verbose,
                                                  tol = tol,
                                                  tol.final = tol.final,
                                                  phase.number.limit = phase.number.limit,
                                                  sub.map.size.diff.limit = sub.map.size.diff.limit,
                                                  info.tail = TRUE,
                                                  reestimate.single.ph.configuration = FALSE,
                                                  high.prec = FALSE))
  
  w2<-system.time(res.mds<-est_rf_hmm_sequential(input.seq = s,
                                                 twopt = all.pairs,
                                                 start.set = start.set,
                                                 thres.twopt = thres.twopt,
                                                 thres.hmm = thres.hmm,
                                                 extend.tail = extend.tail,
                                                 verbose = verbose,
                                                 tol = tol,
                                                 tol.final = tol.final,
                                                 phase.number.limit = phase.number.limit,
                                                 sub.map.size.diff.limit = sub.map.size.diff.limit,
                                                 info.tail = TRUE,
                                                 reestimate.single.ph.configuration = FALSE,
                                                 high.prec = FALSE))
  
  return(list(res.orig = res.orig, res.mds = res.mds, time.twopt = w0, time.mds = w1, time.orig = w2))
}