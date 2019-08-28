#! /usr/bin/env Rscript

# Analysis of real tetraploid potato SNP data set:
# Genetic map of the B2721 population which is a 
# cross between tetraploid potato varieties Atlantic 
# and B1829-5.

#################
# Loading MAPpoly
#################
require(mappoly)
setwd("~/repos/Autopolyploid_Linkage/src/solcap_map_construction/")

####################
# Loading Data and Two-point analysis
####################
load(file = "SolCAP.rda")
all.mrk<-make_seq_mappoly(solcap.dat, "all")
print(solcap.dat, detailed = TRUE)
counts<-get_cache_two_pts_from_web(m = solcap.dat$m)
all.pairs<-est_pairwise_rf(input.seq = all.mrk,
                           count.cache = counts,
                           n.clusters = 16,
                           verbose=TRUE)
mat.full<-rf_list_to_matrix(input.twopt=all.pairs)
plot(mat.full)

###########
## Grouping
###########
lgs<-group_mappoly(input.mat = mat.full,
                   input.seq = all.mrk,
                   expected.groups = 12,
                   comp.mat = TRUE,
                   inter = FALSE)
print(lgs, detailed = TRUE)

###########
## Preparing MDS order and two-point info for sequential phasing
###########
id.ch <- order(c(11, 9, 7, 5, 4, 1, 8, 2, 10, 6, 12, 3)) 
n<-list(c(84, 325, 109, 175, 208, 114, 268, 326), #1
        c(2, 3, 1, 10, 21, 24, 27, 41, 74, 61 , 59, 107, 108, 104, 102), #2
        c(4, 221, 124, 257, 122, 289, 7, 131, 206, 253, 251, 80), #3
        c(43, 231, 230, 325), #4
        c(153, 197, 255, 50, 47), #5
        c(183, 184, 185, 237, 274, 40, 262, 391, 390, 340, 392, 380, 352, 324, 304, 254), #6
        c(177, 178, 152), #7
        c(21, 202, 218, 241, 14), #8
        c(25, 23, 26,  24, 31, 80, 169, 173, 179, 168, 166, 171, 167, 165, 172, 248), #9
        c(8, 45, 46, 47, 48, 49, 50, 97, 132, 189, 205), #10
        c(219, 215, 113, 74, 95, 94, 172, 43, 98, 101, 70, 63), #11
        c(1, 58, 41, 70, 157, 132)) #12
W<-vector("list", 12)
for(ch in 1:12){
  cat("\n~~~~~~~~~~~Ch", ch, "~~~~~~~~~~~\n")
  lg <- make_seq_mappoly(lgs, id.ch[ch])
  lg <- make_seq_mappoly(solcap.dat, lg$seq.num[which(lg$sequence==ch)])
  p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
  m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6,
                         thresh.LOD.rf = 6, thresh.rf = 0.5)
  lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6,
                           thresh.LOD.rf = 6, thresh.rf = 0.5,
                           thresh.perc = 0.05)
  p.filt <- make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
  m.filt <- rf_list_to_matrix(p.filt)
  mds.ord <- mds_mappoly(input.mat = m.filt, n = n[[ch]])
  cat("\n")
  mds.seq <- make_seq_mappoly(mds.ord)
  o<-mds.seq$seq.num
  if(cor(1:length(o), rev(o)) > cor(1:length(o), o))
    o<-rev(o)
  map.seq<-make_seq_mappoly(solcap.dat, o)
  twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
  W[[ch]] <- list(map.seq = map.seq,
                  twopt = twopt,
                  chr = paste0("Ch_", ch))
  pdf(file = paste0("Ch_", ch, "_mds_rec_frac_mat.pdf"), width = 20, height = 20)
  plot(m.filt, ord = W[[ch]]$map.seq$seq.mrk.names)
  dev.off()
}

###########
# Parallel sequential phasing
###########
map_contrsuction <- function(X){
  map.denovo <- est_rf_hmm_sequential(input.seq = X[[1]],
                                      start.set = 10,
                                      thres.twopt = 10, 
                                      thres.hmm = 10,
                                      extend.tail = 100,
                                      twopt = X[[2]],
                                      tol = 10e-3,
                                      tol.final = 10e-4,
                                      sub.map.size.diff.limit = 10)
  pdf(file = paste0(X[[3]], "_denovo.pdf"))
  plot(map.denovo)
  dev.off()
  v <- make_seq_mappoly(solcap.dat, map.denovo$maps[[1]]$seq.num)
  geno.ord <- get_genomic_order(v)
  map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = rownames(geno.ord))
  map.genomic <- est_rf_hmm_sequential(input.seq = map.seq.ord,
                                       start.set = 10,
                                       thres.twopt = 10,
                                       thres.hmm = 10, 
                                       extend.tail = 100,
                                       twopt = X[[2]],
                                       tol = 10e-3,
                                       tol.final = 10e-4,sub.map.size.diff.limit = 10)
  pdf(file = paste0(X[[3]], "_genomic.pdf"))
  plot(map.genomic)
  dev.off()
  denovo.dist <- est_full_hmm_with_prior_dist(input.map = map.denovo, dat.dist = solcap.dat, phase.config = 1, tol = 10e-4, verbose = FALSE)
  genomic.dist <- est_full_hmm_with_prior_dist(input.map = map.genomic, dat.dist = solcap.dat, phase.config = 1, tol = 10e-4, verbose = FALSE)
  denovo.error <- est_full_hmm_with_global_error(input.map = map.denovo, error = 0.05, tol = 10e-4)
  genomic.error <- est_full_hmm_with_global_error(input.map = map.genomic, error = 0.05, tol = 10e-4)
  return(list(map.denovo,
              map.genomic,
              denovo.dist,
              genomic.dist,
              denovo.error,
              genomic.error))
}

cl <- parallel::makeCluster(12)
parallel::clusterEvalQ(cl, require(mappoly))
parallel::clusterExport(cl, "solcap.dat")
system.time(MAPs <- parallel::parLapply(cl, W, map_contrsuction))
parallel::stopCluster(cl)

###########
# Saving results
###########
save.image("solcap_map_results.RData")
