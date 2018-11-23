require(mappoly)
setwd("~/repos/Autopolyploid_Linkage/src/solcap_map_constructiomn/")
#####
solcap.dat<-read_geno_dist(file.in = "~/repos/Autopolyploid_Linkage/src/solcap_map_constructiomn/snp_calling/SolCAP", prob.thres = 0.95)
solcap.dat$sequence[solcap.dat$sequence==0]<-NA
all.mrk<-make_seq_mappoly(solcap.dat, "all")
print(solcap.dat, detailed = TRUE)
plot(solcap.dat)
counts<-get_cache_two_pts_from_web(m = solcap.dat$m)
# all.pairs<-est_pairwise_rf(input.seq = all.mrk,
#                           count.cache = counts,
#                           n.clusters = 16,
#                           verbose=TRUE)
#
# save(all.pairs, file = "~/repos/tutorials/solcap/all_pairs_fittetra.rda", compress = TRUE)
load("~/repos/tutorials/solcap/all_pairs_fittetra.rda")
mat.full<-rf_list_to_matrix(input.twopt=all.pairs)
#plot(mat.full)
###########
## Grouping
###########
lgs<-group_mappoly(input.mat = mat.full,
                   input.seq = all.mrk,
                   expected.groups = 12,
                   comp.mat = TRUE,
                   inter = FALSE)
print(lgs, detailed = TRUE)

all.maps<-vector("list", 12)
##########
# Group 1
##########
ch <- 1
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6,
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6,
                         thresh.LOD.rf = 6, thresh.rf = 0.5,
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n = c(94, 335, 6, 119, 185, 218, 124, 336))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
map.seq<-make_seq_mappoly(solcap.dat, mds.seq$seq.num[-which(is.na(mds.seq$sequence))])
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 10,
                                  thres.hmm = 10,
                                  extend.tail = 100,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE,
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 10,
                                   extend.tail = 100,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE,
                                   high.prec = FALSE)
all.maps[[ch]]<-list(map.denovo,
                    map.genomic)
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
##########
# Group 2
##########
ch <- 2
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6,
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6,
                         thresh.LOD.rf = 6, thresh.rf = 0.5,
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n= c(100, 190, 187, 113, 223,
                                              180, 188, 201,  12,
                                              16, 229, 181, 11, 107,
                                              26, 39, 17, 233, 192,
                                              175, 31))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
map.seq<-make_seq_mappoly(solcap.dat, mds.seq$seq.num[-which(is.na(mds.seq$sequence))])
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  start.set = 10,
                                  thres.twopt = 10,
                                  thres.hmm = 20,
                                  extend.tail = 100,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE,
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10,
                                   start.set = 50,
                                   thres.hmm = 25,
                                   extend.tail = 100,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE,
                                   high.prec = FALSE)
all.maps[[ch]]<-list(map.denovo,
                     map.genomic)
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
###########
## Group 3
###########
ch <- 3
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6, 
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6, 
                         thresh.LOD.rf = 6, thresh.rf = 0.5, 
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n = c(4, 221,124,257, 122, 289, 7, 131, 206, 2100, 253, 251, 80))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
map.seq<-make_seq_mappoly(solcap.dat, mds.seq$seq.num[-which(is.na(mds.seq$sequence))])
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 10,
                                  thres.hmm = 10,
                                  extend.tail = 100,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE, 
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 10,
                                   extend.tail = 100,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE, 
                                   high.prec = FALSE)
all.maps[[ch]]<-list(map.denovo, 
                     map.genomic)
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
###########
## Group 4
###########
ch <- 4
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6, 
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6, 
                         thresh.LOD.rf = 6, thresh.rf = 0.5, 
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n = c(54,242, 241, 336))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
map.seq<-make_seq_mappoly(solcap.dat, mds.seq$seq.num[-which(is.na(mds.seq$sequence))])
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 10,
                                  thres.hmm = 10,
                                  extend.tail = 100,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE, 
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 10,
                                   extend.tail = 100,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE, 
                                   high.prec = FALSE)
all.maps[[ch]]<-list(map.denovo, 
                     map.genomic)
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
###########
## Group 5
###########
ch <- 5
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6, 
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6, 
                         thresh.LOD.rf = 6, thresh.rf = 0.5, 
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n = c(159,203,261, 56, 53))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
o<-mds.seq$seq.num[-which(is.na(mds.seq$sequence))]
if(cor(1:length(o), rev(o)) > cor(1:length(o), o))
  o<-rev(o)
map.seq<-make_seq_mappoly(solcap.dat, o)
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 10,
                                  thres.hmm = 10,
                                  extend.tail = 50,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE, 
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 10,
                                   extend.tail = 50,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE, 
                                   high.prec = FALSE)

all.maps[[ch]]<-list(map.denovo, 
                     map.genomic)
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
###########
## Group 6
###########
ch <- 6
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6, 
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6, 
                         thresh.LOD.rf = 6, thresh.rf = 0.5, 
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n = c(241, 278, 44, 266, 395, 394, 344, 396, 384, 356, 328, 308, 258))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
o<-mds.seq$seq.num[-which(is.na(mds.seq$sequence))]
if(cor(1:length(o), rev(o)) > cor(1:length(o), o))
  o<-rev(o)
map.seq<-make_seq_mappoly(solcap.dat, o)
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 25,
                                  thres.hmm = 10,
                                  extend.tail = 100,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE, 
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v, decreasing = TRUE))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 20,
                                   extend.tail = 50,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE, 
                                   high.prec = FALSE)

all.maps[[ch]]<-list(map.denovo, 
                     rev_map(map.genomic))
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
###########
## Group 7
###########
ch <- 7
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6, 
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6, 
                         thresh.LOD.rf = 6, thresh.rf = 0.5, 
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n = c(183, 184, 158))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
o<-mds.seq$seq.num[-which(is.na(mds.seq$sequence))]
if(cor(1:length(o), rev(o)) > cor(1:length(o), o))
o<-rev(o)
map.seq<-make_seq_mappoly(solcap.dat, o)
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 10,
                                  thres.hmm = 10,
                                  extend.tail = 100,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE, 
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 10,
                                   extend.tail = 100,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE, 
                                   high.prec = FALSE)

all.maps[[ch]]<-list(map.denovo, 
                     map.genomic)
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
###########
## Group 8
###########
ch <- 8
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6, 
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6, 
                         thresh.LOD.rf = 6, thresh.rf = 0.5, 
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n = c(28, 209, 225, 248, 21))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
o<-mds.seq$seq.num[-which(is.na(mds.seq$sequence))]
if(cor(1:length(o), rev(o)) > cor(1:length(o), o))
  o<-rev(o)
map.seq<-make_seq_mappoly(solcap.dat, o)
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 10,
                                  thres.hmm = 10,
                                  extend.tail = 100,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE, 
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 10,
                                   extend.tail = 100,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE, 
                                   high.prec = FALSE)

all.maps[[ch]]<-list(map.denovo, 
                     map.genomic)
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
###########
## Group 9
###########
ch <- 9
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6, 
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6, 
                         thresh.LOD.rf = 6, thresh.rf = 0.5, 
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n=c(104, 105, 94, 270, 214, 208, 11, 248, 181, 268, 135,  266))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
o<-mds.seq$seq.num[-which(is.na(mds.seq$sequence))]
if(cor(1:length(o), rev(o)) > cor(1:length(o), o))
  o<-rev(o)
map.seq<-make_seq_mappoly(solcap.dat, o)
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 10,
                                  thres.hmm = 10,
                                  extend.tail = 200,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE, 
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 10,
                                   extend.tail = 200,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE, 
                                   high.prec = FALSE)

all.maps[[ch]]<-list(map.denovo, 
                     map.genomic)
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
###########
## Group 10
###########
ch <- 10
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6, 
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6, 
                         thresh.LOD.rf = 6, thresh.rf = 0.5, 
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n= c(100, 192, 208, 53))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
o<-mds.seq$seq.num[-which(is.na(mds.seq$sequence))]
if(cor(1:length(o), rev(o)) > cor(1:length(o), o))
  o<-rev(o)
map.seq<-make_seq_mappoly(solcap.dat, o)
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 10,
                                  thres.hmm = 10,
                                  extend.tail = 100,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE, 
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 10,
                                   extend.tail = 100,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE, 
                                   high.prec = FALSE)

all.maps[[ch]]<-list(map.denovo, 
                     map.genomic)
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
###########
## Group 11
###########
ch <- 11
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6, 
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6, 
                         thresh.LOD.rf = 6, thresh.rf = 0.5, 
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n = c(222,218, 115, 76, 97, 174, 45))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
o<-mds.seq$seq.num[-which(is.na(mds.seq$sequence))]
if(cor(1:length(o), rev(o)) > cor(1:length(o), o))
  o<-rev(o)
map.seq<-make_seq_mappoly(solcap.dat, o)
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 10,
                                  thres.hmm = 10,
                                  extend.tail = 200,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE, 
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 10,
                                   extend.tail = 200,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE, 
                                   high.prec = FALSE)

all.maps[[ch]]<-list(map.denovo, 
                     map.genomic)
save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
###########
## Group 12
###########
ch <- 12
lg <- make_seq_mappoly(lgs, ch)
lg <- make_seq_mappoly(solcap.dat, lg$seq.num[-which(lg$sequence!=ch)])
p <- make_pairs_mappoly(input.twopt = all.pairs, input.seq = lg)
m <- rf_list_to_matrix(input.twopt = p, thresh.LOD.ph = 6, 
                       thresh.LOD.rf = 6, thresh.rf = 0.5)
plot(m, index = FALSE)
lg.filt <- rf_snp_filter(input.twopt = p, thresh.LOD.ph = 6, 
                         thresh.LOD.rf = 6, thresh.rf = 0.5, 
                         thresh.perc = 0.05)
plot(m, index = FALSE, ord = lg.filt$seq.mrk.names)
p.filt<-make_pairs_mappoly(input.seq = lg.filt, input.twopt = all.pairs)
m.filt<-rf_list_to_matrix(p.filt)
plot(m.filt, index = FALSE)
mds.ord<-mds_mappoly(input.mat = m.filt, n = c(60, 43, 72, 159, 134))
plot(mds.ord)
mds.seq<-make_seq_mappoly(mds.ord)
o<-mds.seq$seq.num[-which(is.na(mds.seq$sequence))]
if(cor(1:length(o), rev(o)) > cor(1:length(o), o))
  o<-rev(o)
map.seq<-make_seq_mappoly(solcap.dat, o)
plot(m, ord = map.seq$seq.mrk.names)
twopt <- make_pairs_mappoly(input.seq = map.seq, input.twopt = all.pairs)
map.denovo<-est_rf_hmm_sequential(input.seq = map.seq,
                                  thres.twopt = 10, start.set = 10,
                                  thres.hmm = 10,
                                  extend.tail = 100,
                                  twopt = twopt,
                                  verbose = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4,
                                  phase.number.limit = 40,
                                  sub.map.size.diff.limit = 10,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE, 
                                  high.prec = FALSE)
v<-map.denovo$maps[[1]]$seq.num
map.seq.ord <- make_seq_mappoly(input.obj = solcap.dat, arg = sort(v))
map.genomic<-est_rf_hmm_sequential(input.seq = map.seq.ord,
                                   thres.twopt = 10, start.set = 10,
                                   thres.hmm = 10,
                                   extend.tail = 100,
                                   twopt = twopt,
                                   verbose = TRUE,
                                   tol = 10e-3,
                                   tol.final = 10e-4,
                                   phase.number.limit = Inf,
                                   sub.map.size.diff.limit = Inf,
                                   info.tail = TRUE,
                                   reestimate.single.ph.configuration = FALSE, 
                                   high.prec = FALSE)

all.maps[[ch]]<-list(map.denovo, 
                     map.genomic)

save.image("~/repos/tutorials/solcap/fittreta_and_mappoly.rda")
pdf(file = paste0("ch",ch,".pdf"))
plot(map.denovo)
plot(map.genomic)
dev.off()
try(dev.off())