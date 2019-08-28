# Load MAPpoly
require(mappoly)
# Get length and loglike of the map
map_summary<-function(map){
  len <- sum(imf_h(map$maps[[1]]$seq.rf))
  loglike <- map$maps[[1]]$loglike
  return(c(len, loglike))
}
# Convert map from MAPpoly into dataframes
extract_maps<-function(x, dat, type = NA, lg = NA){
    d <- solcap.dat$sequence.pos[x$maps[[1]]$seq.num]/1e6
    if(cor(d, seq_along(d)) < 0) d<-rev(d)
    w<-data.frame(mrk = dat$mrk.names[x$maps[[1]]$seq.num],
               map.pos = round(cumsum(c(0,imf_h(x$maps[[1]]$seq.rf))),1),  
               geno.pos = d,
               ch = solcap.dat$sequence[x$maps[[1]]$seq.num],
               type = type, lg = lg)
  return(w)
}
MAPs.df<-NULL
final.maps<-NULL
all.loglike<-NULL
all.map.length<-NULL
##Gather maps from 12 chromosomes
load("~/repos/Autopolyploid_Linkage/src/solcap_map_construction/solcap_map_results.RData")
MAPs.temp<-vector("list", 12)
for(i in 1:12){
  cat("ch: ", i, "\n")
  MAPs.df<-rbind(MAPs.df, rbind(extract_maps(MAPs[[i]][[1]], dat = solcap.dat, type = "denovo", lg = i),
                                        extract_maps(MAPs[[i]][[2]], dat = solcap.dat, type = "genomic", lg = i),
                                        extract_maps(MAPs[[i]][[3]], dat = solcap.dat, type = "denovo.dist", lg = i),
                                        extract_maps(MAPs[[i]][[4]], dat = solcap.dat, type = "genomic.dist", lg = i),
                                        extract_maps(MAPs[[i]][[5]], dat = solcap.dat, type = "denovo.error", lg = i),
                                        extract_maps(MAPs[[i]][[6]], dat = solcap.dat, type = "genomic.error", lg = i)))
  
  all.map.length<-rbind(all.map.length, data.frame(denovo = sum(imf_h(MAPs[[i]][[1]]$maps[[1]]$seq.rf)),
                                                   genomic = sum(imf_h(MAPs[[i]][[2]]$maps[[1]]$seq.rf)),
                                                   denovo.dist = sum(imf_h(MAPs[[i]][[3]]$maps[[1]]$seq.rf)),
                                                   genomic.dist = sum(imf_h(MAPs[[i]][[4]]$maps[[1]]$seq.rf)),
                                                   denovo.error = sum(imf_h(MAPs[[i]][[5]]$maps[[1]]$seq.rf)),
                                                   genomic.error = sum(imf_h(MAPs[[i]][[6]]$maps[[1]]$seq.rf))))
  
  all.loglike<-rbind(all.loglike, data.frame(denovo = MAPs[[i]][[1]]$maps[[1]]$loglike,
                                             genomic = MAPs[[i]][[2]]$maps[[1]]$loglike,
                                             denovo.dist = MAPs[[i]][[3]]$maps[[1]]$loglike,
                                             genomic.dist = MAPs[[i]][[4]]$maps[[1]]$loglike,
                                             denovo.error = MAPs[[i]][[5]]$maps[[1]]$loglike,
                                             genomic.error = MAPs[[i]][[6]]$maps[[1]]$loglike))
  
  final.maps<-rbind(final.maps, data.frame(extract_maps(MAPs[[i]][[1]], dat = solcap.dat, type = "denovo", lg = i),
                                           cbind(ph_list_to_matrix(L = MAPs[[i]][[1]]$maps[[1]]$seq.ph$P, m = 4),
                                                 ph_list_to_matrix(L = MAPs[[i]][[1]]$maps[[1]]$seq.ph$Q, m = 4))))
}
head(final.maps)
final.maps <- final.maps[, -c(4,5)]
colnames(final.maps)[5:12] <- paste0("H", 1:8)  

#Check if denovo and genomic maps contain the same SNPs in order to make the likelihood comparable
for(i in 1:12)
  print(identical(sort(MAPs[[1]][[1]]$maps[[1]]$seq.num), sort(MAPs[[1]][[2]]$maps[[1]]$seq.num)))

# Loglikes in base 10
all.loglike<-round(all.loglike, 1)
apply(all.loglike[,1:2], 1, diff)
apply(all.loglike[,3:4], 1, diff)
apply(all.loglike[,5:6], 1, diff)

# Map lengths
all.map.length<-round(all.map.length,1)
nmrk<-sapply(MAPs, function(x) sapply(x, function(x) x$info$n.mrk))[1,]

# Table S8.1 from Mollinari and Garcia (2019) - Map lengths and LOD Scores
tab<-cbind(all.map.length[,c(1,2)],
           apply(all.loglike[,1:2], 1, diff),
           all.map.length[,c(3,4)],
           apply(all.loglike[,3:4], 1, diff),
           all.map.length[,c(5,6)],
           apply(all.loglike[,5:6], 1, diff), 
           nmrk)
colnames(tab)<-c("len", "len", "LOD", "len", "len", "LOD", "len", "len", "LOD", "nmrk")
xtable::xtable(tab, digits = 1)

##Phased map
head(final.maps)
nrow(final.maps)
write.csv(x = final.maps, file =  "~/repos/Autopolyploid_Linkage/src/solcap_map_construction/phased_SolCAP_map.csv")

##Percentage of Single and Multidose markers
multiplex<-sum(apply(final.maps[,c("H1", "H2", "H3", "H4")], 1, sum)==2  |
                 apply(final.maps[,c("H5", "H6", "H7", "H8")], 1, sum)==2)
simplex<-nrow(final.maps)-multiplex
simplex;multiplex
round(100*simplex/nrow(final.maps))
round(100*multiplex/nrow(final.maps))

##Figure S8.1 from Mollinari and Garcia (2019) - Map versus genome
head(MAPs.df)
require(ggplot2)
p <- ggplot(MAPs.df, aes(x = geno.pos, y=map.pos)) + geom_point()
p + aes(color=type) + facet_wrap(~ lg) +  
  scale_colour_brewer("Mapping procedure", palette="Paired") + 
  xlab("Genome position (Mb)") + ylab("Map position (cM)")

##Phase is the same for MDS and genomic
ph.compare<-NULL
for(i in 1:12){
  p1<-MAPs[[i]][[1]]$maps[[1]]$seq.ph$P
  q1<-MAPs[[i]][[1]]$maps[[1]]$seq.ph$Q
  p2<-MAPs[[i]][[2]]$maps[[1]]$seq.ph$P
  q2<-MAPs[[i]][[2]]$maps[[1]]$seq.ph$Q    
  a<-intersect(names(p1), names(p2))
  ph.compare<-rbind(ph.compare, data.frame(P = compare_haplotypes(m = 4, h1 = p1[a], h2 = p2[a])$is.same.haplo, 
                                           Q = compare_haplotypes(m = 4, h1 = q1[a], h2 = q2[a])$is.same.haplo))
}
ph.compare
#end of file