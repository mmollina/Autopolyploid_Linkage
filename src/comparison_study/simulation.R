simulate<-function(m = 6, map.length = 100, n.ind = 200, n.mrk = 500, n.chr  =  1,
                   prefPairing  =  0, quadrivalents  =  0, prob.dose = NULL,
                   centromere  =  20, natural.pairing  =  0,  seed.for.config  =  NULL,
                   seed.for.pop  =  NULL)
{
  #####
  ##Simulating homologous chromosomes
  #####
  w<-NULL
  for(i in 1:n.chr)
  {
    ph.temp<-sim_homologous(m=m,
                            n.mrk=n.mrk,
                            max.d=m/2,    #max dosage number
                            max.ph=max.ph,
                            restriction = FALSE,
                            prob.dose = prob.dose,
                            seed = seed.for.config)
    php<-ph_list_to_matrix(ph.temp$hom.allele.p, m)
    phq<-ph_list_to_matrix(ph.temp$hom.allele.q, m)
    w<-rbind(w, cbind(php,phq))
  }
  #####
  ## Simnulation population
  #####
  i<-1
  {
    ##Simulating full-sib population
    file.prefix<-paste0("sim_m_", m,
                        "seed_ph", seed.for.config,
                        "seed_pop", seed.for.pop,
                        "_pp_", prefPairing,
                        "_qua_", quadrivalents,
                        "_nch_", n.chr,
                        paste0(prob.dose, collapse = "-"),
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
    write(x = paste("SEED =", seed.for.pop), file = paste0(file.prefix, ".par"), append = TRUE)
    
    #chromosome file
    x <- a <- NULL
    for(j in 1:n.chr){
      a <- rbind(a, c(paste0("Ch",j), map.length,	centromere,		prefPairing,		quadrivalents)) 
      x <- rbind(x, data.frame(paste0("M",(1:n.mrk) + n.mrk * (j-1)), paste0("Ch",j), seq(0, map.length,length.out = n.mrk)))
    }
    
    write.table(a, file =paste0(file.prefix, ".chrom"), quote = FALSE, row.names = FALSE, col.names = c("chromosome",	"length",	"centromere",	"prefPairing", "quadrivalents"))
    write.table(x, file =paste0(file.prefix, ".map"), quote = FALSE, row.names = FALSE, col.names = c("marker", "chromosome", "position"))
    
    ## parent genotypes
    y<-data.frame(as.character(x[,1]), w)    
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
    
    dat$n.ind<-n.ind
    dat$n.mrk<-n.mrk * n.chr
    dat$m<-m
    dat$ind.names<-paste0("Ind_", 1:n.ind)
    dat$mrk.names<-as.character(x[,1])
    dat$dosage.p<-apply(w[,1:m], 1, function(x) sum(x==1))
    dat$dosage.q<-apply(w[,(m+1):(2*m)], 1, function(x) sum(x==1))
    dat$sequence<-rep(1:n.chr, each = n.mrk)
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
    dat.mappoly<-dat
    dat.polymapr<-cbind(dat$dosage.p, dat$dosage.q, dat$geno.dose)
    colnames(dat.polymapr)[1:2]<-c("P1", "P2")
    dat.polymapr <- as.matrix(dat.polymapr)
    # move PedrigreeSim files to adequate dir 
    system(paste0("mv ", file.prefix, "* pedsim_files/"))
  }
  #######################################
  ## Probability distribution of the genotypes
  #######################################
  myfunc<-function(x, m)
  {
    if(x[1] > m)
      return(mappoly::segreg_poly(m = m, dP = x[2], dQ = x[3]))
    else
    {
      y<-rep(0, m+1)
      y[x[1]+1]<-1
      return(y)
    }
  }
  x<-as.data.frame(as.table(as.matrix(dat.mappoly$geno.dose)))
  colnames(x)<-c("mrk", "ind", "dose")
  x$dose.p<-rep(dat.mappoly$dosage.p, dat.mappoly$n.ind)
  x$dose.q<-rep(dat.mappoly$dosage.q, dat.mappoly$n.ind)
  y<-t(apply(x[,-c(1:2)], 1, myfunc, m = dat.mappoly$m))
  colnames(y)<-0:dat.mappoly$m
  z<-cbind(x[,1:2], y)
  dat.mappoly$geno<-z
  return(list(dat.mp=dat.mappoly, dat.pm=dat.polymapr, ph.temp= ph.temp))
}
