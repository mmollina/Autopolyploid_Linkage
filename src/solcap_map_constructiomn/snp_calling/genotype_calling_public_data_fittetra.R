require(fitTetra)
require(mappoly)
setwd("~/repos/tutorials/solcap/snp_calling/fitTetra/")
datxy<-read.csv(file = "solcap.csv", row.names = 1)
datxy[1:10,1:10]
colnames(datxy)
X<-grep(pattern = ".X.Raw", colnames(datxy))
Y<-grep(pattern = ".Y.Raw", colnames(datxy))
Theta<-grep(pattern = ".Theta", colnames(datxy))
R<-grep(pattern = ".R", colnames(datxy))
dat.X.Raw<-datxy[,X]
dat.Y.Raw<-datxy[,Y]
dat.Theta<-datxy[,Theta]
dat.R<-datxy[,seq(4,ncol(datxy),4)]
ind.names<-sapply(strsplit(colnames(dat.X.Raw), split = "X4x_|\\.X.Raw"), function(x) x[2])[-c(2,4)]
dat.X.Raw<-cbind(apply(dat.X.Raw[,1:2], 1, mean), apply(dat.X.Raw[,3:4], 1, mean), dat.X.Raw[,-c(1:4)])
dat.Y.Raw<-cbind(apply(dat.Y.Raw[,1:2], 1, mean), apply(dat.Y.Raw[,3:4], 1, mean), dat.Y.Raw[,-c(1:4)])
dat.Theta<-cbind(apply(dat.Theta[,1:2], 1, mean), apply(dat.Theta[,3:4], 1, mean), dat.Theta[,-c(1:4)])
dat.R<-cbind(apply(dat.R[,1:2], 1, mean), apply(dat.R[,3:4], 1, mean), dat.R[,-c(1:4)])
colnames(dat.R)<-colnames(dat.Theta)<-colnames(dat.Y.Raw)<-colnames(dat.X.Raw)<-ind.names
dat.Y.Raw[1:3,1:3]

dat<-reshape::melt(as.matrix(dat.X.Raw))
colnames(dat)<-c("MarkerName", "SampleName", "X_Raw")
a1<-reshape::melt(as.matrix(dat.Y.Raw))
dat$Y_Raw<-a1$value
a2<-reshape::melt(as.matrix(dat.Theta))
dat$Theta<-a2$value
a3<-reshape::melt(as.matrix(dat.R))
dat$R<-a3$value

###Running fitTetra
df.tetra <- with(dat, data.frame(MarkerName=MarkerName, 
                                 SampleName=SampleName, ratio=X_Raw/(X_Raw+Y_Raw)))

#################################################################################
## Genomic infomation
load("~/repos/tutorials/solcap/blast/genomic_position.RData")
head(blast)

head(df.tetra)
final.df<-NULL
cte2<-0
cte<-0
miss.perc<-25
dose.p<-dose.q<-NULL
for(k in as.character(blast$query_id))
{
  cte<-cte+1
  cat("ch: ", blast[k, c("chr")], "pos: ", blast[k, c("snp.pos")], "--->", cte, "  ")
  i <- match(k, unique(df.tetra$MarkerName))
  if(is.na(i)){
    cat("skip marker\n")
    next
  }
  fit <- try(fitTetra(marker=i, data=df.tetra, try.HW = FALSE, plot="fitted", p.threshold = 0.8), TRUE)
  if(class(fit) == "try-error"){
    cat("skip marker\n")
    next
  } 
  a<-fit$scores
  
  if(length(a)==1){
    cat("skip marker\n")
    next
  }
  
  dp1<-a["population1_P1",15]
  dp2<-a["population1_P2",15]
  
  if(any(is.na(c(dp1, dp2)))){
    cat("skip marker\n")
    next
  }
  cat(" ", k, " ")
  w<-segreg_poly(m = 4, dP = dp1, dQ = dp2)
  if(any(w==1)){
    cat("\n")
    next()
  } 
  names(w)<-0:4
  id <- grep("_P", rownames(a), invert = TRUE)
  v<-a[id,"geno"]
  u <- is.na(v) + v%in%names(w[w==0]) > 0
  y<-a[id,8:12]
  if(any(u)){
    v[u]<-NA
    for(l in which(u))
      y[l,]<-w
  }
  if(sum(is.na(v)) > miss.perc*length(v)/100){
    cat("\n")
    next()
  } 
  x.temp<-table(v)
  x<-rep(0, 5); names(x)<-0:4
  x[names(x.temp)]<-x.temp
  pv<-chisq.test(x[names(w[w!=0])], p = w[names(w[w!=0])])$p.value
  if(pv > 10e-4){
    cte2<-cte2+1
    dose.p<-c(dose.p, dp1)
    dose.q<-c(dose.q, dp2)
    cat(" ---> passed! mrks so far:", cte2, "\n")
    final.df<-rbind(final.df, data.frame(mrk = fit$modeldata$markername, 
                                         ind = rownames(y), round(y,4), row.names = NULL))
  } else cat("\n")
}

rest<-setdiff(unique(df.tetra$MarkerName), as.character(blast$query_id))

for(k in rest)
{
  cte<-cte+1
  cat("ch: ", blast[k, c("chr")], "pos: ", blast[k, c("snp.pos")], "--->", cte, "  ")
  i <- match(k, unique(df.tetra$MarkerName))
  if(is.na(i)){
    cat("skip marker\n")
    next
  }
  fit <- try(fitTetra(marker=i, data=df.tetra, try.HW = FALSE, plot="fitted", p.threshold = 0.8), TRUE)
  if(class(fit) == "try-error"){
    cat("skip marker\n")
    next
  } 
  a<-fit$scores
  
  if(length(a)==1){
    cat("skip marker\n")
    next
  }
  
  dp1<-a["population1_P1",15]
  dp2<-a["population1_P2",15]
  
  if(any(is.na(c(dp1, dp2)))){
    cat("skip marker\n")
    next
  }
  cat(" ", k, " ")
  w<-segreg_poly(m = 4, dP = dp1, dQ = dp2)
  if(any(w==1)){
    cat("\n")
    next()
  } 
  names(w)<-0:4
  id <- grep("_P", rownames(a), invert = TRUE)
  v<-a[id,"geno"]
  u <- is.na(v) + v%in%names(w[w==0]) > 0
  y<-a[id,8:12]
  if(any(u)){
    v[u]<-NA
    for(l in which(u))
      y[l,]<-w
  }
  if(sum(is.na(v)) > miss.perc*length(v)/100){
    cat("\n")
    next()
  } 
  x.temp<-table(v)
  x<-rep(0, 5); names(x)<-0:4
  x[names(x.temp)]<-x.temp
  pv<-chisq.test(x[names(w[w!=0])], p = w[names(w[w!=0])])$p.value
  if(pv > 10e-4){
    cte2<-cte2+1
    dose.p<-c(dose.p, dp1)
    dose.q<-c(dose.q, dp2)
    cat(" ---> passed! mrks so far:", cte2, "\n")
    final.df<-rbind(final.df, data.frame(mrk = fit$modeldata$markername, 
                                         ind = rownames(y), round(y,4), row.names = NULL))
  } else cat("\n")
}

## Ordering by individuals
DF<-final.df[order(final.df$ind),]
head(DF, 10)

#################################################################################
chr<-blast$chr[match(unique(DF$mrk), blast$query_id)]
seq.pos<-blast$snp.pos[match(unique(DF$mrk), blast$query_id)]   

## Writing mappoly input file
indnames<-as.character(unique(DF$ind))
mrknames<-as.character(unique(DF$mrk))
write(paste("ploidy", 4), file="SolCAP")
write(paste("nind", length(indnames)), file="SolCAP", append=TRUE)
write(paste("nmrk", length(mrknames)), file="SolCAP", append=TRUE)
cat("mrknames", mrknames, file="SolCAP", append=TRUE)
cat("\nindnames", indnames, file="SolCAP", append=TRUE)
cat("\ndosageP", dose.p, file="SolCAP", append=TRUE)
cat("\ndosageQ", dose.q, file="SolCAP", append=TRUE)
cat("\nseq", chr, file="SolCAP", append=TRUE)
cat("\nseqpos", seq.pos, file="SolCAP", append=TRUE)
write("\nnphen 0", file="SolCAP", append=TRUE)
write("pheno---------------------------------------", file="SolCAP", append=TRUE)
write("geno---------------------------------------", file="SolCAP", append=TRUE)
write.table(DF, file="SolCAP", append=TRUE, quote=FALSE,
            row.names=FALSE, col.names=FALSE)

