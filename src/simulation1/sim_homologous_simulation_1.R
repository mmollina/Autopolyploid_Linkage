sim_homologous_simulation_1<-function(m, n.mrk, min.d = 0, max.d, prob.dose = NULL,
                                      max.ph, choose.hom.at.random=FALSE,
                                      restriction=TRUE, seed=NULL)
{
  prob.dose<-NULL
  if(!is.null(seed)) set.seed(seed)
  if(choose.hom.at.random)
  {
    hom.allele.q<-hom.allele.p<-vector("list", n.mrk)
    count<-1
    while(count <= n.mrk)
    {
      hom.p.temp<-hom.q.temp<-0
      if(any(is.null(prob.dose)))
        p.temp<-sample(min.d:max.d,1)
      else
        p.temp<-sample(min.d:max.d,1, prob = prob.dose)
      if(all(p.temp!=0))
        hom.p.temp<-sample(1:m, p.temp)
      if(any(is.null(prob.dose)))
        q.temp<-sample(min.d:max.d,1)
      else
        q.temp<-sample(min.d:max.d,1, prob = prob.dose)
      if(all(q.temp!=0))
        hom.q.temp<-sample(1:m, q.temp)
      p <- sum(as.logical(hom.p.temp))
      q <- sum(as.logical(hom.q.temp))
      if(restriction && count > 1)
      {
        if(!any((p+q)==0,
                (p+q)==2*m,
                sum(as.logical(hom.allele.p[[count-1]])) - q == 0,
                sum(as.logical(hom.allele.q[[count-1]])) - p == 0,
                (p==m & q==0),
                (p==0 & q==m)))
        {
          hom.allele.p[[count]] <- hom.p.temp
          hom.allele.q[[count]] <- hom.q.temp
          count<-count+1
        }
      }
      else
      {
        if(!any((p+q)==0,
                (p+q)==2*m,
                (p==m & q==0),
                (p==0 & q==m)))
        {
          hom.allele.p[[count]] <- hom.p.temp
          hom.allele.q[[count]] <- hom.q.temp
          count<-count+1
        }
      }
      p<-unlist(lapply(hom.allele.p, function(x) sum(as.logical(x))))
      q<-unlist(lapply(hom.allele.q, function(x) sum(as.logical(x))))
    }
    return(list(hom.allele.p=hom.allele.p,
                hom.allele.q=hom.allele.q,
                p=p, q=q))
  }
  else{
    p<-q<-rep(0,n.mrk)
    while(any(any((p+q)==0),
              any((p+q)==2*m),
              any((p[-1]+q[-n.mrk])==0),
              any((q[-1]+p[-n.mrk])==0),
              any(p==m & q==0),
              any(p==0 & q==m)))
    {
      hom.allele.q<-hom.allele.p<-vector("list", n.mrk)
      hom.allele.q[]<-hom.allele.p[]<-0
      for(i in 1:n.mrk){
        if(any(is.null(prob.dose)))
          p<-sample(min.d:max.d,1)
        else
          p<-sample(min.d:max.d,1, prob = prob.dose)
        if(all(p!=0)){
          p.add<-p
          if(max.ph < p.add)
            p.add<-max.ph
          ph.p<-sample(0:p.add,1)
          hom.allele.p[[i]]<-c(1:p)+ph.p
        }
        if(any(is.null(prob.dose)))
          q<-sample(min.d:max.d,1)
        else
          q<-sample(min.d:max.d,1, prob = prob.dose)
        if(all(q!=0)){
          q.add<-q
          if(max.ph < q.add)
            q.add<-max.ph
          ph.q<-sample(0:q.add,1)
          hom.allele.q[[i]]<-c(1:q)+ph.q
        }
      }
      p<-unlist(lapply(hom.allele.p, function(x) sum(as.logical(x))))
      q<-unlist(lapply(hom.allele.q, function(x) sum(as.logical(x))))
    }
    return(list(hom.allele.p=hom.allele.p,
                hom.allele.q=hom.allele.q,
                p=p, q=q))
  }
}
