## Auxiliary codes for the Mollinari and Garcia 2019
## ------------------------------------------------------
## These codes are useful to clarify the method
## described in the article.
## For analysis in real situations, use the program MAPpoly,
## available at https://github.com/mmollina/mappoly
## ------------------------------------------------------
## Algorithm 1
##------------
## This is the Algorithm 1 in Mollinari and Garcia 2019
## For instance, in Table 2, b_17 can be obtained using
## bol_vec(6, 17)
## This function is particularly important for higher
## ploidy levels:
## bol_vec(12, 324)
bol_vec<-function(m, i)
{
  s0<-0
  s1<-1
  increment<-0
  sentinel<-0
  b<-numeric(m)
  while(sentinel < m/2)
  {
    temp<-choose(m-s1,(m/2)-s0-1)
    if(i > temp + increment)
    {
      b[s1]<-0
      increment<-increment+temp
    }
    else
    {
      b[s1]<-1
      s0<-s0+1
    }
    sentinel<-sentinel+b[s1]
    s1<-s1+1
  }
  return(b)
}

## This function returns the number of recombinant bivalents
## given two indices, j and j_prime, which indicate the genotypic
## state in positions k and k+1, respectively. For example
## in Table 10, row 17 (j=17), column 14 (j_prime=14), lP and lQ
## can be obtained using
## ll_from_jj(17, 14, 4)
## the result is lP=0 and lQ=2. Thus the transition is (1-r)^2 * r^2
ll_from_jj<-function(j, j_prime, m)
{
    i=1 + (j-1) %/% choose(m,m/2)
    h=choose(m,m/2) + j - i*choose(m,m/2)
    pk<-bol_vec(m,i)
    qk<-bol_vec(m,h)
    ip=1 + (j_prime-1) %/% choose(m,m/2)
    hp=choose(m,m/2) + j_prime - ip*choose(m,m/2)
    pkp<-bol_vec(m,ip)
    qkp<-bol_vec(m,hp)
    return(c(lP=sum(abs(pk - pkp))/2, lQ=sum(abs(qk - qkp))/2))
}

## This function computes zeta, defined in equation 20

## D_k and D_k_prime  are the genotypic states in the reduced dimension
## l_P and l_Q are the number of recombinant bivalents
## phi_* is the linkage phase configuration between
## markers k and k_prime in parents P and Q
## For example, for the following linkage phase configuration
##
## -o----      -o----
## -o----      ----o-
## ----o-  X   ------
## ------      ------
##
##
## phi_1_P = c(1,1,0,0)
## phi_2_P = c(0,0,1,0)
## phi_1_Q = c(1,0,0,0)
## phi_2_Q = c(0,1,0,0)
## m is the ploidy level
##
zeta<-function(D_k,
               D_k_prime,
               l_P,
               l_Q,
               phi_k_P,
               phi_k_prime_P,
               phi_k_Q,
               phi_k_prime_Q,
               m)
{
  v_k<-kronecker(apply(combn(phi_k_P, m/2), MARGIN = 2, sum),
                 apply(combn(phi_k_Q, m/2), MARGIN = 2, sum), FUN = "+")
  v_k_prime<-kronecker(apply(combn(phi_k_prime_P, m/2), MARGIN = 2, sum),
                       apply(combn(phi_k_prime_Q, m/2), MARGIN = 2, sum), FUN = "+")
  ct<-0
  tau_k<-which(v_k==D_k)
  tau_k_prime<-which(v_k_prime==D_k_prime)
  for(j in tau_k)
  {
    for(j_prime in tau_k_prime)
      {
        x<-ll_from_jj(j,j_prime,m)
        if(x["lP"] == l_P && x["lQ"]==l_Q)
          ct<-ct+1
      }
  }
  ct<-ct/choose(m, m/2)^2
  ct
}

##This is to obtain Table 5, Appendix C
m<-4

## Linkage phase
phi_1_P = c(1,1,0,0)
phi_2_P = c(0,0,1,0)
phi_1_Q = c(1,0,0,0)
phi_2_Q = c(0,1,0,0)

## all dosage combinations
D1<-rep(0:m, each=m+1)
D2<-rep(0:m, m+1)

## all possible combinations of l_P and l_Q
l<-expand.grid(0:(m/2), 0:(m/2))[,2:1]

res<-matrix(NA, nrow=length(D1), ncol=nrow(l))
for(i in 1:nrow(l))
{
  for(j in 1:length(D1))
  {
    res[j,i]<-zeta(D_k = D1[j],
                   D_k_prime = D2[j],
                   l_P = l[i,1],
                   l_Q = l[i,2],
                   phi_k_P = phi_1_P,
                   phi_k_prime_P = phi_2_P,
                   phi_k_Q = phi_1_Q,
                   phi_k_prime_Q = phi_2_Q,
                   m)
  }
}
colnames(res)<-apply(l, 1, paste, collapse="-")
rownames(res)<-paste(D1, D2, sep="-")
res
require(MASS)
fractions(res)
